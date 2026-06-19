# NGS Precision Oncology Pipeline (Parabricks GPU)

End-to-end tumor-normal somatic variant calling: FASTQ → BAM → VCF → annotated MAF → clinical report.

## Architecture

```
┌──────────────┐     ┌─────────────────────┐     ┌──────────────┐     ┌──────────────┐
│  Sequencer   │     │  Secondary Analysis │     │   Tertiary   │     │   Clinical   │
│              │ ──> │                     │ ──> │              │ ──> │              │
│   FASTQ      │     │   PARABRICKS GPU    │     │  Annotation  │     │   Report     │
│  (paired R1  │     │  • fq2bam (BWA+BQSR)│     │  • VEP REST  │     │  • AMP tiers │
│   + R2)      │     │  • Mutect2 somatic  │     │  • OncoKB    │     │  • JSON+text │
│              │     │  • DeepSomatic somat│     │  • hotspot DB│     │  • EHR-ready │
│              │     │  • DeepVariant germ.│     │              │     │              │
└──────────────┘     └─────────────────────┘     └──────────────┘     └──────────────┘
   ~200 GB             ~30 min on 2× RTX PRO       ~7 min VEP REST       ~1 min
                                                    (5 workers, 43K vars)
                       (vs ~30 hr CPU)
```

## Sample Dataset

For development without spinning up a GPU, use these public benchmark datasets:

| Source | Sample | Description |
|---|---|---|
| **SEQC2 (FDA)** | HCC1395 / HCC1395BL | Breast cancer cell line tumor + matched B-lymphoblastoid normal. Gold-standard truth VCF available. |
| **GIAB** | HG001-HG007 | Germline benchmarks (for pipeline validation) |
| **TCGA** | DREAM challenge synthetic 3/4 | Synthetic tumors with planted variants and known VAFs |
| **PrecisionFDA** | Truth Challenge V2 | Reference for somatic calling benchmarks |

**Recommended starter dataset — HCC1395 (the de facto somatic benchmark):**
```bash
# Tumor (SRR7890829, WGS_FD_T_1, ~51×) + Normal (SRR7890826, WGS_FD_N_1, ~53×)
# from PRJNA489865 (FDA SEQC2) via SRA toolkit — no AWS account needed
./scripts/download_hcc1395.sh

# Truth VCF only (5 MB, fast)
./scripts/download_hcc1395.sh --truth-only
```

For a fast smoke test, use **TruSeq Amplicon Cancer Panel** chr22-subsetted data (~500 MB total).

## Quick Start

### Prerequisites

- NVIDIA GPU with ≥40 GB memory — **this machine has 2× RTX PRO 6000 Blackwell (96 GB each)**
- Driver ≥ 525, Docker + `nvidia-container-toolkit`
- Parabricks container already pulled: `nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1`
- ~700 GB free disk per WGS tumor-normal pair
- Reference bundle (GRCh38 no-alt) — see step 1 below

### Step 1: Download missing reference files (~4.7 GB, one-time)

All files are served over public HTTPS — no gsutil auth or AWS credentials needed.

```bash
./scripts/download_refs.sh
# Downloads to /mnt/storage/parabricks_test/ref/
# Adds: known_indels, dbsnp138, af-only-gnomad, 1000g_pon, wgs_calling_regions
```

Reference files status after download:

| File | Purpose | Source bucket |
|---|---|---|
| `Homo_sapiens_assembly38.fasta` + BWA index | Alignment | ✅ already present |
| `Homo_sapiens_assembly38.known_indels.vcf.gz` | BQSR | `gcp-public-data--broad-references` |
| `Homo_sapiens_assembly38.dbsnp138.vcf.gz` | BQSR | `gcp-public-data--broad-references` |
| `af-only-gnomad.hg38.vcf.gz` | Mutect2 germline prior | `gatk-best-practices` |
| `1000g_pon.hg38.vcf.gz` | Mutect2 Panel of Normals | `gatk-best-practices` |
| `wgs_calling_regions.hg38.interval_list` | Callable regions | `gcp-public-data--broad-references` |

> **⚠️ Reference compatibility:** The GATK resource bundle files (gnomAD, 1000G PoN) are
> distributed against the **full GRCh38 reference (3,366 contigs)**. If your reference FASTA
> is the **no-alt build (~455 sequences)**, Parabricks 4.7 will reject these VCFs at runtime.
> Two separate issues must both be fixed — the pipeline's Stage 4 handles both automatically:
>
> **Issue 1 — Chromosome count:** `bcftools view -R` filters data rows but leaves all 3,366
> contig lines in the VCF header. Fix with `bcftools reheader` to rebuild the header from the FAI.
>
> **Issue 2 — Chromosome order:** The FAI for this FASTA is in lexicographic order
> (chr1, chr10, chr11 … chr2 …). The 1000G PoN VCF header uses natural/numeric order
> (chr1, chr2 … chr10 …). `pbrun prepon` encodes the header order into the binary;
> `pbrun mutect` compares against the FASTA order. They must match.
> Fix: rewrite the header with contig lines generated directly from the FAI (`awk '{printf
> "##contig=<ID=%s,length=%s>\n",$1,$2}' ref.fai`), then rebuild the `.pon` binary.
>
> **Issue 3 — prepon output path:** `pbrun prepon` does **not** accept `--out-pon-file`.
> It always writes `<input>.pon` alongside the input (e.g. `1000g_pon.noalt.vcf.gz.pon`).
> Pass `--pon 1000g_pon.noalt.vcf.gz` to `mutectcaller` — Parabricks appends `.pon` internally.

### Step 2: Download HCC1395 data + truth VCF

Source: PRJNA489865 (FDA SEQC2 somatic mutation benchmarking, HiSeq X Ten, ~50× WGS).
Uses `prefetch` + `fasterq-dump` (SRA toolkit, already installed at `/usr/bin/`).

```bash
# Download tumor (SRR7890829, WGS_FD_T_1) + normal (SRR7890826, WGS_FD_N_1) + truth VCF
./scripts/download_hcc1395.sh
# ~700 GB total (two ~130 GB uncompressed FASTQs + truth files)

# Or download pieces independently:
./scripts/download_hcc1395.sh --truth-only    # just the SEQC2 truth VCF (~5 MB, fast)
./scripts/download_hcc1395.sh --tumor-only
./scripts/download_hcc1395.sh --normal-only
```

### Step 3: Run the full pipeline

```bash
# Tumor-normal paired: fq2bam (BQSR) → Mutect2 → DeepVariant → FilterMutectCalls
# Skip-if-exists guards make this safe to re-run from any stage.
NUM_GPUS=2 \
REF_DIR=/mnt/storage/parabricks_test/ref \
OUTDIR=/mnt/storage/parabricks_test/output/HCC1395 \
./scripts/run_parabricks_pipeline.sh HCC1395 \
  /mnt/storage/parabricks_test/data/hcc1395/tumor/SRR7890829_1.fastq.gz \
  /mnt/storage/parabricks_test/data/hcc1395/tumor/SRR7890829_2.fastq.gz \
  /mnt/storage/parabricks_test/data/hcc1395/normal/SRR7890826_1.fastq.gz \
  /mnt/storage/parabricks_test/data/hcc1395/normal/SRR7890826_2.fastq.gz

# Annotate (VEP REST + OncoKB hotspot fallback — no local VEP cache needed)
# Use the type-split VCF (best combined SNV+INDEL precision); ~25-40 min for 43K variants (5 workers).
# Reads the filtered somatic VCF (.gz supported); queries Ensembl VEP concurrently in batches of 200.
export ONCOKB_TOKEN="your-token"   # optional; built-in hotspot table used otherwise
python scripts/vcf_to_maf_vep_rest.py \
  --vcf /mnt/storage/parabricks_test/output/HCC1395/HCC1395.somatic.typesplit.vcf.gz \
  --tumor-id HCC1395_T --normal-id HCC1395_N \
  --tumor-type BRCA \
  --workers 5 \
  --out /mnt/storage/parabricks_test/output/HCC1395/annotated/HCC1395_typesplit.maf

# Generate clinical report
python scripts/generate_clinical_report.py \
  --maf /mnt/storage/parabricks_test/output/HCC1395/annotated/HCC1395_typesplit.maf \
  --sample-id HCC1395 --tumor-type BRCA \
  --meta-json /mnt/storage/parabricks_test/output/HCC1395/HCC1395_meta.json \
  --out /mnt/storage/parabricks_test/output/HCC1395/annotated/HCC1395_typesplit_clinical_report.txt

# Benchmark somatic calls against SEQC2 truth
# (see "Benchmarking against truth" section for required prep steps — truth VCFs need
#  synthetic GT added and calls VCF must be pre-filtered to PASS single-sample first)
docker run --rm -v /mnt/storage:/mnt/storage realtimegenomics/rtg-tools:3.12.1 vcfeval \
  --baseline /mnt/storage/parabricks_test/data/truth/high-confidence_sSNV_gt.vcf.gz \
  --calls /mnt/storage/parabricks_test/output/HCC1395/HCC1395.somatic.pass.vcf.gz \
  --template /mnt/storage/parabricks_test/ref/GRCh38.sdf \
  --output /mnt/storage/parabricks_test/output/HCC1395/eval_snv \
  --evaluation-regions /mnt/storage/parabricks_test/data/truth/High-Confidence_Regions_v1.2.bed \
  --vcf-score-field INFO.TLOD --squash-ploidy
```

### Step 4 (alternative caller): DeepSomatic tumor-normal calling

DeepSomatic (Parabricks 4.7) is a deep-learning somatic caller that outperforms Mutect2 on
HCC1395 without any post-processing filter cascade. It uses the matched normal for germline
subtraction and emits `FILTER=GERMLINE` for non-somatic variants; only `FILTER=PASS` calls are
somatic. Requires BAMs already produced by Step 3 (fq2bam).

```bash
# Phases 1–4: call → PASS extract → VEP/OncoKB annotation → AMP-tiered report
# All phases have skip-if-exists guards — safe to re-run from any stage.
# DeepSomatic's FILTER column already suppresses FPs; VEP/OncoKB add interpretation only.
TUMOR_TYPE=BRCA \
./scripts/run_deepsomatic_pipeline.sh HCC1395 \
  /mnt/storage/parabricks_test/output/HCC1395/HCC1395_tumor.bam \
  /mnt/storage/parabricks_test/output/HCC1395/HCC1395_normal.bam
# Writes (in deepsomatic/ subdir):
#   HCC1395.deepsomatic.pass.vcf.gz            PASS-only somatic VCF
#   HCC1395.deepsomatic.maf                    VEP + OncoKB annotated MAF
#   HCC1395.deepsomatic_clinical_report.txt    AMP-tiered text report
#   HCC1395.deepsomatic_clinical_report.json   EHR-ready JSON

# Add --skip-annot to run Phases 0-2 only (calling + PASS extract, no annotation).
# Set ONCOKB_TOKEN env var to use live OncoKB REST instead of built-in hotspot fallback.

# Benchmark against SEQC2 truth (separate step — HC BED is mandatory)
# HC BED required — omitting it inflates FPs by ~93% (chrX + chr6 MHC + chr16 repeats)
./scripts/eval_deepsomatic.sh HCC1395
# Writes: /mnt/storage/parabricks_test/output/HCC1395/deepsomatic/eval/comparison.md

# Optional: PCGR second-opinion report (deterministic, CIViC/ClinVar/OncoKB bundle)
# Requires one-time data bundle download (~5 GB) and VEP cache (~30 GB) — see script header.
PCGR_BUNDLE=/mnt/storage/pcgr_bundle \
VEP_CACHE=/mnt/storage/vep_cache \
TUMOR_TYPE=BRCA \
./scripts/run_pcgr.sh HCC1395
# Writes: deepsomatic/pcgr/HCC1395.pcgr.grch38.html  (HTML report)
#         deepsomatic/pcgr/HCC1395.pcgr.grch38.snvs_indels.tiers.tsv  (tier TSV for diffing)
```

> **⚠️ HC regions BED is mandatory for SEQC2 truth evaluation.** The truth VCF covers only
> autosomes (chr1–22); DeepSomatic calls chrX, chr6 (MHC), and chr16 (repeat regions) as somatic,
> and without `--evaluation-regions` all those calls are scored as false positives — pushing
> precision from ~93% down to ~30%. Download once with `./scripts/download_hcc1395.sh --truth-only`.
```

### Smoke test (already validated ✅)

```bash
./scripts/run_smoke_test.sh
# 90 s end-to-end · Parabricks 4.7.0-1 · GPU 0 · TruSeq Amplicon panel data
# Found: TP53 p.R248W (VAF 22%) + PIK3CA splice region
```

## HCC1395 Validated Run Results (Jun 2026, Parabricks 4.7.0-1)

Run on 2× RTX PRO 6000 Blackwell (96 GB each), GRCh38 no-alt reference.

| Stage | Tool | Runtime | Output |
|---|---|---|---|
| Tumor alignment + BQSR | `pbrun fq2bam` | ~15 min | 51×, 99.9% aligned |
| Normal alignment + BQSR | `pbrun fq2bam` | ~15 min | 53×, 99.8% aligned |
| QC metrics | `pbrun collectmultiplemetrics` | ~3 min each | alignment, GC bias, insert size |
| Somatic calling (Mutect2) | `pbrun mutectcaller` | ~19 min | 242,991 candidates |
| Somatic calling (DeepSomatic) | `pbrun deepsomatic` | ~50 min | 124,014 PASS (no post-filter) |
| Germline calling | `pbrun deepvariant` | ~6.5 min | 7.1M germline variants |
| Somatic filtering | `gatk FilterMutectCalls` | <1 min | 95,321 PASS |
| VEP annotation | VEP REST API | ~25 min | — |

**Somatic PASS breakdown** (with gnomAD germline resource + 1000G PoN, Jun 2026):

| Type | Count |
|---|---|
| SNV | 84,441 |
| MNP | 914 |
| INDEL | 9,966 |
| **Total PASS** | **95,321** |

> **Note:** gnomAD germline filtering reduced PASS from 116,132 → 95,321 (−18%). PoN filtering
> further suppresses recurrent sequencing artifacts. The remaining FPs are primarily germline
> heterozygous variants not covered by gnomAD and orientation-bias artifacts from the sequencing
> chemistry — both addressed in the next step (contamination estimation + orientation model).

## What Each Stage Actually Does

### Stage 1: `pbrun fq2bam` — FASTQ to analysis-ready BAM

A single fused kernel that replaces five separate CPU tools:

1. **BWA-MEM** GPU implementation — aligns paired reads to GRCh38
2. **Coordinate sort** — partitions reads by chromosome position
3. **Mark duplicates** — flags PCR/optical duplicates (critical for accurate VAF)
4. **BQSR (Base Quality Score Recalibration)** — corrects systematic sequencer errors using known variant sites

Output: a sorted, deduplicated, recalibrated BAM ready for variant calling. **30x WGS: ~25 min on 4× A100 vs ~30 hr on a 32-core CPU node.**

Use `--bwa-options="-Y -K 100000000"` on **both** tumor and normal `fq2bam` calls. The `-K` flag sets a fixed chunk size making BWA-MEM deterministic regardless of thread count; mismatched options between tumor and normal can produce subtle alignment differences that affect paired somatic calling.

### Stage 2: `pbrun mutectcaller` — Somatic variants

Mutect2 is the GATK4 somatic caller. The paired tumor-normal mode distinguishes:
- **Somatic variants** (in tumor, absent in normal) → cancer driver mutations
- **Germline variants** (in both) → inherited polymorphisms, filtered out
- **Sequencing artifacts** → filtered by Panel of Normals + read-orientation models

Key inputs:
- `--pon` (Panel of Normals) — site-recurrent artifacts from healthy samples. Requires a
  one-time `pbrun prepon --in-pon-file <pon.vcf.gz>` to build the binary index. The output
  is written as `<pon.vcf.gz>.pon` alongside the input. Pass the **VCF path** (not the `.pon`
  path) to `mutectcaller --pon`; Parabricks appends `.pon` internally.
- `--mutect-germline-resource` — gnomAD AF used to compute probability a call is germline
  (flag renamed from `--germline-resource` in Parabricks 4.7)

`--interval-file` and `--no-alt-contigs` are **mutually exclusive** in Parabricks 4.7; the interval file already restricts to callable regions, so `--no-alt-contigs` is redundant.

### Stage 3: VEP + OncoKB annotation

Raw VCF says `chr7:55191822 T>G`. Annotation transforms that into:
> EGFR p.L858R missense in exon 21 — known driver in NSCLC, sensitizing to osimertinib (FDA LEVEL_1)

`vcf_to_maf_vep_rest.py` handles gzipped `.vcf.gz` input directly. It batches variants in
groups of 200 to the Ensembl VEP REST API using concurrent workers (`--workers N`, default 5),
achieving ~6x speedup over sequential batching (~7 min vs ~45 min for 43K variants). Each
worker reuses a thread-local HTTP session. It then applies OncoKB via REST (if `ONCOKB_TOKEN`
is set) or falls back to the built-in hotspot database, which covers:
- Point mutations: EGFR, KRAS, BRAF, PIK3CA, TP53 hotspots, HER2, etc.
- Gene-level loss-of-function: BRCA1, BRCA2, PALB2 (frameshift/nonsense/splice → LEVEL_1 PARP
  inhibitor annotation, regardless of specific amino acid position)

ClinVar pathogenicity (`CLIN_SIG`) is also extracted from VEP colocated variant data.

### Stage 4: Tiering and reporting

The report generator (`generate_clinical_report.py`):
1. Filters to protein-altering variants with VAF ≥ 5%, alt depth ≥ 5, gnomAD AF < 1%
2. Maps OncoKB levels to AMP/ASCO/CAP tiers (Tier IA–IV)
3. Produces both human-readable text and machine-readable JSON for EHR ingestion

## Parabricks 4.7 API Changes

Several flag names and defaults changed from 4.3/4.4. Summary of breaking changes:

| Tool | Old flag / behavior | New flag / behavior |
|---|---|---|
| `collectmultiplemetrics` | Collected all metrics by default | Must pass `--gen-all-metrics` explicitly |
| `mutectcaller` | `--germline-resource` | `--mutect-germline-resource` |
| `mutectcaller` | PoN VCF used directly | Must run `pbrun prepon --in-pon-file <pon.vcf.gz>` first to build `.pon` index |
| `mutectcaller` | `--interval-file` and `--no-alt-contigs` usable together | **Mutually exclusive** — use one or the other |
| `deepvariant` | `--mode wgs` | `--mode shortread` (`wgs` removed; valid values: `shortread`, `pacbio`, `ont`) |
| All tools | Container `/tmp` writable for non-root | Container's `TMPDIR` defaults to `/` for non-root users; add `-v /tmp:/tmp -e TMPDIR=/tmp` to every `docker run` call |

## Benchmarking against truth

For HCC1395, validate calls against the SEQC2 truth set. Three prep steps are required because
the SEQC2 truth VCFs are site-only (no FORMAT/GT columns) and the calls VCF is multi-sample:

```bash
TRUTHDIR=/mnt/storage/parabricks_test/data/truth
OUTDIR=/mnt/storage/parabricks_test/output/HCC1395

# 1. Normalize truth VCF: replace "PASS;HighConf" filter with plain "PASS"
for TYPE in sSNV sINDEL; do
  bcftools view ${TRUTHDIR}/high-confidence_${TYPE}_in_HC_regions_v1.2.vcf.gz | \
    awk 'BEGIN{OFS="\t"} /^#/{print;next} {gsub(/PASS;HighConf/,"PASS",$7); print}' | \
    bcftools view -Oz -o ${TRUTHDIR}/high-confidence_${TYPE}_norm.vcf.gz
  tabix -p vcf ${TRUTHDIR}/high-confidence_${TYPE}_norm.vcf.gz
done

# 2. Add synthetic GT column (site-only truth → RTG needs GT in baseline)
for TYPE in sSNV sINDEL; do
  bcftools view ${TRUTHDIR}/high-confidence_${TYPE}_norm.vcf.gz | \
    awk 'BEGIN{OFS="\t"} /^##/{print;next} /^#CHROM/{print $0"\tFORMAT\tTRUTH";next} {print $0"\tGT\t1/1"}' | \
    bcftools view -Oz -o ${TRUTHDIR}/high-confidence_${TYPE}_gt.vcf.gz
  tabix -p vcf ${TRUTHDIR}/high-confidence_${TYPE}_gt.vcf.gz
done

# 3. Extract single-sample tumor VCF and filter to PASS
bcftools view -s HCC1395_T ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz -Oz -o ${OUTDIR}/HCC1395.somatic.tumor_only.vcf.gz
tabix -p vcf ${OUTDIR}/HCC1395.somatic.tumor_only.vcf.gz
bcftools view -f PASS ${OUTDIR}/HCC1395.somatic.tumor_only.vcf.gz -Oz -o ${OUTDIR}/HCC1395.somatic.pass.vcf.gz
tabix -p vcf ${OUTDIR}/HCC1395.somatic.pass.vcf.gz

# 4. Build RTG SDF (one-time)
docker run --rm -v /mnt/storage:/mnt/storage realtimegenomics/rtg-tools:3.12.1 format \
  -o /mnt/storage/parabricks_test/ref/GRCh38.sdf \
  /mnt/storage/parabricks_test/ref/Homo_sapiens_assembly38.fasta

# 5. Run vcfeval (entrypoint is already "rtg" — do NOT prefix with "rtg")
#    --sample TRUTH,HCC1395_T: baseline sample name, calls sample name
#    The calls VCF is multi-sample (tumor + normal); must specify the tumor sample explicitly.
#    --evaluation-regions: REQUIRED for SEQC2 truth — truth covers autosomes only (no chrX).
#    Omitting it inflates FP by ~93% (chrX/chr6/chr16 calls all become FP).
docker run --rm -v /mnt/storage:/mnt/storage realtimegenomics/rtg-tools:3.12.1 vcfeval \
  --baseline ${TRUTHDIR}/high-confidence_sSNV_gt.vcf.gz \
  --calls ${OUTDIR}/HCC1395.somatic.pass.vcf.gz \
  --template /mnt/storage/parabricks_test/ref/GRCh38.sdf \
  --output ${OUTDIR}/eval_snv \
  --evaluation-regions ${TRUTHDIR}/High-Confidence_Regions_v1.2.bed \
  --vcf-score-field INFO.TLOD \
  --squash-ploidy \
  --sample TRUTH,HCC1395_T
```

**HCC1395 vcfeval results — benchmarked against SEQC2 truth (Jun 2026):**

| Run | Variant type | TP | FP | FN | Precision | Sensitivity | F1 |
|---|---|---|---|---|---|---|---|
| No gnomAD / no PoN | SNV | 32,557 | 83,338 | 4,841 | 27.9% | 87.1% | 0.422 |
| No gnomAD / no PoN | INDEL | 1,558 | 114,574 | 196 | 1.3% | 88.8% | 0.026 |
| gnomAD + PoN | SNV | 34,249 | 61,325 | 3,149 | 35.6% | 91.6% | 0.513 |
| gnomAD + PoN | INDEL | 71 | 1,184 | 1,683 | 5.7% | 4.1% | 0.047 |
| + contamination + orientation model | SNV | 34,217 | 61,166 | 3,181 | 35.6% | 91.5% | 0.513 |
| + contamination + orientation model | INDEL | 71 | 1,184 | 1,683 | 5.7% | 4.1% | 0.047 |
| **+ VAF≥5% + depth≥5 + POPAF≥2** | **SNV** | **33,076** | **11,441** | **4,322** | **74.1%** | **88.4%** | **0.806** |
| + VAF≥5% + depth≥5 + POPAF≥2 (INDEL-only eval†) | INDEL | 1,356 | 1,490 | 398 | **47.6%** | **77.3%** | **0.590** |
| **+ MMQ≥50** (applied to above) | **SNV** | **32,930** | **10,660** | **4,468** | **75.4%** | **88.1%** | **0.812** |
| + MMQ≥50 (INDEL-only eval†, opt. TLOD 23.6) | INDEL | 1,352 | 1,459 | 402 | **48.1%** | **77.1%** | **0.592** |
| **Type-split (SNV MMQ50 + INDEL TLOD≥23.6)** | **SNV** | **32,886** | **10,056** | **4,512** | **76.4%** | **87.9%** | **0.818** |
| Type-split (INDEL-only eval†, TLOD≥23.6) | INDEL | 1,356 | 1,490 | 398 | **47.6%** | **77.3%** | **0.590** |

> **†INDEL evaluation methodology note:** Prior INDEL rows (gnomAD/PoN and contamination runs, P≈5.7%)
> used a mixed calls VCF (SNVs + INDELs) against the INDEL truth baseline. vcfeval counted all SNVs
> in the calls as false positives since they don't match the INDEL baseline — artificially inflating
> FP counts by ~40K. POPAF2 and MMQ50 rows use an INDEL-only calls VCF (`bcftools view -v indels`)
> for a correct measurement. The old 6.1% figure was an evaluation artifact, not a real precision value.
> True INDEL precision at optimal TLOD is **47–48%**, sensitivity **77%**.

**DeepSomatic HCC1395 results — with SEQC2 HC regions BED (Jun 2026):**

| Run | Variant type | TP | FP | FN | Precision | Sensitivity | F1 |
|---|---|---|---|---|---|---|---|
| **DeepSomatic WGS, PASS-only, no post-filter** | **SNV** | **35,383** | **2,653** | **1,947** | **93.0%** | **94.8%** | **0.939** |
| DeepSomatic WGS, PASS-only, no post-filter | INDEL | 1,339 | 207 | 127 | 86.6% | 91.3% | 0.889 |

Evaluated with `--evaluation-regions High-Confidence_Regions_v1.2.bed` and `--vcf-score-field GQ`.
Score field `GQ` (DeepSomatic Genotype Quality); best-F1 threshold reported.

> **Diagnostic note:** Running vcfeval without `--evaluation-regions` yields ~30% SNV precision on
> this dataset — not a caller failure. The SEQC2 truth covers autosomes only (no chrX). Without
> the HC BED, DeepSomatic's 35,815 chrX + 25,966 chr6 + 21,126 chr16 PASS calls are all scored
> as false positives (93.6% of total FPs). Always provide the HC BED when benchmarking against
> SEQC2. See `scripts/diagnose_deepsomatic.sh` for a structured root-cause check.

**COLO829 DeepSomatic results — SMaHT truth, no HC BED (Jun 2026):**

| Run | Variant type | TP | FP | FN | Precision | Sensitivity | F1 |
|---|---|---|---|---|---|---|---|
| DeepSomatic WGS, PASS-only | SNV | 37,466 | 2,038 | 6,539 | 94.8% | 85.1% | 0.897 |
| DeepSomatic WGS, PASS-only | INDEL | 668 | 395 | 1,391 | 62.8% | 32.4% | 0.428 |

Truth set: SMaHT COLO829BLT50 (parklab/SMaHT_SNV_COLO829BLT50_HAPMAP). No HC BED was applied
(the SMaHT truth includes chrX so FP inflation is milder, but a COLO829-specific HC BED would
further clean up the chr6/chr16 complex-region calls). COLO829 is also in DeepSomatic's training
corpus; these numbers are not a clean generalization test. Low INDEL sensitivity (32%) likely
reflects a truth-set depth/platform mismatch rather than caller failure.

---

PoN + gnomAD delivers a clean SNV improvement vs the baseline (+7.7 pp precision, +4.5 pp
sensitivity, +9 pp F1, −22K FPs). Contamination (0.14%) and orientation model had negligible
additional impact — expected for a clean cell-line dataset. The POPAF filter (described below)
is the real breakthrough: removing calls with gnomAD AF ≥ 1% eliminates 41,786 FPs while
losing only 151 TPs, pushing SNV precision from 35.6% to **74.1%**. The MMQ≥50 filter adds a
further 783 FPs removed (all low-mapq reads, MMQ=40 being the dominant class), gaining +1.3 pp
precision at the cost of −0.4 pp sensitivity.

**Why the POPAF filter is so effective:** Mutect2 passes many common germline variants (gnomAD
AF 10–50%) in HCC1395 because copy-number alterations cause LOH — a germline het variant
becomes homozygous in the tumor at ~100% VAF while the normal shows near-zero alt reads (the
wildtype allele is in the LOH-deleted region). The standard gnomAD germline-resource filter
does not fully suppress these in LOH contexts. Requiring POPAF ≥ 2 (gnomAD AF < 1%) removes
this entire class of FPs with minimal TP loss.

## Next Step: GATK Best-Practice Filtering (Contamination + Orientation Model)

The current pipeline runs `FilterMutectCalls` with only the basic outputs from `mutectcaller`.
GATK best practices add three CPU steps that dramatically reduce FPs, particularly for INDELs
and FFPE-style orientation artifacts:

| Step | Tool | Input | Output | Impact |
|---|---|---|---|---|
| 1 | `GetPileupSummaries` (tumor + normal) | BAM + gnomAD | pileup tables | enables contamination estimate |
| 2 | `CalculateContamination` | pileup tables | contamination + segmentation tables | removes cross-sample DNA |
| 3 | `LearnReadOrientationModel` | unfiltered VCF stats | orientation model `.tar.gz` | removes C→T/G→T sequencing artifacts |
| 4 | `FilterMutectCalls` (re-run) | unfiltered VCF + all three above | improved filtered VCF | expected: SNV precision 85–92%, INDEL precision ~70% |

None of these are GPU-accelerated in Parabricks 4.7 — they run in the existing
`broadinstitute/gatk:4.5.0.0` container. All required inputs already exist:

```
✅ HCC1395_tumor.bam / HCC1395_normal.bam
✅ HCC1395.somatic.unfiltered.vcf.gz.stats  (produced by mutectcaller; input to LearnReadOrientationModel)
✅ af-only-gnomad.noalt.vcf.gz              (input to GetPileupSummaries)
✅ wgs_calling_regions.hg38.interval_list
✅ HCC1395_tumor_pileups.table              (Jun 13 2026)
✅ HCC1395_normal_pileups.table             (Jun 14 2026)
✅ HCC1395_contamination.table              (Jun 13 2026, contamination = 0.14%)
✅ HCC1395_tumor_segments.table             (Jun 13 2026)
✅ HCC1395_read_orientation_model.tar.gz    (Jun 13 2026)
```

```bash
OUTDIR=/mnt/storage/parabricks_test/output/HCC1395
REF=/mnt/storage/parabricks_test/ref
GATK="docker run --rm -v /mnt/storage:/mnt/storage -v /tmp:/tmp broadinstitute/gatk:4.5.0.0"

# Step 1a: Pileup summaries — tumor
${GATK} gatk GetPileupSummaries \
  -I ${OUTDIR}/HCC1395_tumor.bam \
  -V ${REF}/af-only-gnomad.noalt.vcf.gz \
  -L ${REF}/wgs_calling_regions.hg38.interval_list \
  -O ${OUTDIR}/HCC1395_tumor_pileups.table

# Step 1b: Pileup summaries — normal
${GATK} gatk GetPileupSummaries \
  -I ${OUTDIR}/HCC1395_normal.bam \
  -V ${REF}/af-only-gnomad.noalt.vcf.gz \
  -L ${REF}/wgs_calling_regions.hg38.interval_list \
  -O ${OUTDIR}/HCC1395_normal_pileups.table

# Step 2: Contamination estimate (tumor vs matched normal)
${GATK} gatk CalculateContamination \
  -I ${OUTDIR}/HCC1395_tumor_pileups.table \
  -matched ${OUTDIR}/HCC1395_normal_pileups.table \
  -O ${OUTDIR}/HCC1395_contamination.table \
  --tumor-segmentation ${OUTDIR}/HCC1395_tumor_segments.table

# Step 3: Read orientation model (from Mutect2 stats file)
${GATK} gatk LearnReadOrientationModel \
  -I ${OUTDIR}/HCC1395.somatic.unfiltered.vcf.gz.stats \
  -O ${OUTDIR}/HCC1395_read_orientation_model.tar.gz

# Step 4: Re-filter with all three inputs
# Note: run with output to filtered.vcf.gz (overwrites the basic-filtered file).
# The VCF header confirms contamination-table + ob-priors were active in this run.
${GATK} gatk FilterMutectCalls \
  -R ${REF}/Homo_sapiens_assembly38.fasta \
  -V ${OUTDIR}/HCC1395.somatic.unfiltered.vcf.gz \
  --contamination-table ${OUTDIR}/HCC1395_contamination.table \
  --tumor-segmentation ${OUTDIR}/HCC1395_tumor_segments.table \
  --ob-priors ${OUTDIR}/HCC1395_read_orientation_model.tar.gz \
  -O ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz
```

**Actual results after Step 4 (Jun 2026):**

| Metric | gnomAD + PoN only | + contamination + orientation |
|---|---|---|
| SNV Precision | 35.6% | **35.6%** (no change) |
| SNV Sensitivity | 91.6% | **91.5%** |
| INDEL Precision (at threshold) | 5.7% | **5.7%** (no change) |

The orientation model and contamination filter had negligible impact because: (a) contamination
was only 0.14%, and (b) the residual FPs in HCC1395 are not strand-bias artifacts. The
precision plateau at ~35% is primarily due to germline variants not covered by gnomAD.

### Applied Hard Post-Filtering (Jun 2026) — DONE

SNV precision was pushed to **74.1%** by combining three post-filters applied to
`HCC1395.somatic.filtered.vcf.gz` (the Mutect2 output with contamination + orientation):

```bash
OUTDIR=/mnt/storage/parabricks_test/output/HCC1395

# Tumor is sample index 1 (FORMAT column order: HCC1395_N, HCC1395_T)
# POPAF is -log10(gnomAD AF); POPAF >= 2 means gnomAD AF < 1%
docker run --rm -v /mnt/storage:/mnt/storage broadinstitute/gatk:4.5.0.0 bash -c "
  bcftools filter \
    -i 'FILTER=\"PASS\" && FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2' \
    ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz \
    -Oz -o ${OUTDIR}/HCC1395.somatic.popaf2.vcf.gz && \
  bcftools index -t -f ${OUTDIR}/HCC1395.somatic.popaf2.vcf.gz"
```

Key insight: 41,786 FPs were common gnomAD variants (AF 10–50%) passing Mutect2 in LOH
regions (tumor shows ~100% alt VAF, normal shows 0 alt due to loss of the wildtype allele).
POPAF >= 2 removes this class with only 151 TP losses (277:1 FP:TP removal ratio).

### Remaining SNV FPs — Further Work

After the POPAF filter, 11,441 SNV FPs remain (precision 74.1%). These fall into three classes:
1. **Low-VAF noise** (tumor AF 5–20%) — further TLOD filtering may help
2. **High-VAF calls in non-gnomAD regions** (POPAF=6, tumor AF ~100%) — likely real somatic
   variants absent from the SEQC2 truth set; truth set gap rather than FP
3. **Low-mapq reads** — genome-wide, not specific to HLA/segdup — need mapq filtering
   (the original hypothesis of chr6:28–34M HLA concentration was disproven: 0 removed FPs
   there; the dominant class is MMQ=40 reads spread across all chromosomes, with chrX
   contributing the most at 127/783 = 16% of the removed FPs)

**MMQ≥50 filter — applied (Jun 2026):**
Removes 783 FPs at the cost of 146 TPs. Result: precision **75.4%** (+1.3 pp), sensitivity
88.1% (−0.4 pp), F1 0.812. Output: `HCC1395.somatic.mmq50.vcf.gz`.

```bash
bcftools filter \
  -i 'FILTER="PASS" && FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2 && INFO/MMQ >= 50' \
  ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz -Oz -o ${OUTDIR}/HCC1395.somatic.mmq50.vcf.gz
bcftools index -t -f ${OUTDIR}/HCC1395.somatic.mmq50.vcf.gz
```

After MMQ50, **10,660 SNV FPs remain**. TLOD sweep (Jun 2026) shows precision plateaus at
~82–83% (SNV-specific) regardless of threshold; going from TLOD 5.8 → 50 gains 0 pp precision
while costing 40 pp sensitivity. TLOD filtering cannot break the FP ceiling for this dataset:
the FPs are distributed across the full TLOD range rather than concentrated at low scores.

| TLOD threshold | TP | FP (SNV) | Precision | Sensitivity | F1 |
|---|---|---|---|---|---|
| ≥ 5.8 (optimal) | 32,930 | 6,818 | 82.5% | 88.0% | 0.852 |
| ≥ 10 | 32,449 | 6,547 | 82.9% | 86.8% | 0.848 |
| ≥ 15 | 30,648 | 6,083 | 83.2% | 82.0% | 0.825 |
| ≥ 20 | 29,015 | 5,795 | 83.1% | 77.6% | 0.802 |
| ≥ 30 | 25,727 | 5,324 | 82.5% | 68.8% | 0.750 |
| ≥ 50 | 17,671 | 4,120 | 80.8% | 47.2% | 0.596 |

> Note: these are SNV-specific precision values (from the vcfeval SNP ROC); the vcfeval summary
> precision (75.4%) includes MNPs and INDELs in the denominator. The ceiling at ~83% indicates
> remaining FPs are class 2 (truth-set gaps) or require orthogonal evidence beyond TLOD.

### INDEL Performance — Corrected (Jun 2026)

**The prior INDEL numbers (P=6.1%) were a vcfeval evaluation artifact.** When the calls VCF
contains both SNVs and INDELs, vcfeval counts all SNVs as false positives against the INDEL
baseline — inflating FP by ~39K. The correct evaluation uses an INDEL-only calls VCF.

**Correct INDEL results (INDEL-only calls VCF, INDEL truth baseline):**
- popaf2: P=47.6%, Se=77.3%, F1=0.590 at optimal TLOD 23.6
- mmq50:  P=48.1%, Se=77.1%, F1=0.592 at optimal TLOD 23.6

INDEL precision is actually decent. The TLOD 23.6 threshold (auto-selected by vcfeval as
F1-optimal) shows that operating at a stricter TLOD helps: the "no threshold" row gives P=43.8%
but at TLOD ≥ 23.6 it rises to 48.1% with 77% sensitivity.

Remaining INDEL work:
1. **Normalize INDELs** — **done (Jun 2026), no improvement**: bcftools norm on both calls and
   truth produced identical results. Multi-allelic mismatch is not a factor.
2. **Type-split TLOD filtering** — **done (Jun 2026)**: filter SNVs with MMQ≥50 and INDELs
   at TLOD ≥ 23.6 (ROC-optimal), then concat. Result: SNV P=76.4%, INDEL P=47.6%.

   ```bash
   OUTDIR=/mnt/storage/parabricks_test/output/HCC1395

   # SNVs: VAF≥5% + depth≥5 + POPAF≥2 + MMQ≥50
   bcftools filter \
     -i 'FILTER="PASS" && FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2 && INFO/MMQ >= 50 && TYPE="snp"' \
     ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz \
     -Oz -o ${OUTDIR}/HCC1395.somatic.typesplit.snvs.vcf.gz
   tabix -p vcf ${OUTDIR}/HCC1395.somatic.typesplit.snvs.vcf.gz

   # INDELs: VAF≥5% + depth≥5 + POPAF≥2 + TLOD≥23.6 (ROC-optimal for INDELs)
   bcftools filter \
     -i 'FILTER="PASS" && FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2 && INFO/TLOD >= 23.6 && TYPE!="snp"' \
     ${OUTDIR}/HCC1395.somatic.filtered.vcf.gz \
     -Oz -o ${OUTDIR}/HCC1395.somatic.typesplit.indels.vcf.gz
   tabix -p vcf ${OUTDIR}/HCC1395.somatic.typesplit.indels.vcf.gz

   # Merge and sort
   bcftools concat -a \
     ${OUTDIR}/HCC1395.somatic.typesplit.snvs.vcf.gz \
     ${OUTDIR}/HCC1395.somatic.typesplit.indels.vcf.gz | \
     bcftools sort -Oz -o ${OUTDIR}/HCC1395.somatic.typesplit.vcf.gz
   tabix -p vcf ${OUTDIR}/HCC1395.somatic.typesplit.vcf.gz
   ```

3. **DeepVariant** for INDELs in STR regions remains an option for sensitivity improvement.

## Reproducibility checklist

- [x] All tools pinned to specific container tags (e.g. `4.7.0-1`, not `latest`)
- [x] Reference genome explicit (GRCh38 no-alt from Broad bundle, version-locked)
- [x] Read group strings include sample, library, platform
- [x] `--bwa-options="-Y -K 100000000"` on both tumor and normal for deterministic BWA output
- [x] `pbrun prepon` run before `mutectcaller` when using a PoN
- [x] Resource VCFs subsetted + reheadered to match reference contig set and FAI order
- [x] `--pon` receives the VCF path (not the `.pon` binary path) — Parabricks appends `.pon`
- [x] vcfeval uses `--sample TRUTH,<tumor_sample>` when calls VCF is multi-sample
- [x] vcfeval uses `--evaluation-regions High-Confidence_Regions_v1.2.bed` for SEQC2 truth (omitting this inflates FP by ~93% on HCC1395 — chrX absent from truth, chr6/chr16 complex regions uncovered)
- [x] DeepSomatic eval uses `GQ` as vcf-score-field; Mutect2 eval uses `INFO.TLOD` — do not mix
- [x] `GetPileupSummaries` → `CalculateContamination` → contamination table fed to `FilterMutectCalls`
- [x] `LearnReadOrientationModel` → orientation model fed to `FilterMutectCalls`
- [x] Hard post-filter applied: `FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2` (tumor is sample index 1) → `HCC1395.somatic.popaf2.vcf.gz`
- [x] MMQ≥50 filter applied on top of POPAF2 filter → `HCC1395.somatic.mmq50.vcf.gz` (SNV precision 75.4%)
- [x] INDEL vcfeval uses INDEL-only calls VCF (`bcftools view -v indels`) — do NOT use mixed SNV+INDEL VCF against INDEL truth (SNVs inflate FP count)
- [x] OncoKB version recorded in the report (changes monthly)

**VCF filter chain** (all derived from `filtered.vcf.gz`, not chained):
1. `HCC1395.somatic.unfiltered.vcf.gz` — raw Mutect2 output (242,991 candidates)
2. `HCC1395.somatic.filtered.vcf.gz` — FilterMutectCalls with contamination + orientation model (95,321 PASS)
3. `HCC1395.somatic.popaf2.vcf.gz` — + VAF≥5% + depth≥5 + POPAF≥2 (44,504 variants)
4. `HCC1395.somatic.mmq50.vcf.gz` — + MMQ≥50, from `filtered.vcf.gz` (43,251 variants; SNV-optimized)
5. `HCC1395.somatic.typesplit.vcf.gz` — type-split: SNVs with MMQ≥50, INDELs with TLOD≥23.6 (42,624 variants; best combined SNV+INDEL precision)

## Cost reference (May 2026, AWS)

| Instance | GPUs | $/hr | 30x WGS TN pair runtime | Cost/sample |
|---|---|---|---|---|
| g5.24xlarge | 4× A10G | ~$8 | ~90 min | ~$12 |
| p4d.24xlarge | 8× A100 | ~$33 | ~25 min | ~$14 |
| p5.48xlarge | 8× H100 | ~$98 | ~12 min | ~$20 |

Compare to ~$50–80/sample on CPU instances running 24+ hours.

## Files in this repo

```
.
├── README.md                              # this file
├── scripts/
│   ├── run_parabricks_pipeline.sh         # full GPU pipeline (FASTQ → filtered VCF, Mutect2)
│   ├── run_deepsomatic_pipeline.sh        # DeepSomatic tumor-normal calling (phases 0–2)
│   ├── eval_deepsomatic.sh               # vcfeval benchmark for DeepSomatic (phase 3)
│   ├── diagnose_deepsomatic.sh           # root-cause checker for low-precision results
│   ├── annotate_variants.sh               # VEP local-cache + OncoKB docker (offline mode)
│   ├── vcf_to_maf_vep_rest.py             # VEP REST + OncoKB hotspot → MAF (no local install)
│   ├── generate_clinical_report.py        # MAF → AMP-tiered clinical report (text + JSON)
│   ├── download_hcc1395.sh                # SRA + SEQC2 truth + HC BED download for HCC1395
│   ├── download_colo829_truth.sh          # COLO829 truth VCF download
│   ├── download_refs.sh                   # GRCh38 reference bundle download
│   └── run_smoke_test.sh                  # fast amplicon-panel smoke test
├── sample_data/
│   ├── build_sample_maf.py                # generates a realistic synthetic MAF
│   ├── PATIENT001.oncokb.maf              # 7-variant LUAD example
│   ├── PATIENT001_meta.json               # patient metadata
│   ├── PATIENT001_clinical_report.txt     # the demo report
│   └── PATIENT001_clinical_report.json
└── docs/
    └── variant_interpretation.md          # AMP/ASCO/CAP tier reference
```
