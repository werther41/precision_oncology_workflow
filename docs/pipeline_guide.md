# Unified Precision Oncology Pipeline — User Guide

`run_pipeline.sh` + `pipeline.conf` is the single entry point for all four stages of somatic
variant analysis: alignment, variant calling, functional annotation, and clinical reporting.

---

## Table of Contents

1. [Pipeline Overview](#1-pipeline-overview)
2. [One-Time Setup](#2-one-time-setup)
3. [Configuration (`pipeline.conf`)](#3-configuration-pipelineconf)
4. [Running the Pipeline](#4-running-the-pipeline)
   - [From FASTQs (full run)](#41-from-fastqs-full-run)
   - [From existing BAMs (skip alignment)](#42-from-existing-bams-skip-alignment)
   - [From a pre-called VCF (annotation only)](#43-from-a-pre-called-vcf-annotation-only)
   - [Dry-run mode](#44-dry-run-mode)
5. [CLI Flag Reference](#5-cli-flag-reference)
6. [Batch Mode](#6-batch-mode)
7. [Output Layout](#7-output-layout)
8. [Idempotency and Stage Skipping](#8-idempotency-and-stage-skipping)
9. [PCGR Preprocessing Detail](#9-pcgr-preprocessing-detail)
10. [Validated Results — HCC1395](#10-validated-results--hcc1395)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Pipeline Overview

```
FASTQ R1/R2          BAM (BQSR)           PASS VCF             Reports
tumor + normal
    │                    │                    │                    │
    ▼                    ▼                    ▼                    ▼
┌──────────┐       ┌──────────────┐    ┌────────────┐     ┌─────────────┐
│ FQ2BAM   │──────>│ DeepSomatic  │───>│ Open CRAVAT│────>│    PCGR     │
│ Parabricks│       │ Parabricks  │    │ CIViC      │     │ ClinVar     │
│ BWA+BQSR │       │ GPU deep-   │    │ OncoKB     │     │ COSMIC      │
│ ~50 min  │       │ learning TN │    │ ~10 min    │     │ gnomAD      │
└──────────┘       │ caller      │    └────────────┘     │ VEP 112     │
                   │ ~50 min GPU │                        │ ~10 min     │
                   └──────────────┘                       └─────────────┘

run_pipeline.sh ──────────────────────────────────────────────────────────>
Entry point auto-detected from inputs:
  --tumor-r1/r2 → FQ2BAM  |  --tumor-bam → DeepSomatic  |  --vcf → OC+PCGR
```

**Stage summary:**

| Stage | Tool | Input → Output | Runtime |
|---|---|---|---|
| 1. FQ2BAM | Parabricks `pbrun fq2bam` | FASTQ.gz → BAM+BQSR | ~50 min (RTX PRO 6000) |
| 2. DeepSomatic | Parabricks `pbrun deepsomatic` | BAM pair → PASS VCF + MAF | ~50 min GPU |
| 3. Open CRAVAT | `oc run` (CIViC, OncoKB) | PASS VCF → annotated TSV | ~10 min |
| 4. PCGR | Docker `sigven/pcgr:2.2.5` | PASS VCF → HTML/XLSX report | ~10 min |

Each stage is **idempotent** — if its sentinel output already exists it is skipped. A failed
run can be safely re-submitted from the top; it resumes where it left off.

---

## 2. One-Time Setup

All dependencies are already installed on this machine. The table below documents what is
needed for a fresh install.

| Dependency | Status | Location / Command |
|---|---|---|
| Docker + nvidia-container-toolkit | ✅ installed | — |
| Parabricks 4.7.0-1 | ✅ pulled | `docker pull nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1` |
| Open CRAVAT 3.1.1 | ✅ installed | `/mnt/storage/open-cravat-reports/oc-venv/bin/oc` |
| Open CRAVAT annotators (CIViC, OncoKB) | ✅ installed | `oc module install civic oncokb` |
| PCGR 2.2.5 | ✅ pulled | `docker pull sigven/pcgr:2.2.5` |
| PCGR reference bundle (~5 GB) | ✅ downloaded | `/mnt/storage/pcgr_bundle` |
| VEP cache GRCh38 v112 (~25 GB) | ✅ downloaded | `/mnt/storage/vep_cache` |
| GRCh38 reference + resource bundle | ✅ present | `/mnt/storage/parabricks_test/ref/` |
| OncoKB API token | **configure** | `pipeline.conf` → `ONCOKB_TOKEN=<your_token>` |

### Installing PCGR bundle (if missing)

```bash
mkdir -p /mnt/storage/pcgr_bundle
curl -L https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz \
  | tar -xz -C /mnt/storage/pcgr_bundle
```

### Installing VEP cache (if missing)

```bash
mkdir -p /mnt/storage/vep_cache
docker run --rm \
  -v /mnt/storage/vep_cache:/vep_cache \
  sigven/pcgr:2.2.5 \
  bash -c "vep_install -a cf -s homo_sapiens -y GRCh38 --CACHE_VERSION 112 \
           --CACHEDIR /vep_cache --NO_BIOPERL --NO_HTSLIB --NO_TEST"
```

### OncoKB token setup

Register a free account at [oncokb.org](https://www.oncokb.org/account/register).
Add your token to `pipeline.conf`:

```bash
ONCOKB_TOKEN=your_token_here
```

Without a token, Open CRAVAT's OncoKB annotator returns empty results (pipeline still runs).

---

## 3. Configuration (`pipeline.conf`)

`pipeline.conf` sits next to `run_pipeline.sh` and holds machine-level defaults. Edit it once
per machine. Every value can be overridden per-run via a CLI flag or by exporting an environment
variable before calling the pipeline.

**Precedence (highest → lowest):**
```
exported env var  >  CLI flag (--outdir, --tumor-type, …)  >  pipeline.conf  >  built-in default
```

### Full annotated `pipeline.conf`

```bash
# ── Reference genome ────────────────────────────────────────────────────────────
REF_DIR=/mnt/storage/parabricks_test/ref
REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"
KNOWN_SITES="${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
GNOMAD_VCF="${REF_DIR}/af-only-gnomad.noalt.vcf.gz"   # Mutect2 germline prior
PON_VCF="${REF_DIR}/1000g_pon.noalt.vcf.gz"            # Panel of Normals

# ── Annotation tool paths ────────────────────────────────────────────────────────
PCGR_BUNDLE=/mnt/storage/pcgr_bundle      # PCGR reference data bundle (~5 GB)
VEP_CACHE=/mnt/storage/vep_cache          # VEP cache GRCh38 v112 (~25 GB)
OC_BIN=/mnt/storage/open-cravat-reports/oc-venv/bin/oc   # Open CRAVAT binary

# OncoKB API token — obtain at oncokb.org  DO NOT COMMIT A REAL VALUE
ONCOKB_TOKEN=

# ── GPU ─────────────────────────────────────────────────────────────────────────
# GPU 1 on this machine is reserved for VLLM — always use GPU 0
NUM_GPUS=1
GPU_DEVICE=0

# ── Docker image versions (pinned) ───────────────────────────────────────────────
PARABRICKS_IMAGE=nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
PCGR_IMAGE=sigven/pcgr:2.2.5

# ── Pipeline defaults (all overridable via CLI) ───────────────────────────────────
OUTDIR=/mnt/storage/parabricks_test/output   # per-sample subdirs created automatically

MODEL_TYPE=WGS           # DeepSomatic model: WGS | WES | pacbio | ont
TUMOR_TYPE=LUAD          # OncoTree code — override per sample with --tumor-type
TUMOR_PURITY=0.6         # Estimated tumor purity (0.0–1.0), passed to PCGR
PCGR_ASSAY=WGS           # WGS | WES | TARGETED
WORKERS=5                # Concurrent VEP REST workers (DeepSomatic annotation phase)

OC_ANNOTATORS="civic oncokb"   # Space-separated Open CRAVAT module names
```

---

## 4. Running the Pipeline

### 4.1 From FASTQs (full run)

Runs all four stages: FQ2BAM → DeepSomatic → Open CRAVAT → PCGR.

```bash
bash run_pipeline.sh PATIENT001 \
  --tumor-r1  /data/PATIENT001/tumor_R1.fastq.gz \
  --tumor-r2  /data/PATIENT001/tumor_R2.fastq.gz \
  --normal-r1 /data/PATIENT001/normal_R1.fastq.gz \
  --normal-r2 /data/PATIENT001/normal_R2.fastq.gz \
  --tumor-type LUAD
```

Runtime: ~2 hours total on 1× RTX PRO 6000 Blackwell (96 GB).

### 4.2 From existing BAMs (skip alignment)

Skips FQ2BAM when BAMs are already aligned and BQSR-corrected. Runs: DeepSomatic → Open CRAVAT → PCGR.

```bash
bash run_pipeline.sh HCC1395 \
  --tumor-bam  /mnt/storage/parabricks_test/output/HCC1395/HCC1395_tumor.bam \
  --normal-bam /mnt/storage/parabricks_test/output/HCC1395/HCC1395_normal.bam \
  --tumor-type BRCA
```

### 4.3 From a pre-called VCF (annotation only)

Skips FQ2BAM and DeepSomatic. Runs: Open CRAVAT → PCGR.

The VCF must be:
- Filtered to PASS variants only (no multi-filter records)
- bgzipped and tabix-indexed (`.vcf.gz` + `.vcf.gz.tbi`)
- Called against GRCh38 (hg38)

```bash
bash run_pipeline.sh PATIENT002 \
  --vcf /data/PATIENT002.deepsomatic.pass.vcf.gz \
  --tumor-type PAAD
```

### 4.4 Dry-run mode

Prints all stage commands without executing anything. Use to verify stage detection, sentinel
paths, and Docker invocations before a real run.

```bash
bash run_pipeline.sh HCC1395 \
  --tumor-bam /mnt/storage/parabricks_test/output/HCC1395/HCC1395_tumor.bam \
  --normal-bam /mnt/storage/parabricks_test/output/HCC1395/HCC1395_normal.bam \
  --tumor-type BRCA \
  --dry-run
```

Expected output:
```
[HH:MM:SS] ===== Precision Oncology Pipeline — HCC1395 =====
[HH:MM:SS] Config:       /…/pipeline.conf
[HH:MM:SS] Start stage:  deepsomatic
[HH:MM:SS] ===== Stage 2/4: DeepSomatic =====
[HH:MM:SS] [DRY-RUN] Would run: env SAMPLE=HCC1395 … bash scripts/run_deepsomatic_pipeline.sh …
…
```

### 4.5 Partial pipeline (start-from / stop-after)

Override auto-detection to run a specific range of stages:

```bash
# Re-run only Open CRAVAT and PCGR (skip everything earlier)
bash run_pipeline.sh HCC1395 \
  --vcf /mnt/storage/parabricks_test/output/HCC1395/deepsomatic/HCC1395.deepsomatic.pass.vcf.gz \
  --start-from opencravat --stop-after pcgr \
  --tumor-type BRCA

# Run alignment only (stop before calling)
bash run_pipeline.sh PATIENT001 \
  --tumor-r1 … --tumor-r2 … --normal-r1 … --normal-r2 … \
  --stop-after fq2bam

# Skip annotation stages entirely
bash run_pipeline.sh HCC1395 \
  --tumor-bam … --normal-bam … \
  --stop-after deepsomatic
```

---

## 5. CLI Flag Reference

```
Usage: run_pipeline.sh SAMPLE_ID [OPTIONS]
       run_pipeline.sh --manifest samples.tsv [OPTIONS]
```

### Input flags

| Flag | Description |
|---|---|
| `SAMPLE_ID` (positional) | Sample identifier. Used in all output filenames and subdirectories. |
| `--tumor-r1 FILE` | Tumor R1 FASTQ.gz |
| `--tumor-r2 FILE` | Tumor R2 FASTQ.gz |
| `--normal-r1 FILE` | Normal R1 FASTQ.gz |
| `--normal-r2 FILE` | Normal R2 FASTQ.gz |
| `--tumor-bam FILE` | Pre-aligned tumor BAM (skips FQ2BAM) |
| `--normal-bam FILE` | Pre-aligned normal BAM |
| `--vcf FILE` | Pre-called PASS VCF.gz (skips FQ2BAM + DeepSomatic) |

### Pipeline control flags

| Flag | Default | Description |
|---|---|---|
| `--config FILE` | `./pipeline.conf` | Alternate config file path |
| `--outdir DIR` | from config | Override base output directory |
| `--tumor-type CODE` | `LUAD` | OncoTree code (BRCA, LUAD, PAAD, COAD, …) |
| `--tumor-purity FLOAT` | `0.6` | Estimated purity for PCGR (0.0–1.0) |
| `--num-gpus N` | `1` | Number of GPUs for Parabricks |
| `--gpu-device D` | `0` | GPU device index (GPU 1 reserved for VLLM) |
| `--start-from STAGE` | auto-detect | `fq2bam` \| `deepsomatic` \| `opencravat` \| `pcgr` |
| `--stop-after STAGE` | `pcgr` | Stop pipeline after this stage |
| `--skip-oc` | off | Skip the Open CRAVAT stage entirely |
| `--skip-pcgr` | off | Skip the PCGR stage entirely |
| `--manifest FILE` | — | Batch mode: TSV manifest (see §6) |
| `--dry-run` | off | Print commands; execute nothing |
| `-h`, `--help` | — | Show usage |

### Auto-detection of start stage

When `--start-from` is not given, the pipeline infers the entry point:

```
--vcf provided              →  start from opencravat
--tumor-bam + --normal-bam  →  start from deepsomatic
--tumor-r1/r2 provided      →  start from fq2bam
nothing provided            →  error
```

---

## 6. Batch Mode

Create a tab-delimited manifest. Lines beginning with `#` are skipped. Use `-` for unused fields.

### `samples.tsv` format

```
sample_id  tumor_r1  tumor_r2  normal_r1  normal_r2  tumor_type  tumor_bam  normal_bam  vcf
```

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `tumor_r1` | Tumor R1 FASTQ.gz path, or `-` |
| `tumor_r2` | Tumor R2 FASTQ.gz path, or `-` |
| `normal_r1` | Normal R1 FASTQ.gz path, or `-` |
| `normal_r2` | Normal R2 FASTQ.gz path, or `-` |
| `tumor_type` | OncoTree code, or `-` to use the config default |
| `tumor_bam` | Pre-aligned tumor BAM path, or `-` |
| `normal_bam` | Pre-aligned normal BAM path, or `-` |
| `vcf` | Pre-called PASS VCF.gz path, or `-` |

### Example manifest

```tsv
sample_id	tumor_r1	tumor_r2	normal_r1	normal_r2	tumor_type	tumor_bam	normal_bam	vcf
# Full pipeline from FASTQs
PT001	/data/PT001_tumor_R1.fq.gz	/data/PT001_tumor_R2.fq.gz	/data/PT001_N_R1.fq.gz	/data/PT001_N_R2.fq.gz	LUAD	-	-	-
# Start from existing BAMs
HCC1395	-	-	-	-	BRCA	/mnt/storage/.../HCC1395_tumor.bam	/mnt/storage/.../HCC1395_normal.bam	-
# Annotation only
PT003	-	-	-	-	PAAD	-	-	/data/pt003.pass.vcf.gz
```

### Running the batch

```bash
bash run_pipeline.sh --manifest samples.tsv --outdir /mnt/storage/results
```

Samples run **sequentially** (the GPU is the bottleneck; a single WGS pair saturates one GPU).
Each sample gets its own log at `${OUTDIR}/${SAMPLE}/pipeline_run.log`.

Pass `--dry-run` to preview all commands for all samples before executing:

```bash
bash run_pipeline.sh --manifest samples.tsv --outdir /mnt/storage/results --dry-run
```

---

## 7. Output Layout

```
${OUTDIR}/${SAMPLE}/
│
├── ${SAMPLE}_tumor.bam                        ← FQ2BAM sentinel
├── ${SAMPLE}_tumor.bam.bai
├── ${SAMPLE}_normal.bam
├── ${SAMPLE}_normal.bam.bai
│
├── deepsomatic/
│   ├── ${SAMPLE}.deepsomatic.vcf.gz           Raw DeepSomatic calls (all filters)
│   ├── ${SAMPLE}.deepsomatic.pass.vcf.gz      ← DeepSomatic sentinel (PASS only)
│   ├── ${SAMPLE}.deepsomatic.maf              VEP + OncoKB annotated MAF
│   ├── ${SAMPLE}.deepsomatic_clinical_report.txt   AMP-tiered text report
│   ├── ${SAMPLE}.deepsomatic_clinical_report.json  EHR-ready JSON
│   │
│   └── pcgr/
│       ├── ${SAMPLE}.pcgr_input.vcf.gz        Preprocessed VCF (TDP/TVAF lifted)
│       ├── ${SAMPLE}.pcgr_input.vcf.gz.tbi
│       ├── ${SAMPLE}.pcgr.grch38.html         ← PCGR sentinel (full HTML report, ~16 MB)
│       ├── ${SAMPLE}.pcgr.grch38.xlsx         Excel summary
│       ├── ${SAMPLE}.pcgr.grch38.snv_indel_ann.tsv.gz   Full annotation TSV
│       └── ${SAMPLE}.pcgr.grch38.tmb.tsv      TMB metrics
│
├── opencravat/
│   ├── ${SAMPLE}-oc-report.tsv                ← Open CRAVAT sentinel
│   ├── ${SAMPLE}-oc-report.sqlite             Full OC annotation database
│   ├── ${SAMPLE}-oc-report.variant.tsv        Per-variant flat TSV
│   └── ${SAMPLE}.pass.vcf                     Intermediate plain VCF (cleaned up)
│
└── pipeline_run.log                           Combined log for all stages
```

### Key sentinel files

The pipeline skips a stage if its sentinel file already exists:

| Stage | Sentinel |
|---|---|
| FQ2BAM | `${SAMPLE}_tumor.bam` AND `${SAMPLE}_normal.bam` |
| DeepSomatic | `deepsomatic/${SAMPLE}.deepsomatic.pass.vcf.gz` |
| Open CRAVAT | `opencravat/${SAMPLE}-oc-report.tsv` |
| PCGR | `deepsomatic/pcgr/${SAMPLE}.pcgr.grch38.html` |

---

## 8. Idempotency and Stage Skipping

Every stage checks for its sentinel file at entry. If it exists:

```
[HH:MM:SS] Skipping DeepSomatic — VCF exists: …/HCC1395.deepsomatic.pass.vcf.gz
```

This means:
- A failed run at any stage can be re-submitted from the top — only the failed and subsequent stages will re-run.
- Adding new samples to a batch manifest is safe: previously completed samples skip all stages instantly.
- Re-running after updating `pipeline.conf` (e.g., changing `TUMOR_TYPE`) will still skip completed stages. To force a re-run of a stage, delete its sentinel file.

**Forcing a stage to re-run:**

```bash
# Re-run PCGR only
rm /mnt/storage/parabricks_test/output/HCC1395/deepsomatic/pcgr/HCC1395.pcgr.grch38.html
bash run_pipeline.sh HCC1395 \
  --vcf .../HCC1395.deepsomatic.pass.vcf.gz \
  --tumor-type BRCA

# Re-run Open CRAVAT and PCGR
rm /mnt/storage/.../opencravat/HCC1395-oc-report.tsv
rm /mnt/storage/.../deepsomatic/pcgr/HCC1395.pcgr.grch38.html
bash run_pipeline.sh HCC1395 --vcf .../HCC1395.deepsomatic.pass.vcf.gz --tumor-type BRCA
```

---

## 9. PCGR Preprocessing Detail

### The problem

DeepSomatic stores per-variant depth and allele fraction as VCF **FORMAT** fields:
- `FORMAT/DP` — total read depth at variant position
- `FORMAT/VAF` — variant allele fraction

PCGR 2.2.5 requires depth and VAF as **INFO** fields and rejects `--tumor_dp_tag DP` with:

```
ERROR: Custom INFO tag (tumor_dp_tag) needs another name —
       'DP' is a reserved field in the VCF specification (INFO)
```

### The fix

`scripts/vcf_add_info_dp_vaf.awk` lifts the FORMAT values into new INFO fields with safe names:

```
FORMAT/DP  → INFO/TDP   (Tumor Depth)
FORMAT/VAF → INFO/TVAF  (Tumor VAF)
```

This preprocessing step runs automatically inside the PCGR Docker container (which has
`bcftools` and `bgzip`). The result is written to `${SAMPLE}.pcgr_input.vcf.gz` and used as
input to PCGR with `--tumor_dp_tag TDP --tumor_af_tag TVAF`.

The original PASS VCF is not modified. The preprocessing step is idempotent: if
`${SAMPLE}.pcgr_input.vcf.gz` already exists it is reused.

### AWK script (`scripts/vcf_add_info_dp_vaf.awk`)

```awk
BEGIN { OFS = "\t"; dp_hdr = 0; vaf_hdr = 0 }
/^##FORMAT=<ID=DP,/ {
    if (!dp_hdr) {
        print "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Tumor depth from FORMAT/DP\">"
        dp_hdr = 1
    }
    print; next
}
/^##FORMAT=<ID=VAF,/ {
    if (!vaf_hdr) {
        print "##INFO=<ID=TVAF,Number=A,Type=Float,Description=\"Tumor VAF from FORMAT/VAF\">"
        vaf_hdr = 1
    }
    print; next
}
/^#/ { print; next }
{
    n = split($9, fmt, ":"); split($10, vals, ":")
    dp = ""; vaf = ""
    for (i = 1; i <= n; i++) {
        if (fmt[i] == "DP")  dp  = vals[i]
        if (fmt[i] == "VAF") vaf = vals[i]
    }
    tag = ""
    if (dp  != "" && dp  != ".") tag = tag (tag=="" ? "" : ";") "TDP="  dp
    if (vaf != "" && vaf != ".") tag = tag (tag=="" ? "" : ";") "TVAF=" vaf
    if ($8 == ".") $8 = (tag=="" ? "." : tag)
    else           $8 = $8 (tag=="" ? "" : ";" tag)
    print
}
```

---

## 10. Validated Results — HCC1395

Complete end-to-end run on HCC1395 (SEQC2 FDA breast cancer benchmark), June 2026.
Platform: 1× RTX PRO 6000 Blackwell (97.8 GB), GRCh38 no-alt, Parabricks 4.7.0-1.

### Runtime

| Stage | Tool | Runtime |
|---|---|---|
| FQ2BAM (tumor) | `pbrun fq2bam` | ~15 min |
| FQ2BAM (normal) | `pbrun fq2bam` | ~15 min |
| DeepSomatic | `pbrun deepsomatic` | ~50 min |
| VEP/OncoKB annotation | VEP REST (5 workers) | ~25 min |
| Open CRAVAT | `oc run` (civic, oncokb) | ~10 min |
| PCGR | Docker `sigven/pcgr:2.2.5` | ~8 min |

### Somatic call summary

| Metric | Value |
|---|---|
| Total DeepSomatic PASS calls | 124,014 |
| SNVs | ~122,000 |
| INDELs | ~2,000 |
| TMB (missense only) | 13.8 mut/Mb |
| TMB (coding non-silent) | 15.5 mut/Mb |

### Annotation findings

| Stage | Finding | Detail |
|---|---|---|
| Open CRAVAT — CIViC | TP53 p.R175H | Confirmed oncogenic; lung + breast cancer clinical evidence |
| Open CRAVAT — OncoKB | (requires valid token) | Set `ONCOKB_TOKEN` in `pipeline.conf` |
| PCGR — Tier 1 | TP53 p.R175H | ClinVar Pathogenic; direct evidence |
| PCGR — Tier IA | BRCA2 frameshift | Tumor suppressor LoF → AMP/ESMO Tier IA by rule |
| PCGR — Tier IA | BRCA1 truncating variant | Tumor suppressor LoF |

### Benchmark against SEQC2 truth (DeepSomatic)

Run `scripts/eval_deepsomatic.sh HCC1395` after calling to evaluate against the SEQC2 truth VCF.
**The HC regions BED is mandatory** — omitting it inflates false positives from chrX, chr6 MHC,
and chr16 repeat regions and drops precision from ~93% to ~30%.

```bash
./scripts/eval_deepsomatic.sh HCC1395
# Output: deepsomatic/eval/comparison.md
```

---

## 11. Troubleshooting

### "PCGR_BUNDLE not set or not found"

PCGR is skipped with a warning. Download the bundle:

```bash
mkdir -p /mnt/storage/pcgr_bundle
curl -L https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz \
  | tar -xz -C /mnt/storage/pcgr_bundle
```

Then set `PCGR_BUNDLE=/mnt/storage/pcgr_bundle` in `pipeline.conf`.

### "Open CRAVAT binary not found or not executable"

Open CRAVAT is skipped with a warning. Install it:

```bash
python3 -m venv /opt/oc-venv
/opt/oc-venv/bin/pip install open-cravat
/opt/oc-venv/bin/oc module install civic oncokb
```

Then set `OC_BIN=/opt/oc-venv/bin/oc` in `pipeline.conf`.

### "Cannot determine start stage"

You must provide at least one of: `--tumor-r1/r2` (FASTQs), `--tumor-bam/--normal-bam` (BAMs),
or `--vcf` (PASS VCF). If none are provided the pipeline cannot auto-detect an entry point.

### PCGR Docker fails with "DP is a reserved field"

This error should not occur with the current `run_pcgr.sh` (which uses `--tumor_dp_tag TDP`).
If you call `run_pcgr.sh` directly with an older version, ensure you are using the script at
`scripts/run_pcgr.sh` from this repo (not a local copy from before the preprocessing fix).

### OncoKB returns empty results

Expected when `ONCOKB_TOKEN` is not set. Set it in `pipeline.conf`:

```bash
ONCOKB_TOKEN=your_token_here
```

CIViC annotations (Open CRAVAT) and PCGR ClinVar/COSMIC tiers are unaffected by the token.

### GPU error: "device 1 not available"

`pipeline.conf` defaults to `GPU_DEVICE=0`. If you see a GPU 1 error, check that nothing
is overriding `GPU_DEVICE` to `1`. GPU 1 on this machine is reserved for VLLM.

### Stage re-runs on restart

If you want a stage to re-run after it already completed, delete its sentinel file (see §8).

### Log files

Per-sample log: `${OUTDIR}/${SAMPLE}/pipeline_run.log`  
Check it for exact error messages and timestamps for each stage.

---

*Generated June 2026. Pipeline version: `run_pipeline.sh` + `pipeline.conf` v2.0.*
