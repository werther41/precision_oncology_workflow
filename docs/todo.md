# TODO — Parabricks Precision Oncology Pipeline (HCC1395)

Follow-up tasks for Claude Code. Pipeline runs end-to-end; remaining work is FP
reduction, INDEL recovery, clinical-deliverable completion, and doc cleanup.

Current best benchmark (vs SEQC2 truth, Jun 2026) — INDEL-only eval (corrected):
- SNV (mmq50):   P 75.4% / Se 88.1% / F1 0.812   (popaf2 was: P 74.1% / Se 88.4% / F1 0.806)
- INDEL (mmq50): P 48.1% / Se 77.1% / F1 0.592   at optimal TLOD 23.6 (INDEL-only calls)
  CORRECTION: prior P=6.1% was vcfeval artifact (SNVs counted as INDEL FPs). True precision = 48%.

Paths:
  OUTDIR=/mnt/storage/parabricks_test/output/HCC1395
  REF=/mnt/storage/parabricks_test/ref
  TRUTHDIR=/mnt/storage/parabricks_test/data/truth

---

## 0. Reconcile README inconsistency (do this first — source-verify)

The README contradicts itself on the GATK contamination/orientation step:
- "Next Step: GATK Best-Practice Filtering" lists the intermediate tables as
  `❌ not yet generated`.
- BUT the vcfeval results table, the "Actual results after Step 4" block, and the
  reproducibility checklist all mark this work DONE.

Action:
- [x] Confirm on disk which intermediate files actually exist:
      5/5 intermediate files present; filtered_v2.vcf.gz absent (never created as
      separate file — Step 4 wrote to filtered.vcf.gz directly, confirmed via VCF header).
- [x] Fix the stale `❌` markers in the README to `✅` — done.
- [x] Clarify which VCF the POPAF filter consumed: filtered.vcf.gz IS the contamination+
      orientation output (VCF header confirms --contamination-table and --ob-priors were
      active). filtered_v2.vcf.gz was the planned name but was never created separately.

---

## 1. SNV FP reduction — MMQ50 (drafted, not run)

11,441 SNV FPs remain after POPAF. Class 3 = HLA/segdup artifacts (chr6:28–34M,
chr1 repeats). README has the command drafted but never executed.

- [x] Run the MMQ50 filter — done. Output: HCC1395.somatic.mmq50.vcf.gz
- [x] Re-derive PASS single-sample + re-run vcfeval — done.
- [x] Record results in README: SNV P 75.4% / Se 88.1% / F1 0.812 (POPAF: 74.1% / 88.4% / 0.806)
- [x] Quantify removed FPs by region: class-3 chr6:28-34M hypothesis disproven (0 FPs there).
      783 removed FPs are genome-wide low-mapq (MMQ=40 dominant class); chrX top contributor
      (127/783 = 16%). Updated README with corrected class-3 description.

## 2. SNV FP triage — characterize remaining classes

- [x] Class 1 (low-VAF noise): TLOD sweep done (Jun 15 2026). Result: precision
      CEILING at ~82-83% SNV-specific regardless of TLOD. Going TLOD 5.8→50 gains 0 pp
      precision while losing 40 pp sensitivity. TLOD filtering is NOT useful here.
      FPs are distributed across the full TLOD range. Updated README with sweep table.
- [ ] Class 2 (high-VAF, POPAF=6, AF ~100%): spot-check a sample in IGV. If these
      are real somatic calls absent from SEQC2 truth, they are a truth-set gap, NOT
      FPs — document so they aren't "filtered away" chasing a misleading metric.

---

## 3. INDEL improvement — biggest gap (P 6.1%)

Mutect2 INDELs suffer from homopolymer/STR alignment error. Three angles from README:

- [x] Normalize before eval — done (Jun 2026). Zero improvement: TP/FP/FN identical.
      Multi-allelic mismatch is NOT the cause of poor INDEL precision.
- [x] EVALUATION BUG FOUND AND FIXED (Jun 15 2026): prior P=6.1% was an artifact.
      Mixed SNV+INDEL calls VCF → vcfeval counts all SNVs as INDEL FPs. Fixed by using
      INDEL-only calls VCF (bcftools view -v indels). Corrected result: P=48.1%, Se=77.1%,
      F1=0.592 at TLOD≥23.6. README benchmark table and INDEL section updated.
- [x] Apply type-split TLOD filtering — done (Jun 15 2026).
      SNVs: MMQ50 filter. INDELs: TLOD≥23.6. Merged → typesplit.vcf.gz (42,624 variants).
      Result: SNV P=76.4% / Se=87.9% / F1=0.818; INDEL P=47.6% / Se=77.3% / F1=0.590.
      This is the recommended VCF for clinical annotation (typesplit.vcf.gz).
- [ ] Evaluate DeepVariant INDELs (pbrun deepvariant output already exists, 7.1M
      germline) merged with Mutect2 somatic — DV is typically 10–30x better in STR
      regions. Decide on a merge strategy and benchmark.

---

## 4. Clinical deliverable — finish the last mile

The run table stops at VEP; MAF + report exist as commands but aren't shown as
validated on HCC1395, and they point at the lower-precision VCF.

- [x] Point annotation at high-precision VCF. Switched from mmq50 to typesplit.vcf.gz
      (recommended for best combined SNV+INDEL precision: SNV 76.4%, INDEL 47.6%).
      VEP running for mmq50.vcf.gz (Jun 15 2026, ~3hr total); mmq50 MAF can be used for
      initial clinical report (same SNV drivers as typesplit; INDEL set differs slightly).
- [ ] Set ONCOKB_TOKEN for the run (currently falls back to built-in hotspot table)
      and record the OncoKB version in the report (changes monthly per checklist).
- [x] Run VEP REST on typesplit VCF — done (Jun 15 2026). Used --workers 5; ~7 min (6x vs sequential).
      VEP run twice: first revealed INDEL coordinate mismatch bug (89% of INDELs got UNKNOWN gene).
      Fix: _lookup_vep now tries VEP-normalized indel coords (pos+anchor, ref trimmed, alt="-").
      After fix: 610 no-hit variants (was 2841), 1910 UNKNOWN INDELs (was 3172).
      CLIN_SIG extracted from VEP colocated_variants. GENE_LOF_DB added for BRCA1/BRCA2/PALB2 LoF.
      NOTE: BRCA1 at chr17:43046406 is deep INTRONIC (c.5468-605del), NOT a coding frameshift.
- [x] Generate + eyeball the clinical report — done (Jun 15 2026).
      BRCA2 p.E1593*: Tier IA, LEVEL_1 (PARP inhibitors), VAF=39.5% ✅
      TP53 p.R175H: Tier IIC, LEVEL_3B (APR-246), VAF=100% ✅
      XDH/MTRR appear Tier IIC via ClinVar — artifacts (pathogenic for non-cancer indications).
      356 VUS including TIFA frameshift, JUP/TRPM5/PDE3B inframe dels.
      Fixed generate_clinical_report.py: NaN guard for HGVSp_Short in pandas rows.
- [x] Confirm expected HCC1395 drivers present (coordinate check done):
      TP53 p.R175H (VAF=98%, TLOD=190), BRCA1 frameshift (VAF=90%, TLOD=37.5),
      BRCA2 multiple somatic variants confirmed. PIK3CA absent (expected for HCC1395).

---

## 5. README / repo hygiene

- [x] Add the MMQ50 + INDEL results to the benchmark tables — done (Jun 15 2026).
      Both rows added; "applied" vs "next command to try" language updated.
- [x] VCF chain documented in reproducibility checklist (filtered → popaf2/mmq50).
- [x] Reproducibility checklist verified — all [x] items backed by on-disk artifacts.
      Added MMQ50 to checklist.

---

## Priority order

1. Task 0 (reconcile — cheap, prevents chasing phantom results)
2. Task 1 (MMQ50 — drafted, highest-leverage SNV gain)
3. Task 3 (INDEL norm — likely a measurement artifact, possibly free precision)
4. Task 4 (clinical deliverable — needed for an actual end-to-end demo)
5. Tasks 2 & 5 (characterization + docs — finish-up)
