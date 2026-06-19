# TODO — Parabricks Precision Oncology Pipeline (HCC1395)

Follow-up tasks for Claude Code. Pipeline runs end-to-end and is benchmarked.
Remaining work: truth-set gap verification, INDEL sensitivity, clinical-call
correctness, and doc hygiene. Updated after the Jun 15 2026 work pass.

Current best benchmark (vs SEQC2 truth):
- SNV  (typesplit): P 76.4% / Se 87.9% / F1 0.818
- INDEL (typesplit, INDEL-only eval): P 47.6% / Se 77.3% / F1 0.590
- Recommended clinical-annotation VCF: HCC1395.somatic.typesplit.vcf.gz

Paths:
  OUTDIR=/mnt/storage/parabricks_test/output/HCC1395
  REF=/mnt/storage/parabricks_test/ref
  TRUTHDIR=/mnt/storage/parabricks_test/data/truth

================================================================================
DONE since last pass (Jun 15 2026)
================================================================================
- [x] Task 0 — reconciled: 5/5 contamination/orientation tables on disk;
      filtered_v2 was never separate (Step 4 overwrote filtered.vcf.gz, header
      confirms --contamination-table + --ob-priors). Stale README markers fixed.
- [x] Task 1 — MMQ50 run: SNV P 75.4% / Se 88.1% / F1 0.812. HLA/segdup
      hypothesis DISPROVEN (0 FPs chr6:28-34M); removed FPs are genome-wide
      low-mapq (MMQ=40, chrX top contributor 16%).
- [x] Task 2 Class 1 — TLOD sweep: precision ceiling ~82-83% SNV-specific,
      flat across TLOD 5.8->50. TLOD filtering is a dead end here.
- [x] Task 3 — INDEL "gap" was a vcfeval artifact (mixed SNV+INDEL VCF vs
      INDEL truth inflated FP by ~39K). Corrected: P 48.1% / Se 77.1% at
      TLOD>=23.6. Normalization gave zero change. Type-split VCF built.
- [x] Task 4 (partial) — VEP concurrency (~7 min, 5 workers) + INDEL coord-norm
      bug fixed (UNKNOWN INDELs 89% -> ~38%). Clinical report generated;
      BRCA2 p.E1593* Tier IA and TP53 p.R175H Tier IIC verified.
- [x] Task 5 — README benchmark tables, VCF chain, and reproducibility
      checklist updated + verified against on-disk artifacts.

================================================================================
OPEN — substantive
================================================================================

## A. Class-2 truth-set-gap verification  (HIGHEST VALUE — do first)
The ~82-83% SNV precision ceiling may be partly truth-set gaps, not real FPs.
- [ ] IGV spot-check a sample of high-VAF / POPAF=6 calls (tumor AF ~100%,
      absent from SEQC2 truth). Pull a few from the FP set:
      bcftools isec confirms which "FP" calls are POPAF=6 high-VAF.
- [ ] If they're clean somatic calls missing from truth, document them as a
      truth-set limitation and report precision as a LOWER BOUND. This changes
      how the ~83% ceiling is interpreted (you may not be FP-bound at all).
- [ ] Decide: stop filter-tuning SNVs once confirmed (no headroom worth chasing).

## B. DeepVariant INDEL merge  (only remaining INDEL lever)
- [ ] Use the existing pbrun deepvariant output; intersect/merge somatic INDELs
      with Mutect2 in STR regions. Benchmark Se gain vs precision cost.
- [ ] Decide a merge strategy (DV-supports-Mutect2 vs union) and record results.

## C. OncoKB token
- [ ] Set ONCOKB_TOKEN (currently built-in hotspot fallback) and re-annotate.
- [ ] Record OncoKB data version in the report (changes monthly).

================================================================================
VERIFY — clinical-call correctness (NEW — surfaced from report review)
================================================================================

## D. BRCA1: possible false-actionable
VEP note: chr17:43046406 is deep INTRONIC (c.5468-605del), NOT coding. But the
driver check reports "BRCA1 frameshift (VAF 90%)" and GENE_LOF_DB adds BRCA1 LoF.
- [ ] Confirm whether the BRCA1 LoF annotation is firing on the intronic variant.
- [ ] If so, GENE_LOF_DB is mislabeling an intronic call as actionable LoF —
      fix the LoF gate to require a coding consequence. A spurious Tier IA/PARP
      call in a clinical report is a hard error.

## E. ClinVar non-cancer leakage into tiering
XDH/MTRR surfaced Tier IIC via ClinVar but are pathogenic for non-cancer
indications (caught manually this run).
- [ ] Add disease-context filtering so generate_clinical_report.py doesn't tier
      ClinVar hits whose pathogenicity is unrelated to oncology. Prevents
      recurrence on the next sample.

================================================================================
DOC HYGIENE
================================================================================
- [ ] VEP runtime stated inconsistently: ~7 min (architecture + Stage 3),
      ~25 min (validated-run table), ~25-40 min (Quick Start comment).
      Pick the real concurrent number (~7 min) and make all four agree.
- [ ] GATK "Next Step" table still lists expected "SNV 85-92% / INDEL ~70%"
      impact that the actual run (no change) disproved in the section directly
      below. Reword to "expected (not realized — see results below)".

================================================================================
Priority order
================================================================================
1. A (truth-set gap) — reframes the SNV ceiling; cheap; prevents wasted tuning
2. D (BRCA1 false-actionable) — clinical correctness, hard error if real
3. E (ClinVar leakage) — clinical correctness, recurs every sample
4. C (OncoKB token) — needed for a real clinical run
5. B (DeepVariant INDELs) — sensitivity upside, larger effort
6. Doc hygiene — finish-up