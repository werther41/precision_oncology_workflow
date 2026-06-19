#!/usr/bin/env bash
#
# RTG vcfeval benchmark for DeepSomatic calls (Phase 3)
# Reuses SEQC2 truth VCFs already prepared at /mnt/storage/parabricks_test/data/truth/
#
# Evaluates:
#   - SNV:   deepsomatic.pass.vcf.gz  vs  high-confidence_sSNV_gt.vcf.gz
#   - INDEL: deepsomatic.pass.indels.vcf.gz  vs  high-confidence_sINDEL_gt.vcf.gz
#
# Writes: <outdir>/deepsomatic/eval/comparison.md
#
# Usage:
#   ./eval_deepsomatic.sh <sample_id>
#
# Optional env overrides:
#   TRUTHDIR=/path/truth   SEQC2 truth VCF directory
#   OUTDIR=/path/out       pipeline output root (deepsomatic/ subdir must exist)
#   SDF=/path/ref.sdf      RTG SDF template for GRCh38
#   HC_BED=/path/hc.bed    SEQC2 High-Confidence regions BED (required for correct precision)

set -euo pipefail

SAMPLE="${1:-HCC1395}"
TRUTHDIR="${TRUTHDIR:-/mnt/storage/parabricks_test/data/truth}"
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/${SAMPLE}}"
DS_OUTDIR="${OUTDIR}/deepsomatic"
EVAL_DIR="${DS_OUTDIR}/eval"
SDF="${SDF:-/mnt/storage/parabricks_test/ref/GRCh38.sdf}"
HC_BED="${HC_BED:-${TRUTHDIR}/High-Confidence_Regions_v1.2.bed}"

RTG_IMAGE="realtimegenomics/rtg-tools:3.12.1"
BCFTOOLS_IMAGE="staphb/bcftools:1.19"

mkdir -p "${EVAL_DIR}"

# ============================================================
# Read sample name detected in Phase 2
# ============================================================
SAMPLE_NAME_FILE="${DS_OUTDIR}/${SAMPLE}.sample_name.txt"
if [[ ! -f "${SAMPLE_NAME_FILE}" ]]; then
  echo "ERROR: Sample name file not found: ${SAMPLE_NAME_FILE}"
  echo "       Run run_deepsomatic_pipeline.sh first to complete Phase 2."
  exit 1
fi
TUMOR_SAMPLE=$(cat "${SAMPLE_NAME_FILE}")
echo "[$(date +%H:%M:%S)] Tumor sample name: ${TUMOR_SAMPLE}"

# ============================================================
# Input VCF paths
# Truth VCFs can be overridden via SNV_TRUTH / INDEL_TRUTH env vars.
# Default: SEQC2 HCC1395 truth. For COLO829 pass:
#   SNV_TRUTH=/mnt/.../colo829/colo829_sSNV_gt.vcf.gz
#   INDEL_TRUTH=/mnt/.../colo829/colo829_sINDEL_gt.vcf.gz
# ============================================================
DS_PASS="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.vcf.gz"
DS_INDELS="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.indels.vcf.gz"
SNV_TRUTH="${SNV_TRUTH:-${TRUTHDIR}/high-confidence_sSNV_gt.vcf.gz}"
INDEL_TRUTH="${INDEL_TRUTH:-${TRUTHDIR}/high-confidence_sINDEL_gt.vcf.gz}"
TRUTH_LABEL="${TRUTH_LABEL:-SEQC2 HCC1395 v1.2 (FDA benchmark)}"

for f in "${DS_PASS}" "${DS_INDELS}" "${SNV_TRUTH}" "${INDEL_TRUTH}" "${SDF}"; do
  [[ -e "${f}" ]] || { echo "ERROR: Required file/dir not found: ${f}"; exit 1; }
done

# HC_BED is optional: use it if found, warn and run without it if not.
HC_BED_FLAG=""
if [[ -n "${HC_BED}" && -e "${HC_BED}" ]]; then
  HC_BED_FLAG="--evaluation-regions ${HC_BED}"
  echo "[$(date +%H:%M:%S)] HC regions BED: ${HC_BED}"
elif [[ -n "${HC_BED}" ]]; then
  echo "WARNING: HC_BED set but not found (${HC_BED}) — running without --evaluation-regions."
  echo "         Calls outside the truth's HC regions will inflate false-positive counts."
else
  echo "WARNING: HC_BED not set — running without --evaluation-regions."
  echo "         Calls outside the truth's HC regions will inflate false-positive counts."
fi

# ============================================================
# SNV evaluation
# ============================================================
SNV_EVAL="${EVAL_DIR}/eval_snv"
if [[ -d "${SNV_EVAL}" ]]; then
  echo "[$(date +%H:%M:%S)] Skipping SNV vcfeval — output dir exists: ${SNV_EVAL}"
else
  echo "[$(date +%H:%M:%S)] Running SNV vcfeval..."
  docker run --rm -v /mnt/storage:/mnt/storage ${RTG_IMAGE} vcfeval \
    --baseline "${SNV_TRUTH}" \
    --calls "${DS_PASS}" \
    --template "${SDF}" \
    --output "${SNV_EVAL}" \
    ${HC_BED_FLAG} \
    --vcf-score-field GQ \
    --squash-ploidy \
    --sample "TRUTH,${TUMOR_SAMPLE}"
  echo "[$(date +%H:%M:%S)] SNV vcfeval done → ${SNV_EVAL}"
fi

# ============================================================
# INDEL evaluation (INDEL-only calls to avoid mixed-VCF FP inflation)
# ============================================================
INDEL_EVAL="${EVAL_DIR}/eval_indel"
if [[ -d "${INDEL_EVAL}" ]]; then
  echo "[$(date +%H:%M:%S)] Skipping INDEL vcfeval — output dir exists: ${INDEL_EVAL}"
else
  echo "[$(date +%H:%M:%S)] Running INDEL vcfeval..."
  docker run --rm -v /mnt/storage:/mnt/storage ${RTG_IMAGE} vcfeval \
    --baseline "${INDEL_TRUTH}" \
    --calls "${DS_INDELS}" \
    --template "${SDF}" \
    --output "${INDEL_EVAL}" \
    ${HC_BED_FLAG} \
    --vcf-score-field GQ \
    --squash-ploidy \
    --sample "TRUTH,${TUMOR_SAMPLE}"
  echo "[$(date +%H:%M:%S)] INDEL vcfeval done → ${INDEL_EVAL}"
fi

# ============================================================
# Parse vcfeval summary.txt → extract best F1 row
# ============================================================
parse_vcfeval_summary() {
  local summary_file="$1"
  # summary.txt columns: Threshold | True-pos-baseline | True-pos-call | False-pos | False-neg | Precision | Sensitivity | F-measure
  # RTG marks best-F1 row with leading * in full ROC output; when only one threshold exists
  # the row is not marked. Prefer any non-None row over the catch-all None row; among
  # non-None rows pick the highest F-measure (col 8).
  local best
  best=$(grep '\*' "${summary_file}" 2>/dev/null | head -1)
  if [[ -z "${best}" ]]; then
    # No * marker — sort non-None rows by F-measure (col 8) descending
    best=$(awk 'NR>2 && $1!="None"' "${summary_file}" | sort -k8 -rn | head -1)
  fi
  # Fallback: if still empty (only None row present), use the None row
  [[ -z "${best}" ]] && best=$(tail -1 "${summary_file}")
  echo "${best}"
}

SNV_ROW=$(parse_vcfeval_summary "${SNV_EVAL}/summary.txt")
INDEL_ROW=$(parse_vcfeval_summary "${INDEL_EVAL}/summary.txt")

# vcfeval summary columns: Threshold | True-pos-baseline | True-pos-call | False-pos | False-neg | Precision | Sensitivity | F-measure
# Best-F1 row is prefixed with * — strip it for awk field counting
SNV_TP=$(echo "${SNV_ROW}"  | awk '{gsub(/^\*/,""); print $2}')
SNV_FP=$(echo "${SNV_ROW}"  | awk '{gsub(/^\*/,""); print $4}')
SNV_FN=$(echo "${SNV_ROW}"  | awk '{gsub(/^\*/,""); print $5}')
SNV_P=$(echo "${SNV_ROW}"   | awk '{gsub(/^\*/,""); print $6}')
SNV_Se=$(echo "${SNV_ROW}"  | awk '{gsub(/^\*/,""); print $7}')
SNV_F1=$(echo "${SNV_ROW}"  | awk '{gsub(/^\*/,""); print $8}')

INDEL_TP=$(echo "${INDEL_ROW}" | awk '{gsub(/^\*/,""); print $2}')
INDEL_FP=$(echo "${INDEL_ROW}" | awk '{gsub(/^\*/,""); print $4}')
INDEL_FN=$(echo "${INDEL_ROW}" | awk '{gsub(/^\*/,""); print $5}')
INDEL_P=$(echo "${INDEL_ROW}"  | awk '{gsub(/^\*/,""); print $6}')
INDEL_Se=$(echo "${INDEL_ROW}" | awk '{gsub(/^\*/,""); print $7}')
INDEL_F1=$(echo "${INDEL_ROW}" | awk '{gsub(/^\*/,""); print $8}')

# ============================================================
# Mutect2 best results (type-split) for comparison
# Hard-coded from validated benchmark run on HCC1395 (README, Jun 2026)
# Only shown in comparison table when running on HCC1395 (same sample).
# For other samples, Mutect2 was not run — rows are omitted.
# ============================================================
M2_SNV_TP=32886; M2_SNV_FP=10056; M2_SNV_FN=4512
M2_SNV_P="76.4%"; M2_SNV_Se="87.9%"; M2_SNV_F1=0.818
M2_INDEL_TP=2604; M2_INDEL_FP=2862; M2_INDEL_FN=758
M2_INDEL_P="47.6%"; M2_INDEL_Se="77.3%"; M2_INDEL_F1=0.590

# ============================================================
# Write comparison.md
# ============================================================
OUT_MD="${EVAL_DIR}/comparison.md"

# Build table content as a variable to avoid empty heredoc lines
TABLE_HEADER="| Caller | Variant type | TP | FP | FN | Precision | Sensitivity | F1 |
|--------|-------------|---:|---:|---:|----------:|------------:|---:|"
if [[ "${SAMPLE}" == "HCC1395" ]]; then
  TABLE_BODY="| Mutect2 (type-split best) | SNV   | ${M2_SNV_TP} | ${M2_SNV_FP} | ${M2_SNV_FN} | ${M2_SNV_P} | ${M2_SNV_Se} | ${M2_SNV_F1} |
| DeepSomatic (${MODEL_TYPE:-WGS})    | SNV   | ${SNV_TP} | ${SNV_FP} | ${SNV_FN} | ${SNV_P} | ${SNV_Se} | ${SNV_F1} |
| Mutect2 (type-split best) | INDEL | ${M2_INDEL_TP} | ${M2_INDEL_FP} | ${M2_INDEL_FN} | ${M2_INDEL_P} | ${M2_INDEL_Se} | ${M2_INDEL_F1} |
| DeepSomatic (${MODEL_TYPE:-WGS})    | INDEL | ${INDEL_TP} | ${INDEL_FP} | ${INDEL_FN} | ${INDEL_P} | ${INDEL_Se} | ${INDEL_F1} |"
  M2_NOTE="- Mutect2 result uses the type-split filter (SNV: MMQ≥50; INDEL: TLOD≥23.6) — the best achieved after a multi-step filter cascade (FilterMutectCalls → POPAF≥2 → MMQ≥50 / TLOD split)."
else
  TABLE_BODY="| DeepSomatic (${MODEL_TYPE:-WGS})    | SNV   | ${SNV_TP} | ${SNV_FP} | ${SNV_FN} | ${SNV_P} | ${SNV_Se} | ${SNV_F1} |
| DeepSomatic (${MODEL_TYPE:-WGS})    | INDEL | ${INDEL_TP} | ${INDEL_FP} | ${INDEL_FN} | ${INDEL_P} | ${INDEL_Se} | ${INDEL_F1} |"
  M2_NOTE="- Mutect2 was only benchmarked on HCC1395; no Mutect2 numbers available for ${SAMPLE}."
fi

cat > "${OUT_MD}" << EOF
# DeepSomatic vs Mutect2 — Somatic Variant Calling Benchmark

**Sample:** ${SAMPLE}
**DeepSomatic model:** ${MODEL_TYPE:-WGS}
**Truth set:** ${TRUTH_LABEL}
**Evaluated:** $(date '+%Y-%m-%d %H:%M')

---

## ⚠ CRITICAL CAVEAT

> **${SAMPLE} is in DeepSomatic's training data — this score measures memorization, not
> generalization. Do not cite these numbers as evidence of clinical performance.**
>
> HCC1395 and COLO829 are both confirmed in DeepSomatic's training corpus.
> These evaluations are presented for pipeline validation only.

---

## Results

${TABLE_HEADER}
${TABLE_BODY}

**Notes:**
${M2_NOTE}
- DeepSomatic result uses PASS-only FILTER, no post-processing filter cascade applied.
- INDEL evaluation uses INDEL-only calls VCF to avoid the mixed-VCF FP-inflation artifact
  where all SNVs in a mixed VCF are counted as INDEL FPs by vcfeval.
- vcfeval score field: \`GQ\` (DeepSomatic Genotype Quality) for ROC threshold sweep; results
  shown at the GQ threshold that maximises F-measure (reported as the threshold value in the
  vcfeval summary row). Numbers at Threshold=None (no GQ filter) are also available in
  \`eval_snv/summary.txt\` and \`eval_indel/summary.txt\`.

---

## Methodology

\`\`\`
Truth set:   ${TRUTH_LABEL}
HC regions:  ${HC_BED_FLAG:+$(basename "${HC_BED}") (evaluation restricted to HC regions)}${HC_BED_FLAG:-none — calls outside truth HC regions counted as FP}
Evaluator:   RTG vcfeval 3.12.1, --squash-ploidy${HC_BED_FLAG:+ --evaluation-regions}
Calls:       PASS-only VCF from DeepSomatic Parabricks 4.7.0-1
Reference:   GRCh38 no-alt (455 contigs)
\`\`\`

EOF

echo "[$(date +%H:%M:%S)] Comparison written: ${OUT_MD}"
echo ""
echo "========================================================"
echo "RESULTS SUMMARY"
echo "========================================================"
cat "${OUT_MD}"
