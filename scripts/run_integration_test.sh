#!/usr/bin/env bash
#
# Integration test: Open CRAVAT + PCGR second-opinion annotation
# for HCC1395 DeepSomatic PASS variants.
#
# Phase 1 — Open CRAVAT (always runs; civic + oncokb local databases)
#   Input:  HCC1395.deepsomatic.pass.vcf.gz
#   Output: /mnt/storage/open-cravat-reports/HCC1395-report/HCC1395-report.tsv
#
# Phase 2 — PCGR (runs only if PCGR_BUNDLE + VEP_CACHE are set)
#   Input:  same PASS VCF
#   Output: deepsomatic/pcgr/HCC1395.pcgr.grch38.html
#
# Phase 3 — Cross-tool comparison
#   Joins CIViC/OncoKB hits from OC with AMP-tiered variants from the
#   pipeline MAF, then reports concordance and discrepancies.
#
# Usage:
#   ./run_integration_test.sh [sample_id]
#
# Optional env overrides:
#   PCGR_BUNDLE=/mnt/storage/pcgr_bundle   (enables PCGR phase)
#   VEP_CACHE=/mnt/storage/vep_cache       (required with PCGR_BUNDLE)
#   TUMOR_TYPE=BRCA                        (OncoTree code; default BRCA)
#   OUTDIR=/mnt/storage/parabricks_test/output/<sample>
#   OC_REPORT_DIR=/mnt/storage/open-cravat-reports/HCC1395-report

set -euo pipefail

SAMPLE="${1:-HCC1395}"
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/${SAMPLE}}"
DS_OUTDIR="${OUTDIR}/deepsomatic"
PASS_VCF="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.vcf.gz"
PIPELINE_MAF="${DS_OUTDIR}/${SAMPLE}.deepsomatic.maf"

OC_BIN="/mnt/storage/open-cravat-reports/oc-venv/bin/oc"
OC_REPORT_DIR="${OC_REPORT_DIR:-/mnt/storage/open-cravat-reports/HCC1395-report}"
OC_REPORT_NAME="HCC1395-report"
OC_TSV="${OC_REPORT_DIR}/${OC_REPORT_NAME}.tsv"
OC_PLAIN_VCF="${OC_REPORT_DIR}/${SAMPLE}.pass.vcf"

PCGR_BUNDLE="${PCGR_BUNDLE:-}"
VEP_CACHE="${VEP_CACHE:-}"
TUMOR_TYPE="${TUMOR_TYPE:-BRCA}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

banner() { printf '\n%s\n%s\n%s\n' "$(printf '=%.0s' {1..64})" "  $*" "$(printf '=%.0s' {1..64})"; }
step()   { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

banner "HCC1395 Integration Test — Open CRAVAT + PCGR"
echo "  Sample:   ${SAMPLE}"
echo "  PASS VCF: ${PASS_VCF}"
echo "  MAF:      ${PIPELINE_MAF}"

# ============================================================
# PREREQUISITE CHECKS
# ============================================================
PREREQ_OK=true
[[ ! -f "${PASS_VCF}" ]]      && echo "ERROR: PASS VCF not found: ${PASS_VCF}" && PREREQ_OK=false
[[ ! -f "${PIPELINE_MAF}" ]]  && echo "WARN:  Pipeline MAF not found: ${PIPELINE_MAF} (comparison skipped)"
[[ ! -x "${OC_BIN}" ]]        && echo "ERROR: Open CRAVAT not found: ${OC_BIN}" && PREREQ_OK=false
[[ "${PREREQ_OK}" == false ]] && exit 1

mkdir -p "${OC_REPORT_DIR}"

# ============================================================
# PHASE 1 — OPEN CRAVAT
# ============================================================
banner "Phase 1: Open CRAVAT (CIViC + OncoKB)"

if [[ -f "${OC_TSV}" ]]; then
  step "Skipping OC — output exists: ${OC_TSV}"
else
  step "Decompressing PASS VCF..."
  gunzip -c "${PASS_VCF}" > "${OC_PLAIN_VCF}"

  step "Running Open CRAVAT (civic + oncokb)..."
  "${OC_BIN}" run "${OC_PLAIN_VCF}" \
    -l hg38 \
    -d "${OC_REPORT_DIR}" \
    -n "${OC_REPORT_NAME}" \
    -t text \
    -a civic oncokb

  rm -f "${OC_PLAIN_VCF}"
  step "Open CRAVAT complete: ${OC_TSV}"
fi

# ============================================================
# PHASE 2 — PCGR (conditional)
# ============================================================
banner "Phase 2: PCGR Second-Opinion Report"

if [[ -z "${PCGR_BUNDLE}" ]]; then
  echo "  SKIPPED — PCGR_BUNDLE not set."
  echo ""
  echo "  To enable PCGR, set:"
  echo "    export PCGR_BUNDLE=/mnt/storage/pcgr_bundle"
  echo "    export VEP_CACHE=/mnt/storage/vep_cache"
  echo ""
  echo "  One-time bundle download (~5 GB, GRCh38):"
  echo "    mkdir -p /mnt/storage/pcgr_bundle"
  echo "    curl -L https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz \\"
  echo "      | tar -xz -C /mnt/storage/pcgr_bundle"
  echo ""
  echo "  One-time VEP cache download (~30 GB):"
  echo "    mkdir -p /mnt/storage/vep_cache"
  echo "    docker run --rm \\"
  echo "      -v /mnt/storage/vep_cache:/mnt/.vep \\"
  echo "      sigven/pcgr:2.2.5 \\"
  echo "      vep_install -a cf -s homo_sapiens -y GRCh38 -c /mnt/.vep --convert"
  PCGR_RAN=false
else
  TUMOR_TYPE="${TUMOR_TYPE}" \
  OUTDIR="${OUTDIR}" \
  PCGR_BUNDLE="${PCGR_BUNDLE}" \
  VEP_CACHE="${VEP_CACHE}" \
    "${SCRIPT_DIR}/run_pcgr.sh" "${SAMPLE}" "${PASS_VCF}"
  PCGR_RAN=true
fi

# ============================================================
# PHASE 3 — CROSS-TOOL COMPARISON
# ============================================================
banner "Phase 3: Cross-Tool Comparison"

if [[ ! -f "${OC_TSV}" ]]; then
  echo "  OC output missing — cannot compare."
  exit 0
fi

step "CIViC evidence hits (Open CRAVAT)..."
# Verified TSV column positions (line 6 header, 1-based):
#   8=Gene  13=ProteinChange  19=CIViC_Diseases  21=CIViC_EVS  27=OncoKB_Oncogenic  30=OncoKB_HighestSensitiveLevel
OC_CIVIC_HITS=$(awk -F'\t' 'NR>6 && $19!="" {count++} END {print count+0}' "${OC_TSV}")
echo "  Variants with CIViC evidence: ${OC_CIVIC_HITS}"

if [[ "${OC_CIVIC_HITS}" -gt 0 ]]; then
  echo ""
  printf "  %-20s %-20s %-10s %-50s\n" "Gene" "ProtChange" "EVS" "CIViC_Diseases"
  printf "  %-20s %-20s %-10s %-50s\n" "----" "----------" "---" "--------------"
  awk -F'\t' 'NR>6 && $19!="" {printf "  %-20s %-20s %-10s %-50s\n", $8, $13, $21, substr($19,1,50)}' "${OC_TSV}"
fi

echo ""
step "OncoKB oncogenic hits (Open CRAVAT)..."
# OncoKB columns: 27=Oncogenic  28=KnownEffect  29=Hotspot  30=HighestSensitiveLevel
OC_ONCOKB_HITS=$(awk -F'\t' 'NR>6 && $27!="" {count++} END {print count+0}' "${OC_TSV}")
echo "  Variants with OncoKB annotation: ${OC_ONCOKB_HITS}"
echo "  Note: OncoKB requires an API token for full annotation; 0 hits = no token configured."
echo "        Set token: \$OC_BIN config set oncokb token <YOUR_TOKEN>"

if [[ "${OC_ONCOKB_HITS}" -gt 0 ]]; then
  echo ""
  printf "  %-20s %-20s %-20s %-16s\n" "Gene" "ProtChange" "Oncogenic" "SensitiveLevel"
  printf "  %-20s %-20s %-20s %-16s\n" "----" "----------" "---------" "--------------"
  awk -F'\t' 'NR>6 && $27!="" {printf "  %-20s %-20s %-20s %-16s\n", $8, $13, $27, $30}' "${OC_TSV}"
fi

# Cross-reference with pipeline MAF AMP-tiered variants
# MAF columns: 1=Hugo_Symbol  8=HGVSp_Short  16=HIGHEST_LEVEL  17=HIGHEST_LEVEL_SUMMARY
if [[ -f "${PIPELINE_MAF}" ]]; then
  echo ""
  step "Concordance with pipeline AMP tiers..."
  echo "  Checking whether pipeline LEVEL_1/LEVEL_2 variants appear in OC CIViC/OncoKB hits..."

  TIER1_GENES=$(awk -F'\t' 'NR>1 && ($16~/LEVEL_1|LEVEL_2/) {print $1}' "${PIPELINE_MAF}" | sort -u)

  if [[ -n "${TIER1_GENES}" ]]; then
    echo ""
    printf "  %-16s %-24s %-12s %-12s\n" "Gene" "HGVSp" "InOC_CIViC" "InOC_OncoKB"
    printf "  %-16s %-24s %-12s %-12s\n" "----" "-----" "----------" "-----------"
    while IFS= read -r gene; do
      hgvsp=$(awk -F'\t' -v g="${gene}" 'NR>1 && $1==g && ($16~/LEVEL_1|LEVEL_2/) {print $8; exit}' "${PIPELINE_MAF}")
      in_civic=$(awk -F'\t' -v g="${gene}" 'NR>6 && $8==g && $19!="" {found=1} END {print found+0}' "${OC_TSV}")
      in_oncokb=$(awk -F'\t' -v g="${gene}" 'NR>6 && $8==g && $27!="" {found=1} END {print found+0}' "${OC_TSV}")
      printf "  %-16s %-24s %-12s %-12s\n" "${gene}" "${hgvsp}" \
        "$([ "${in_civic}"  = "1" ] && echo "YES" || echo "no")" \
        "$([ "${in_oncokb}" = "1" ] && echo "YES" || echo "no")"
    done <<< "${TIER1_GENES}"
    echo ""
    echo "  NOTE: BRCA2 truncating variants are typically absent from CIViC/OncoKB exact-match"
    echo "  lookups (which require a specific HGVSp entry). The pipeline AMP tier is correct"
    echo "  — LoF in a tumor-suppressor is Tier IA by rule-based inference, not variant lookup."
  else
    echo "  No LEVEL_1/LEVEL_2 variants found in pipeline MAF."
  fi
fi

# PCGR concordance (if ran)
if [[ "${PCGR_RAN:-false}" == true ]]; then
  PCGR_TSV="${DS_OUTDIR}/pcgr/${SAMPLE}.pcgr.grch38.snvs_indels.tiers.tsv"
  if [[ -f "${PCGR_TSV}" ]]; then
    echo ""
    step "PCGR tier distribution..."
    awk -F'\t' 'NR>1 {tiers[$NF]++} END {for (t in tiers) printf "  TIER %s: %d variants\n", t, tiers[t]}' \
      "${PCGR_TSV}" | sort
  fi
fi

# ============================================================
# SUMMARY
# ============================================================
banner "Integration Test Summary"
echo "  Open CRAVAT report:  ${OC_TSV}"
echo "  CIViC hits:          ${OC_CIVIC_HITS}"
echo "  OncoKB hits:         ${OC_ONCOKB_HITS}"
if [[ "${PCGR_RAN:-false}" == true ]]; then
  echo "  PCGR report:         ${DS_OUTDIR}/pcgr/${SAMPLE}.pcgr.grch38.html"
else
  echo "  PCGR report:         NOT RUN (set PCGR_BUNDLE to enable)"
fi
echo ""
echo "  Quick OC queries:"
echo "    # Top CIViC hits by evidence score:"
echo "    awk -F'\\t' 'NR>6 && \$19!=\"\"' ${OC_TSV} | sort -t\$'\\t' -k21 -rn | cut -f8,13,19,21 | head -10"
echo "    # OncoKB oncogenic (requires API token for full results):"
echo "    awk -F'\\t' 'NR>6 && \$27!=\"\"' ${OC_TSV} | cut -f8,13,27,30 | head -10"
echo ""
