#!/usr/bin/env bash
#
# Download and prepare COLO829 somatic truth VCFs for DeepSomatic held-out validation
#
# Source: Alioto et al. 2015 (Nature Communications) / ICGC
# These truth calls are used in the DeepSomatic paper (Chen et al. 2024) as an
# independent benchmark. NOTE: Some COLO829 data may overlap with DeepSomatic's
# training set — see DeepSomatic supplementary materials for training data provenance.
#
# Truth set used: Boutros lab COLO829 benchmark
# (https://github.com/PapenfussLab/COLO829-benchmarking)
#
# Prep steps (mirrors HCC1395 truth prep):
#   1. Download VCF from GitHub / Zenodo
#   2. Normalize FILTER field
#   3. Add synthetic GT column (site-only → FORMAT GT=1/1)
#   4. bgzip + tabix index
#
# Usage:
#   ./download_colo829_truth.sh [--snv-only | --indel-only]
#
# Output:
#   /mnt/storage/parabricks_test/data/truth/colo829/
#     colo829_sSNV_gt.vcf.gz
#     colo829_sINDEL_gt.vcf.gz

set -euo pipefail

TRUTHDIR="${TRUTHDIR:-/mnt/storage/parabricks_test/data/truth/colo829}"
SNV_ONLY=false
INDEL_ONLY=false

for arg in "$@"; do
  case "$arg" in
    --snv-only)   SNV_ONLY=true ;;
    --indel-only) INDEL_ONLY=true ;;
  esac
done

mkdir -p "${TRUTHDIR}"

# ============================================================
# Source URLs
# Boutros lab COLO829 benchmark (Zenodo)
# Primary paper: https://doi.org/10.1038/s41467-020-16182-z
# ============================================================
ZENODO_BASE="https://zenodo.org/record/3979165/files"
SNV_URL="${ZENODO_BASE}/COLO829.somatic.snvs.vcf.gz"
INDEL_URL="${ZENODO_BASE}/COLO829.somatic.indels.vcf.gz"

BCFTOOLS_IMAGE="staphb/bcftools:1.19"
BCFTOOLS_RUN="docker run --rm -v /mnt/storage:/mnt/storage ${BCFTOOLS_IMAGE} bcftools"

download_and_prep() {
  local url="$1"
  local raw_out="$2"
  local prep_out="$3"
  local label="$4"

  if [[ -f "${prep_out}" && -f "${prep_out}.tbi" ]]; then
    echo "[$(date +%H:%M:%S)] Skipping ${label} — already prepared: ${prep_out}"
    return
  fi

  # Download raw VCF
  if [[ ! -f "${raw_out}" ]]; then
    echo "[$(date +%H:%M:%S)] Downloading ${label} truth VCF..."
    wget -q --show-progress -O "${raw_out}" "${url}" || \
      curl -fL -o "${raw_out}" "${url}"
    # Try tabix in case it's already bgzipped but lacks index
    tabix -p vcf "${raw_out}" 2>/dev/null || true
  fi

  # Prep: normalize FILTER + add synthetic GT column
  echo "[$(date +%H:%M:%S)] Preparing ${label} truth VCF..."
  ${BCFTOOLS_RUN} view "${raw_out}" | \
    awk 'BEGIN{OFS="\t"}
      /^##/{print;next}
      /^#CHROM/{print $0"\tFORMAT\tTRUTH";next}
      {gsub(/PASS;HighConf/,"PASS",$7); print $0"\tGT\t1/1"}
    ' | \
    ${BCFTOOLS_RUN} view -Oz -o "${prep_out}"

  ${BCFTOOLS_RUN} index -t "${prep_out}"
  echo "[$(date +%H:%M:%S)] ${label} ready: ${prep_out}"
}

[[ "${INDEL_ONLY}" == false ]] && \
  download_and_prep "${SNV_URL}" \
    "${TRUTHDIR}/COLO829.somatic.snvs.raw.vcf.gz" \
    "${TRUTHDIR}/colo829_sSNV_gt.vcf.gz" \
    "SNV"

[[ "${SNV_ONLY}" == false ]] && \
  download_and_prep "${INDEL_URL}" \
    "${TRUTHDIR}/COLO829.somatic.indels.raw.vcf.gz" \
    "${TRUTHDIR}/colo829_sINDEL_gt.vcf.gz" \
    "INDEL"

echo ""
echo "COLO829 truth VCFs ready in ${TRUTHDIR}/"
echo "To evaluate DeepSomatic on COLO829, set:"
echo "  TRUTHDIR=${TRUTHDIR} SAMPLE=COLO829 ./eval_deepsomatic.sh COLO829"
echo ""
echo "NOTE: Verify COLO829 is not in DeepSomatic training data before citing"
echo "      these numbers as a generalization result. Check:"
echo "      https://github.com/google/deepsomatic — training data section"
