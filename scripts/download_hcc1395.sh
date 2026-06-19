#!/usr/bin/env bash
#
# Download HCC1395 / HCC1395BL SEQC2 benchmark data from NCBI SRA.
#
# Source: PRJNA489865 (FDA SEQC2 somatic mutation benchmarking)
# Sample: HCC1395 (breast cancer tumor) + HCC1395BL (matched B-lymphoblastoid normal)
# Platform: Illumina HiSeq X Ten, 150 bp PE, ~50–55x WGS coverage
#
# Using the FD-lab runs (Foundation Diagnostics) for a matched tumor-normal pair:
#   Tumor  SRR7890829  WGS_FD_T_1  ~52x
#   Tumor  SRR7890824  WGS_FD_T_2  ~53x   (combined ≈ 105x — use one or both)
#   Normal SRR7890826  WGS_FD_N_1  ~50x
#
# Tools required: prefetch + fasterq-dump (SRA toolkit) — both already installed.
#
# Estimated download: ~50 GB/run (compressed .sra); conversion to FASTQ ~130 GB each.
# Plan for ~700 GB total scratch space for a tumor-normal pair.
#
# Usage:
#   ./download_hcc1395.sh [--tumor-only | --normal-only | --truth-only]

set -euo pipefail

DATA_DIR="${DATA_DIR:-/mnt/storage/parabricks_test/data/hcc1395}"
TRUTH_DIR="${TRUTH_DIR:-/mnt/storage/parabricks_test/data/truth}"
PREFETCH_OPTS="--max-size 200GB --progress"

MODE="${1:-all}"

mkdir -p "${DATA_DIR}/tumor" "${DATA_DIR}/normal" "${TRUTH_DIR}"

# ============================================================
# Truth VCF (SEQC2 high-confidence sSNVs and sINDELs)
# ============================================================
download_truth() {
  echo "[truth] Downloading SEQC2 truth VCFs..."
  BASE="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest"
  cd "${TRUTH_DIR}"
  for f in \
    "high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz" \
    "high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz" \
    "high-confidence_sINDEL_in_HC_regions_v1.2.vcf.gz" \
    "High-Confidence_Regions_v1.2.bed" \
    "High-Confidence_Regions_v1.2.1.bed"; do
    [[ -f "$f" ]] && echo "[skip] $f" && continue
    wget -q --show-progress "${BASE}/${f}" 2>/dev/null \
      || echo "[warn] ${f} not found on server (may not exist for this release)"
  done
  # Index VCFs for rtg vcfeval
  for vcf in *.vcf.gz; do
    [[ -f "${vcf}.tbi" ]] || tabix -p vcf "${vcf}" 2>/dev/null || true
  done
  echo "[truth] Done."
}

# ============================================================
# FASTQ download via SRA toolkit
# ============================================================
sra_download() {
  local accession="$1" outdir="$2" label="$3"
  echo "[${label}] Fetching ${accession} with prefetch..."
  prefetch ${PREFETCH_OPTS} --output-directory "${outdir}" "${accession}"

  echo "[${label}] Converting to FASTQ..."
  fasterq-dump \
    --outdir "${outdir}" \
    --threads 8 \
    --split-files \
    --progress \
    "${outdir}/${accession}/${accession}.sra" \
    2>&1 | tail -5

  echo "[${label}] Compressing FASTQ..."
  pigz -p 8 "${outdir}/${accession}_1.fastq" "${outdir}/${accession}_2.fastq" 2>/dev/null \
    || gzip "${outdir}/${accession}_1.fastq" "${outdir}/${accession}_2.fastq"

  rm -rf "${outdir}/${accession}"
  echo "[${label}] Done → ${outdir}/${accession}_{1,2}.fastq.gz"
}

# ============================================================
# Main
# ============================================================
[[ "${MODE}" == "--truth-only" ]] && download_truth && exit 0

if [[ "${MODE}" == "all" || "${MODE}" == "--tumor-only" ]]; then
  echo "=== Tumor: SRR7890829 (HCC1395, WGS_FD_T_1, ~52x) ==="
  if [[ ! -f "${DATA_DIR}/tumor/SRR7890829_1.fastq.gz" ]]; then
    sra_download SRR7890829 "${DATA_DIR}/tumor" "tumor"
  else
    echo "[skip] SRR7890829 already present"
  fi
fi

if [[ "${MODE}" == "all" || "${MODE}" == "--normal-only" ]]; then
  echo "=== Normal: SRR7890826 (HCC1395BL, WGS_FD_N_1, ~50x) ==="
  if [[ ! -f "${DATA_DIR}/normal/SRR7890826_1.fastq.gz" ]]; then
    sra_download SRR7890826 "${DATA_DIR}/normal" "normal"
  else
    echo "[skip] SRR7890826 already present"
  fi
fi

[[ "${MODE}" == "all" ]] && download_truth

echo ""
echo "=== HCC1395 data ready ==="
ls -lh "${DATA_DIR}/tumor/" "${DATA_DIR}/normal/" "${TRUTH_DIR}/"
echo ""
echo "Run the pipeline with:"
echo "  NUM_GPUS=2 GPU_DEVICE=0,1 \\"
echo "  REF_DIR=/mnt/storage/parabricks_test/ref \\"
echo "  ./scripts/run_parabricks_pipeline.sh HCC1395 \\"
echo "    ${DATA_DIR}/tumor/SRR7890829_1.fastq.gz \\"
echo "    ${DATA_DIR}/tumor/SRR7890829_2.fastq.gz \\"
echo "    ${DATA_DIR}/normal/SRR7890826_1.fastq.gz \\"
echo "    ${DATA_DIR}/normal/SRR7890826_2.fastq.gz"
