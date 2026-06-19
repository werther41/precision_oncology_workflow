#!/usr/bin/env bash
#
# DeepSomatic GPU somatic variant calling pipeline
# Reuses existing aligned BAMs — do NOT re-run fq2bam.
#
# Phases 0-2:
#   0. Preflight: BAM/index presence, SM-tag distinctness, GPU check
#   1. pbrun deepsomatic call (GPU-accelerated)
#   2. PASS extraction + INDEL-only split; sample-column auto-detection
#
# Usage:
#   ./run_deepsomatic_pipeline.sh <sample_id> <tumor_bam> <normal_bam>
#
# Optional env overrides:
#   NUM_GPUS=2           number of GPUs for DeepSomatic (default 2 — uses both GPUs)
#   GPU_DEVICE=all       GPU device(s) exposed to docker (default "all"; use "0" or "1" to restrict)
#   REF_DIR=/path/ref    reference bundle directory
#   OUTDIR=/path/out     output root (deepsomatic/ subdir created here)
#   MODEL_TYPE=WGS       DeepSomatic model: WGS (default) | WES | pacbio | ont
#                        WGS → --mode shortread (default model)
#                        WES → --mode shortread --use-wes-model
#                        pacbio/ont → --mode pacbio / --mode ont
#
# Flags:
#   --preflight-only     run Phase 0 only, then exit

set -euo pipefail

# ============================================================
# PARSE ARGUMENTS
# ============================================================
PREFLIGHT_ONLY=false
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --preflight-only) PREFLIGHT_ONLY=true ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]:-}"

SAMPLE="${1:-HCC1395}"
TUMOR_BAM="${2:-}"
NORMAL_BAM="${3:-}"

NUM_GPUS="${NUM_GPUS:-2}"
GPU_DEVICE="${GPU_DEVICE:-all}"
MODEL_TYPE="${MODEL_TYPE:-WGS}"

REF_DIR="${REF_DIR:-/mnt/storage/parabricks_test/ref}"
REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/${SAMPLE}}"
DS_OUTDIR="${OUTDIR}/deepsomatic"

TUMOR_BAM="${TUMOR_BAM:-${OUTDIR}/${SAMPLE}_tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:-${OUTDIR}/${SAMPLE}_normal.bam}"

PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
BCFTOOLS_IMAGE="staphb/bcftools:1.19"

# ============================================================
# PHASE 0: PREFLIGHT
# ============================================================
echo "========================================================"
echo "[$(date +%H:%M:%S)] PHASE 0: Preflight"
echo "========================================================"

# BAMs and indices
PREFLIGHT_OK=true
for bam in "${TUMOR_BAM}" "${NORMAL_BAM}"; do
  if [[ ! -f "${bam}" ]]; then
    echo "ERROR: BAM not found: ${bam}"
    PREFLIGHT_OK=false
  elif [[ ! -f "${bam}.bai" ]]; then
    echo "ERROR: BAM index not found: ${bam}.bai"
    PREFLIGHT_OK=false
  else
    SIZE=$(du -sh "${bam}" | cut -f1)
    echo "  OK  ${bam} (${SIZE})"
  fi
done

# Reference
if [[ ! -f "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA not found: ${REF_FASTA}"
  PREFLIGHT_OK=false
elif [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "ERROR: Reference FASTA index not found: ${REF_FASTA}.fai"
  PREFLIGHT_OK=false
else
  echo "  OK  ${REF_FASTA}"
fi

[[ "${PREFLIGHT_OK}" == false ]] && { echo "Preflight FAILED"; exit 1; }

# SM tag distinctness — DeepSomatic requires tumor != normal
TUMOR_SM=$(samtools view -H "${TUMOR_BAM}" | grep "^@RG" | grep -oP "SM:\K[^\t]+" | sort -u | tr '\n' ',' | sed 's/,$//')
NORMAL_SM=$(samtools view -H "${NORMAL_BAM}" | grep "^@RG" | grep -oP "SM:\K[^\t]+" | sort -u | tr '\n' ',' | sed 's/,$//')
echo "  Tumor  SM tag: ${TUMOR_SM}"
echo "  Normal SM tag: ${NORMAL_SM}"
if [[ "${TUMOR_SM}" == "${NORMAL_SM}" ]]; then
  echo "ERROR: Tumor and normal share the same SM read-group tag (${TUMOR_SM})."
  echo "       DeepSomatic requires distinct SM tags to differentiate tumor from normal."
  exit 1
fi

# GPU check
echo ""
echo "  GPU status:"
nvidia-smi --query-gpu=index,name,memory.free,memory.total --format=csv,noheader | \
  awk -F', ' '{printf "    GPU %s: %s  free=%s / total=%s\n", $1, $2, $3, $4}'
ACTIVE_PIDS=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader | wc -l)
if [[ "${ACTIVE_PIDS}" -gt 0 ]]; then
  echo "  WARNING: ${ACTIVE_PIDS} compute process(es) currently using GPUs."
  echo "           GPU_DEVICE=${GPU_DEVICE}; confirm required GPU(s) are free."
fi

echo ""
echo "  Sample:     ${SAMPLE}"
echo "  Model type: ${MODEL_TYPE}"
echo "  Tumor BAM:  ${TUMOR_BAM}"
echo "  Normal BAM: ${NORMAL_BAM}"
echo "  Output dir: ${DS_OUTDIR}"
echo "  GPUs:       ${NUM_GPUS} (device: ${GPU_DEVICE})"
echo ""
echo "Preflight PASSED"

[[ "${PREFLIGHT_ONLY}" == true ]] && exit 0

# ============================================================
# PHASE 1: DeepSomatic call
# ============================================================
echo "========================================================"
echo "[$(date +%H:%M:%S)] PHASE 1: DeepSomatic variant calling"
echo "========================================================"

mkdir -p "${DS_OUTDIR}"
DS_VCF="${DS_OUTDIR}/${SAMPLE}.deepsomatic.vcf.gz"

if [[ -f "${DS_VCF}" && -f "${DS_VCF}.tbi" ]]; then
  echo "[$(date +%H:%M:%S)] Skipping DeepSomatic call — output exists: ${DS_VCF}"
else
  DOCKER_RUN="docker run --gpus \"device=${GPU_DEVICE}\" --rm \
    --user $(id -u):$(id -g) \
    -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
    -v /mnt/storage:/mnt/storage \
    -v /tmp:/tmp \
    -e TMPDIR=/tmp \
    ${PB_IMAGE}"

  TIME_START=$(date +%s)
  echo "[$(date +%H:%M:%S)] Starting pbrun deepsomatic..."

  # Translate MODEL_TYPE to pbrun flags
  # --mode: shortread (Illumina WGS/WES), pacbio, ont
  # WES shortread uses --use-wes-model; WGS shortread uses default model
  case "${MODEL_TYPE}" in
    WGS)    PB_MODE="shortread"; PB_MODEL_EXTRA="" ;;
    WES)    PB_MODE="shortread"; PB_MODEL_EXTRA="--use-wes-model" ;;
    pacbio) PB_MODE="pacbio";    PB_MODEL_EXTRA="" ;;
    ont)    PB_MODE="ont";       PB_MODEL_EXTRA="" ;;
    *)      echo "ERROR: Unknown MODEL_TYPE '${MODEL_TYPE}'. Use WGS, WES, pacbio, or ont."; exit 1 ;;
  esac

  ${DOCKER_RUN} pbrun deepsomatic \
    --ref "${REF_FASTA}" \
    --in-tumor-bam "${TUMOR_BAM}" \
    --in-normal-bam "${NORMAL_BAM}" \
    --out-variants "${DS_VCF}" \
    --mode "${PB_MODE}" \
    ${PB_MODEL_EXTRA} \
    --num-gpus "${NUM_GPUS}" \
    --tmp-dir /tmp

  ELAPSED=$(( $(date +%s) - TIME_START ))
  echo "[$(date +%H:%M:%S)] DeepSomatic completed in ${ELAPSED}s ($(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s)"
  echo "  Runtime | DeepSomatic (${MODEL_TYPE}) | ${ELAPSED}s" >> "${DS_OUTDIR}/runtime.log"

  # Ensure tabix index exists
  if [[ ! -f "${DS_VCF}.tbi" ]]; then
    docker run --rm -v /mnt/storage:/mnt/storage ${BCFTOOLS_IMAGE} tabix -p vcf "${DS_VCF}"
  fi
fi

# ============================================================
# PHASE 2: PASS extraction + sample-column detection
# ============================================================
echo "========================================================"
echo "[$(date +%H:%M:%S)] PHASE 2: PASS extraction + sample detection"
echo "========================================================"

# Auto-detect sample structure
N_SAMPLES=$(docker run --rm -v /mnt/storage:/mnt/storage ${BCFTOOLS_IMAGE} \
  bcftools view -h "${DS_VCF}" | grep "^#CHROM" | awk '{print NF-9}')
TUMOR_SAMPLE=$(docker run --rm -v /mnt/storage:/mnt/storage ${BCFTOOLS_IMAGE} \
  bcftools view -h "${DS_VCF}" | grep "^#CHROM" | awk '{print $NF}')
echo "  DeepSomatic VCF: ${N_SAMPLES} sample column(s), tumor sample name = '${TUMOR_SAMPLE}'"

# Verify SM tag matches VCF sample name (critical for vcfeval --sample flag)
if [[ "${TUMOR_SAMPLE}" != "${TUMOR_SM}" ]]; then
  echo "  WARNING: VCF sample name '${TUMOR_SAMPLE}' differs from BAM SM tag '${TUMOR_SM}'."
  echo "           eval_deepsomatic.sh will use VCF sample name '${TUMOR_SAMPLE}' for vcfeval."
fi

# Export for eval script to consume
echo "${TUMOR_SAMPLE}" > "${DS_OUTDIR}/${SAMPLE}.sample_name.txt"
echo "${N_SAMPLES}"    > "${DS_OUTDIR}/${SAMPLE}.n_samples.txt"

DS_PASS="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.vcf.gz"
DS_INDELS="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.indels.vcf.gz"

BCFTOOLS_RUN="docker run --rm -v /mnt/storage:/mnt/storage ${BCFTOOLS_IMAGE} bcftools"

# PASS-only VCF
if [[ ! -f "${DS_PASS}" ]]; then
  ${BCFTOOLS_RUN} view -f PASS "${DS_VCF}" -Oz -o "${DS_PASS}"
  ${BCFTOOLS_RUN} index -t "${DS_PASS}"
  echo "[$(date +%H:%M:%S)] PASS VCF written: ${DS_PASS}"
else
  echo "[$(date +%H:%M:%S)] Skipping PASS extract — exists: ${DS_PASS}"
fi

# INDEL-only subset for separate INDEL evaluation (avoids mixed-VCF FP inflation in vcfeval)
if [[ ! -f "${DS_INDELS}" ]]; then
  ${BCFTOOLS_RUN} view -f PASS --type indels "${DS_VCF}" -Oz -o "${DS_INDELS}"
  ${BCFTOOLS_RUN} index -t "${DS_INDELS}"
  echo "[$(date +%H:%M:%S)] INDEL-only VCF written: ${DS_INDELS}"
else
  echo "[$(date +%H:%M:%S)] Skipping INDEL extract — exists: ${DS_INDELS}"
fi

# FILTER distribution
echo ""
echo "  FILTER value distribution:"
${BCFTOOLS_RUN} query -f '%FILTER\n' "${DS_VCF}" | sort | uniq -c | sort -rn | head -10

# Variant count
PASS_COUNT=$(${BCFTOOLS_RUN} view -f PASS "${DS_VCF}" | grep -v "^#" | wc -l)
echo ""
echo "  PASS variants: ${PASS_COUNT}"
echo ""
echo "[$(date +%H:%M:%S)] Phases 0-2 complete. Run eval_deepsomatic.sh to benchmark."
