#!/usr/bin/env bash
#
# Smoke test: TruSeq Amplicon Cancer Panel data, tumor-normal paired
#   Tumor  = SIP0001  (onc089eec2c-SIP0001_S1_L001_R{1,2}_001.fastq.gz)
#   Normal = SIP0002  (onc089eec2c-SIP0002_S2_L001_R{1,2}_001.fastq.gz)
#
# Hardware assumptions (this machine):
#   GPU 0 — RTX PRO 6000 Blackwell, ~96 GB, free
#   GPU 1 — same, occupied by VLLM — DO NOT USE
#
# Missing Broad bundle files (BQSR, PoN, gnomAD) are detected at runtime;
# the pipeline degrades gracefully: alignment runs without BQSR,
# Mutect2 runs without PoN/germline-resource (lower specificity, fine for smoke test).
#
# Expected runtime: ~3–5 min total on a single RTX PRO 6000

set -euo pipefail

SAMPLE_ID="${1:-SMOKETEST_TN}"

PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
GPU_DEVICE="0"
NUM_GPUS="1"

DATA_DIR="/mnt/storage/parabricks_test/data/pillar_test_data"
REF_DIR="/mnt/storage/parabricks_test/ref"
OUTDIR="/mnt/storage/parabricks_test/output/smoketest/${SAMPLE_ID}"
mkdir -p "${OUTDIR}"

REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"
KNOWN_SITES="${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
PON="${REF_DIR}/1000g_pon.hg38.vcf.gz"
GERMLINE_RESOURCE="${REF_DIR}/af-only-gnomad.hg38.vcf.gz"

TUMOR_R1="${DATA_DIR}/onc089eec2c-SIP0001_S1_L001_R1_001.fastq.gz"
TUMOR_R2="${DATA_DIR}/onc089eec2c-SIP0001_S1_L001_R2_001.fastq.gz"
NORMAL_R1="${DATA_DIR}/onc089eec2c-SIP0002_S2_L001_R1_001.fastq.gz"
NORMAL_R2="${DATA_DIR}/onc089eec2c-SIP0002_S2_L001_R2_001.fastq.gz"

DOCKER_BASE="docker run --rm --gpus \"device=${GPU_DEVICE}\" \
  --user $(id -u):$(id -g) \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v /mnt/storage:/mnt/storage \
  --tmp-dir /tmp \
  ${PB_IMAGE}"

PIPELINE_START=$SECONDS

# ============================================================
# STAGE 1: fq2bam — Tumor (SIP0001)
# BWA-MEM + sort + mark duplicates; BQSR skipped (no known sites)
# ============================================================
echo "[$(date +%H:%M:%S)] Stage 1/3 — Aligning TUMOR (SIP0001)..."
T0=$SECONDS

BQSR_ARGS=""
if [[ -f "${KNOWN_SITES}" ]]; then
  BQSR_ARGS="--knownSites ${KNOWN_SITES} --out-recal-file ${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt"
  echo "[$(date +%H:%M:%S)]   BQSR: enabled"
else
  echo "[$(date +%H:%M:%S)]   BQSR: skipped (${KNOWN_SITES} not found)"
fi

docker run --rm \
  --gpus "device=${GPU_DEVICE}" \
  --user "$(id -u):$(id -g)" \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v /mnt/storage:/mnt/storage \
  "${PB_IMAGE}" \
  pbrun fq2bam \
    --ref "${REF_FASTA}" \
    --in-fq "${TUMOR_R1}" "${TUMOR_R2}" \
            "@RG\tID:${SAMPLE_ID}_T\tLB:lib1\tPL:ILLUMINA\tSM:${SAMPLE_ID}_T\tPU:unit1" \
    --out-bam "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    --tmp-dir /tmp \
    --num-gpus "${NUM_GPUS}" \
    ${BQSR_ARGS}

echo "[$(date +%H:%M:%S)]   Tumor alignment done: $((SECONDS - T0))s"

# ============================================================
# STAGE 2: fq2bam — Normal (SIP0002)
# ============================================================
echo "[$(date +%H:%M:%S)] Stage 2/3 — Aligning NORMAL (SIP0002)..."
T0=$SECONDS

BQSR_ARGS=""
if [[ -f "${KNOWN_SITES}" ]]; then
  BQSR_ARGS="--knownSites ${KNOWN_SITES} --out-recal-file ${OUTDIR}/${SAMPLE_ID}_normal.recal.txt"
fi

docker run --rm \
  --gpus "device=${GPU_DEVICE}" \
  --user "$(id -u):$(id -g)" \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v /mnt/storage:/mnt/storage \
  "${PB_IMAGE}" \
  pbrun fq2bam \
    --ref "${REF_FASTA}" \
    --in-fq "${NORMAL_R1}" "${NORMAL_R2}" \
            "@RG\tID:${SAMPLE_ID}_N\tLB:lib1\tPL:ILLUMINA\tSM:${SAMPLE_ID}_N\tPU:unit1" \
    --out-bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --tmp-dir /tmp \
    --num-gpus "${NUM_GPUS}" \
    ${BQSR_ARGS}

echo "[$(date +%H:%M:%S)]   Normal alignment done: $((SECONDS - T0))s"

# ============================================================
# STAGE 3: mutectcaller — Paired tumor-normal somatic calling
# PoN and gnomAD used only if present in ref dir
# ============================================================
echo "[$(date +%H:%M:%S)] Stage 3/3 — Somatic variant calling (Mutect2 paired)..."
T0=$SECONDS

MUTECT_ARGS=""
[[ -f "${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt" ]]  && MUTECT_ARGS+=" --in-tumor-recal-file ${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt"
[[ -f "${OUTDIR}/${SAMPLE_ID}_normal.recal.txt" ]] && MUTECT_ARGS+=" --in-normal-recal-file ${OUTDIR}/${SAMPLE_ID}_normal.recal.txt"
[[ -f "${PON}" ]]               && MUTECT_ARGS+=" --pon ${PON}"
[[ -f "${GERMLINE_RESOURCE}" ]] && MUTECT_ARGS+=" --germline-resource ${GERMLINE_RESOURCE}"

docker run --rm \
  --gpus "device=${GPU_DEVICE}" \
  --user "$(id -u):$(id -g)" \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v /mnt/storage:/mnt/storage \
  "${PB_IMAGE}" \
  pbrun mutectcaller \
    --ref "${REF_FASTA}" \
    --tumor-name "${SAMPLE_ID}_T" \
    --in-tumor-bam "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    --normal-name "${SAMPLE_ID}_N" \
    --in-normal-bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --out-vcf "${OUTDIR}/${SAMPLE_ID}.somatic.vcf" \
    --tmp-dir /tmp \
    --num-gpus "${NUM_GPUS}" \
    ${MUTECT_ARGS}

echo "[$(date +%H:%M:%S)]   Somatic calling done: $((SECONDS - T0))s"

# ============================================================
# Summary
# ============================================================
VARIANT_COUNT=$(grep -cv "^#" "${OUTDIR}/${SAMPLE_ID}.somatic.vcf" 2>/dev/null || echo "?")

echo ""
echo "=================================================="
echo " Smoke test COMPLETE — total: $((SECONDS - PIPELINE_START))s"
echo "=================================================="
echo "  Tumor BAM:   ${OUTDIR}/${SAMPLE_ID}_tumor.bam"
echo "  Normal BAM:  ${OUTDIR}/${SAMPLE_ID}_normal.bam"
echo "  Somatic VCF: ${OUTDIR}/${SAMPLE_ID}.somatic.vcf  (${VARIANT_COUNT} raw variant records)"
echo ""
echo "  Config: Parabricks ${PB_IMAGE##*:}, GPU device=${GPU_DEVICE}, num_gpus=${NUM_GPUS}"
echo "  BQSR: $([ -f "${KNOWN_SITES}" ] && echo 'enabled' || echo 'skipped — known_indels VCF missing')"
echo "  PoN:  $([ -f "${PON}" ] && echo 'enabled' || echo 'skipped — PoN VCF missing')"
echo "=================================================="
