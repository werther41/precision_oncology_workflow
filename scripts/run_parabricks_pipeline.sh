#!/usr/bin/env bash
#
# Parabricks GPU-Accelerated Precision Oncology Pipeline
# Tumor-Normal paired somatic variant calling: FASTQ -> BAM -> VCF
#
# Requirements:
#   - NVIDIA GPU (RTX PRO 6000 / A100 / H100), driver >= 525
#   - Docker + nvidia-container-toolkit
#   - Parabricks 4.7+ container: nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
#   - BQSR (optional): known indels + dbSNP VCFs from GATK Broad bundle
#   - Somatic filtering (optional): gnomAD AF VCF from Broad bundle
#   - ~80 GB GPU memory for WGS; panel data works on a single GPU
#
# Usage:
#   ./run_parabricks_pipeline.sh <sample_id> <tumor_R1> <tumor_R2> <normal_R1> <normal_R2>
#
# Optional env overrides:
#   NUM_GPUS=2          number of GPUs to use (default 1)
#   GPU_DEVICE=0        which GPU device index to expose (default 0)
#   REF_DIR=/data/ref   path to GRCh38 reference bundle
#
# gnomAD and PoN from the Broad bundle target full GRCh38 (3366 contigs); this pipeline uses
# the no-alt reference (455 sequences). Stage 4 subsets both VCFs with bcftools and builds
# the prepon binary index automatically before the Mutect2 call.

set -euo pipefail

# ============================================================
# CONFIG
# ============================================================
SAMPLE_ID="${1:-${SAMPLE_ID:-PATIENT001}}"
TUMOR_R1="${TUMOR_R1:-${2:-sample_data/tumor_R1.fastq.gz}}"
TUMOR_R2="${TUMOR_R2:-${3:-sample_data/tumor_R2.fastq.gz}}"
NORMAL_R1="${NORMAL_R1:-${4:-sample_data/normal_R1.fastq.gz}}"
NORMAL_R2="${NORMAL_R2:-${5:-sample_data/normal_R2.fastq.gz}}"

NUM_GPUS="${NUM_GPUS:-1}"
GPU_DEVICE="${GPU_DEVICE:-0}"      # use "0,1" for 2-GPU runs

# Reference data (GRCh38 no-alt build, 455 sequences)
REF_DIR="${REF_DIR:-/mnt/storage/parabricks_test/ref}"
REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"
KNOWN_SITES="${REF_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz"
DBSNP="${REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
GERMLINE_RESOURCE="${REF_DIR}/af-only-gnomad.hg38.vcf.gz"
PON_VCF="${REF_DIR}/1000g_pon.hg38.vcf.gz"
GNOMAD_NOALT="${REF_DIR}/af-only-gnomad.noalt.vcf.gz"
PON_NOALT_VCF="${REF_DIR}/1000g_pon.noalt.vcf.gz"
PON_BINARY="${REF_DIR}/1000g_pon.noalt.vcf.gz.pon"
CHROM_LIST="${REF_DIR}/chrom_list.txt"
INTERVALS="${REF_DIR}/wgs_calling_regions.hg38.interval_list"

# Output
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/${SAMPLE_ID}}"
mkdir -p "${OUTDIR}"

# Best-practice filtering artifacts (produced / consumed by GATK CPU stages)
F1R2_TAR="${OUTDIR}/${SAMPLE_ID}_f1r2.tar.gz"
ORIENTATION_MODEL="${OUTDIR}/${SAMPLE_ID}_read_orientation_model.tar.gz"
TUMOR_PILEUPS="${OUTDIR}/${SAMPLE_ID}_tumor_pileups.table"
NORMAL_PILEUPS="${OUTDIR}/${SAMPLE_ID}_normal_pileups.table"
CONTAMINATION_TABLE="${OUTDIR}/${SAMPLE_ID}_contamination.table"
TUMOR_SEGMENTS="${OUTDIR}/${SAMPLE_ID}_tumor_segments.table"

# Parabricks container
PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
DOCKER_RUN="docker run --gpus \"device=${GPU_DEVICE}\" --rm \
  --user $(id -u):$(id -g) \
  -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro \
  -v /mnt/storage:/mnt/storage \
  -v /tmp:/tmp \
  -e TMPDIR=/tmp \
  ${PB_IMAGE}"

# BQSR known-sites args (shared by both fq2bam stages)
BQSR_ARGS=""
if [[ -f "${KNOWN_SITES}" ]]; then
  BQSR_ARGS="--knownSites ${KNOWN_SITES}"
  [[ -f "${DBSNP}" ]] && BQSR_ARGS+=" --knownSites ${DBSNP}"
fi

# ============================================================
# STAGE 1: fq2bam — Tumor sample
# Fused: BWA-MEM alignment + sort + mark duplicates [+ BQSR if known sites present]
# ============================================================
if [[ ! -f "${OUTDIR}/${SAMPLE_ID}_tumor.bam" ]]; then
  echo "[$(date +%H:%M:%S)] Aligning TUMOR ${SAMPLE_ID}..."
  if [[ -n "${BQSR_ARGS}" ]]; then
    echo "[$(date +%H:%M:%S)]   BQSR: enabled (known_indels + dbsnp)"
  else
    echo "[$(date +%H:%M:%S)]   BQSR: skipped — run scripts/download_refs.sh to enable"
  fi
  ${DOCKER_RUN} pbrun fq2bam \
    --ref "${REF_FASTA}" \
    --in-fq "${TUMOR_R1}" "${TUMOR_R2}" \
            "@RG\tID:${SAMPLE_ID}_T\tLB:lib1\tPL:ILLUMINA\tSM:${SAMPLE_ID}_T\tPU:unit1" \
    --out-bam "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    --out-recal-file "${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt" \
    --tmp-dir /tmp \
    --num-gpus "${NUM_GPUS}" \
    --bwa-options="-Y -K 100000000" \
    ${BQSR_ARGS}
else
  echo "[$(date +%H:%M:%S)] Skipping TUMOR alignment — BAM exists"
fi

# ============================================================
# STAGE 2: fq2bam — Normal sample
# ============================================================
if [[ ! -f "${OUTDIR}/${SAMPLE_ID}_normal.bam" ]]; then
  echo "[$(date +%H:%M:%S)] Aligning NORMAL ${SAMPLE_ID}..."
  ${DOCKER_RUN} pbrun fq2bam \
    --ref "${REF_FASTA}" \
    --in-fq "${NORMAL_R1}" "${NORMAL_R2}" \
            "@RG\tID:${SAMPLE_ID}_N\tLB:lib1\tPL:ILLUMINA\tSM:${SAMPLE_ID}_N\tPU:unit1" \
    --out-bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --out-recal-file "${OUTDIR}/${SAMPLE_ID}_normal.recal.txt" \
    --tmp-dir /tmp \
    --num-gpus "${NUM_GPUS}" \
    --bwa-options="-Y -K 100000000" \
    ${BQSR_ARGS}
else
  echo "[$(date +%H:%M:%S)] Skipping NORMAL alignment — BAM exists"
fi

# ============================================================
# STAGE 3: QC — Coverage and alignment metrics
# ============================================================
if [[ ! -d "${OUTDIR}/qc_tumor" ]]; then
  echo "[$(date +%H:%M:%S)] Running QC (tumor)..."
  ${DOCKER_RUN} pbrun collectmultiplemetrics \
    --ref "${REF_FASTA}" \
    --bam "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    --out-qc-metrics-dir "${OUTDIR}/qc_tumor" \
    --gen-all-metrics \
    --tmp-dir /tmp
else
  echo "[$(date +%H:%M:%S)] Skipping QC tumor — metrics exist"
fi

if [[ ! -d "${OUTDIR}/qc_normal" ]]; then
  echo "[$(date +%H:%M:%S)] Running QC (normal)..."
  ${DOCKER_RUN} pbrun collectmultiplemetrics \
    --ref "${REF_FASTA}" \
    --bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --out-qc-metrics-dir "${OUTDIR}/qc_normal" \
    --gen-all-metrics \
    --tmp-dir /tmp
else
  echo "[$(date +%H:%M:%S)] Skipping QC normal — metrics exist"
fi

# ============================================================
# STAGE 4: Subset gnomAD + PoN to no-alt contigs; build prepon index
# The Broad bundle VCFs cover all 3366 GRCh38 contigs; Parabricks rejects them
# against the no-alt reference (455 sequences). bcftools restricts both to the
# chromosomes present in this FASTA, then prepon builds the binary PoN index.
# ============================================================
GATK_RUN="docker run --rm \
  -v ${REF_DIR}:${REF_DIR} \
  -v /tmp:/tmp \
  broadinstitute/gatk:4.5.0.0"

if [[ -f "${GERMLINE_RESOURCE}" && ! -f "${GNOMAD_NOALT}" ]]; then
  echo "[$(date +%H:%M:%S)] Subsetting gnomAD to no-alt contigs..."
  [[ ! -f "${CHROM_LIST}" ]] && \
    awk '{print $1"\t0\t"$2}' "${REF_FASTA}.fai" > "${CHROM_LIST}"
  ${GATK_RUN} bcftools view -R "${CHROM_LIST}" "${GERMLINE_RESOURCE}" -Oz -o "${GNOMAD_NOALT}"
  # bcftools view -R filters data rows but leaves 3366 contig lines in the header.
  # Parabricks mutect compares PoN/germline contig order against the FASTA (FAI order).
  # Rewrite header with contig lines in FAI lexicographic order to match.
  ${GATK_RUN} bash -c "
    bcftools view -h '${GNOMAD_NOALT}' | grep -v '^##contig' | grep -v '^#CHROM' > /tmp/_gnomad_hdr.txt
    awk '{printf \"##contig=<ID=%s,length=%s>\n\",\$1,\$2}' '${REF_FASTA}.fai' >> /tmp/_gnomad_hdr.txt
    bcftools view -h '${GNOMAD_NOALT}' | grep '^#CHROM' >> /tmp/_gnomad_hdr.txt
    bcftools reheader -h /tmp/_gnomad_hdr.txt '${GNOMAD_NOALT}' -o '${GNOMAD_NOALT}.tmp'
    mv '${GNOMAD_NOALT}.tmp' '${GNOMAD_NOALT}'"
  ${GATK_RUN} bcftools index -t -f "${GNOMAD_NOALT}"
elif [[ ! -f "${GERMLINE_RESOURCE}" ]]; then
  echo "[$(date +%H:%M:%S)] gnomAD not found — run download_refs.sh to enable germline filtering"
else
  echo "[$(date +%H:%M:%S)] Skipping gnomAD subset — noalt file exists"
fi

if [[ -f "${PON_VCF}" && ! -f "${PON_NOALT_VCF}" ]]; then
  echo "[$(date +%H:%M:%S)] Subsetting PoN to no-alt contigs..."
  [[ ! -f "${CHROM_LIST}" ]] && \
    awk '{print $1"\t0\t"$2}' "${REF_FASTA}.fai" > "${CHROM_LIST}"
  ${GATK_RUN} bcftools view -R "${CHROM_LIST}" "${PON_VCF}" -Oz -o "${PON_NOALT_VCF}"
  ${GATK_RUN} bash -c "
    bcftools view -h '${PON_NOALT_VCF}' | grep -v '^##contig' | grep -v '^#CHROM' > /tmp/_pon_hdr.txt
    awk '{printf \"##contig=<ID=%s,length=%s>\n\",\$1,\$2}' '${REF_FASTA}.fai' >> /tmp/_pon_hdr.txt
    bcftools view -h '${PON_NOALT_VCF}' | grep '^#CHROM' >> /tmp/_pon_hdr.txt
    bcftools reheader -h /tmp/_pon_hdr.txt '${PON_NOALT_VCF}' -o '${PON_NOALT_VCF}.tmp'
    mv '${PON_NOALT_VCF}.tmp' '${PON_NOALT_VCF}'"
  ${GATK_RUN} bcftools index -t -f "${PON_NOALT_VCF}"
elif [[ ! -f "${PON_VCF}" ]]; then
  echo "[$(date +%H:%M:%S)] PoN VCF not found — run download_refs.sh to enable PoN filtering"
else
  echo "[$(date +%H:%M:%S)] Skipping PoN subset — noalt VCF exists"
fi

if [[ -f "${PON_NOALT_VCF}" && ! -f "${PON_BINARY}" ]]; then
  echo "[$(date +%H:%M:%S)] Building PoN binary index (pbrun prepon)..."
  ${DOCKER_RUN} pbrun prepon \
    --in-pon-file "${PON_NOALT_VCF}" \
    --num-gpus "${NUM_GPUS}" \
    --tmp-dir /tmp
  # prepon appends .pon to the full input filename (e.g. foo.vcf.gz -> foo.vcf.gz.pon)
  PREPON_AUTO="${PON_NOALT_VCF}.pon"
  [[ -f "${PREPON_AUTO}" && "${PREPON_AUTO}" != "${PON_BINARY}" ]] && mv "${PREPON_AUTO}" "${PON_BINARY}"
else
  [[ -f "${PON_BINARY}" ]] && echo "[$(date +%H:%M:%S)] Skipping prepon — PoN binary exists"
fi

# ============================================================
# STAGE 5: Somatic variant calling — Mutect2 (paired)
# ============================================================
if [[ ! -f "${OUTDIR}/${SAMPLE_ID}.somatic.unfiltered.vcf.gz" || ! -f "${F1R2_TAR}" ]]; then
  echo "[$(date +%H:%M:%S)] Calling somatic variants with Mutect2..."
  MUTECT_ARGS=""
  [[ -f "${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt" ]]  && MUTECT_ARGS+=" --in-tumor-recal-file ${OUTDIR}/${SAMPLE_ID}_tumor.recal.txt"
  [[ -f "${OUTDIR}/${SAMPLE_ID}_normal.recal.txt" ]] && MUTECT_ARGS+=" --in-normal-recal-file ${OUTDIR}/${SAMPLE_ID}_normal.recal.txt"
  [[ -f "${INTERVALS}" ]]         && MUTECT_ARGS+=" --interval-file ${INTERVALS}"
  [[ -f "${GNOMAD_NOALT}" ]]      && MUTECT_ARGS+=" --mutect-germline-resource ${GNOMAD_NOALT}"
  [[ -f "${PON_BINARY}" ]]        && MUTECT_ARGS+=" --pon ${PON_NOALT_VCF}"

  ${DOCKER_RUN} pbrun mutectcaller \
    --ref "${REF_FASTA}" \
    --tumor-name "${SAMPLE_ID}_T" \
    --in-tumor-bam "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    --normal-name "${SAMPLE_ID}_N" \
    --in-normal-bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --out-vcf "${OUTDIR}/${SAMPLE_ID}.somatic.unfiltered.vcf.gz" \
    --mutect-f1r2-tar-gz "${F1R2_TAR}" \
    --num-gpus "${NUM_GPUS}" \
    --tmp-dir /tmp \
    ${MUTECT_ARGS}
else
  echo "[$(date +%H:%M:%S)] Skipping Mutect2 — unfiltered VCF and f1r2 tar exist"
fi

# ============================================================
# STAGE 6: Germline variant calling — DeepVariant on normal
# Identifies inherited cancer predisposition variants (BRCA1/2, TP53, etc.)
# ============================================================
if [[ ! -f "${OUTDIR}/${SAMPLE_ID}.germline.vcf.gz" ]]; then
  echo "[$(date +%H:%M:%S)] Calling germline variants with DeepVariant..."
  ${DOCKER_RUN} pbrun deepvariant \
    --ref "${REF_FASTA}" \
    --in-bam "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    --out-variants "${OUTDIR}/${SAMPLE_ID}.germline.vcf.gz" \
    --num-gpus "${NUM_GPUS}" \
    --tmp-dir /tmp \
    --mode shortread
else
  echo "[$(date +%H:%M:%S)] Skipping DeepVariant — germline VCF exists"
fi

# ============================================================
# STAGE 7a: Read orientation model (LearnReadOrientationModel — CPU)
# Requires f1r2 tar produced by mutectcaller (--mutect-f1r2-tar-gz).
# ============================================================
if [[ -f "${F1R2_TAR}" && ! -f "${ORIENTATION_MODEL}" ]]; then
  echo "[$(date +%H:%M:%S)] Learning read orientation model..."
  docker run --rm \
    -v "${OUTDIR}:${OUTDIR}" -v /tmp:/tmp -e TMPDIR=/tmp \
    broadinstitute/gatk:4.5.0.0 \
    gatk LearnReadOrientationModel \
      -I "${F1R2_TAR}" \
      -O "${ORIENTATION_MODEL}"
elif [[ -f "${ORIENTATION_MODEL}" ]]; then
  echo "[$(date +%H:%M:%S)] Skipping LearnReadOrientationModel — model exists"
fi

# ============================================================
# STAGE 7b: Pileup summaries + contamination estimate (CPU)
# ============================================================
GATK_RUN_FULL="docker run --rm -v /mnt/storage:/mnt/storage -v /tmp:/tmp -e TMPDIR=/tmp broadinstitute/gatk:4.5.0.0"

if [[ -f "${GNOMAD_NOALT}" && -f "${INTERVALS}" && ! -f "${TUMOR_PILEUPS}" ]]; then
  echo "[$(date +%H:%M:%S)] GetPileupSummaries — tumor..."
  ${GATK_RUN_FULL} gatk GetPileupSummaries \
    -I "${OUTDIR}/${SAMPLE_ID}_tumor.bam" \
    -V "${GNOMAD_NOALT}" \
    -L "${INTERVALS}" \
    -O "${TUMOR_PILEUPS}"
  awk 'NF == 6 || /^#/' "${TUMOR_PILEUPS}" > "${TUMOR_PILEUPS}.clean" && mv "${TUMOR_PILEUPS}.clean" "${TUMOR_PILEUPS}"
else
  [[ -f "${TUMOR_PILEUPS}" ]] && echo "[$(date +%H:%M:%S)] Skipping tumor pileup — table exists"
fi

if [[ -f "${GNOMAD_NOALT}" && -f "${INTERVALS}" && ! -f "${NORMAL_PILEUPS}" ]]; then
  echo "[$(date +%H:%M:%S)] GetPileupSummaries — normal..."
  ${GATK_RUN_FULL} gatk GetPileupSummaries \
    -I "${OUTDIR}/${SAMPLE_ID}_normal.bam" \
    -V "${GNOMAD_NOALT}" \
    -L "${INTERVALS}" \
    -O "${NORMAL_PILEUPS}"
  awk 'NF == 6 || /^#/' "${NORMAL_PILEUPS}" > "${NORMAL_PILEUPS}.clean" && mv "${NORMAL_PILEUPS}.clean" "${NORMAL_PILEUPS}"
else
  [[ -f "${NORMAL_PILEUPS}" ]] && echo "[$(date +%H:%M:%S)] Skipping normal pileup — table exists"
fi

if [[ -f "${TUMOR_PILEUPS}" && -f "${NORMAL_PILEUPS}" && ! -f "${CONTAMINATION_TABLE}" ]]; then
  echo "[$(date +%H:%M:%S)] Calculating contamination..."
  ${GATK_RUN_FULL} gatk CalculateContamination \
    -I "${TUMOR_PILEUPS}" \
    -matched "${NORMAL_PILEUPS}" \
    -O "${CONTAMINATION_TABLE}" \
    --tumor-segmentation "${TUMOR_SEGMENTS}"
else
  [[ -f "${CONTAMINATION_TABLE}" ]] && echo "[$(date +%H:%M:%S)] Skipping contamination — table exists"
fi

# ============================================================
# STAGE 7c: Filter somatic variants (FilterMutectCalls — CPU, via GATK)
# ============================================================
if [[ ! -f "${OUTDIR}/${SAMPLE_ID}.somatic.filtered.vcf.gz" ]]; then
  echo "[$(date +%H:%M:%S)] Filtering somatic calls..."
  FILTER_ARGS=""
  [[ -f "${CONTAMINATION_TABLE}" ]] && FILTER_ARGS+=" --contamination-table ${CONTAMINATION_TABLE} --tumor-segmentation ${TUMOR_SEGMENTS}"
  [[ -f "${ORIENTATION_MODEL}" ]]   && FILTER_ARGS+=" --ob-priors ${ORIENTATION_MODEL}"
  docker run --rm \
    -v /mnt/storage:/mnt/storage \
    -v /tmp:/tmp -e TMPDIR=/tmp \
    broadinstitute/gatk:4.5.0.0 \
    gatk FilterMutectCalls \
      -R "${REF_FASTA}" \
      -V "${OUTDIR}/${SAMPLE_ID}.somatic.unfiltered.vcf.gz" \
      -O "${OUTDIR}/${SAMPLE_ID}.somatic.filtered.vcf.gz" \
      ${FILTER_ARGS}
else
  echo "[$(date +%H:%M:%S)] Skipping FilterMutectCalls — filtered VCF exists"
fi

# ============================================================
# STAGE 7d: Hard post-filter (VAF + depth + POPAF)
# ============================================================
# Tumor is sample index 1 (FORMAT order: normal, tumor in Mutect2 output)
# POPAF >= 2 means gnomAD AF < 1% — removes LOH-elevated common germline variants
# that pass Mutect2 in copy-number-altered regions (dominant FP class in HCC1395).
FINAL_VCF="${OUTDIR}/${SAMPLE_ID}.somatic.popaf2.vcf.gz"
if [[ ! -f "${FINAL_VCF}" && -f "${OUTDIR}/${SAMPLE_ID}.somatic.filtered.vcf.gz" ]]; then
  echo "[$(date +%H:%M:%S)] Applying hard post-filter (VAF>=5%, depth>=5, POPAF>=2)..."
  docker run --rm \
    -v /mnt/storage:/mnt/storage \
    -v /tmp:/tmp -e TMPDIR=/tmp \
    broadinstitute/gatk:4.5.0.0 bash -c "
      bcftools filter \
        -i 'FILTER=\"PASS\" && FORMAT/AF[1:0] >= 0.05 && FORMAT/AD[1:1] >= 5 && INFO/POPAF >= 2' \
        '${OUTDIR}/${SAMPLE_ID}.somatic.filtered.vcf.gz' \
        -Oz -o '${FINAL_VCF}' && \
      bcftools index -t -f '${FINAL_VCF}'"
else
  echo "[$(date +%H:%M:%S)] Skipping hard post-filter — output exists or filtered VCF missing"
fi

echo "[$(date +%H:%M:%S)] Pipeline complete."
echo "  Somatic VCF (filtered):    ${OUTDIR}/${SAMPLE_ID}.somatic.filtered.vcf.gz"
echo "  Somatic VCF (hard-filter): ${FINAL_VCF}"
echo "  Germline VCF: ${OUTDIR}/${SAMPLE_ID}.germline.vcf.gz"
