#!/usr/bin/env bash
#
# PCGR second-opinion report for a somatic VCF (DeepSomatic or Mutect2).
# Runs PCGR v2.2.5 via Docker and produces an HTML + TSV report annotated
# with CIViC, ClinVar, Open Targets Platform, and OncoKB.
#
# PCGR is a deterministic second-opinion complement to the pipeline's
# generative Nemotron/CIViC report: its tiering is rule-based rather than
# model-driven, making the two outputs useful for diff-based QC.
#
# Usage:
#   ./run_pcgr.sh <sample_id> [vcf_path]
#
# <vcf_path> defaults to:
#   ${OUTDIR}/deepsomatic/${SAMPLE}.deepsomatic.pass.vcf.gz
#
# Required env overrides:
#   PCGR_BUNDLE=/mnt/storage/pcgr_bundle   PCGR reference data bundle (v20250314)
#   VEP_CACHE=/mnt/storage/vep_cache       VEP homo_sapiens GRCh38 cache (~30 GB)
#
# Optional env overrides:
#   TUMOR_TYPE=LUAD      OncoTree code (mapped to PCGR site integer; default LUAD→15)
#   OUTDIR=/path/out     base output dir (default /mnt/storage/parabricks_test/output/<sample>)
#   PCGR_VERSION=2.2.5   Docker image tag (default 2.2.5)
#   PCGR_ASSAY=WGS       WGS | WES | TARGETED (default WGS)
#   TUMOR_PURITY=0.6     estimated tumor purity fraction (default 0.6)
#
# One-time data bundle download (~5 GB, GRCh38):
#   mkdir -p "${PCGR_BUNDLE}"
#   curl -L https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz \
#     | tar -xz -C "${PCGR_BUNDLE}"
#
# One-time VEP cache download (~30 GB):
#   docker run --rm \
#     -v "${VEP_CACHE}":/mnt/.vep \
#     sigven/pcgr:2.2.5 \
#     vep_install -a cf -s homo_sapiens -y GRCh38 -c /mnt/.vep --convert
#
# PCGR tumor_site integer codes (v2.2.5):
#   0=Any  4=Bladder  6=Breast  8=CNS/Brain  9=Colon/Rectum  10=Esoph/Stomach
#   12=Head_and_Neck  13=Kidney  14=Liver  15=Lung  16=Lymphoid  17=Myeloid
#   18=Ovary  19=Pancreas  23=Prostate  24=Skin  25=Soft_Tissue  28=Thyroid  29=Uterus

set -euo pipefail

SAMPLE="${1:-HCC1395}"
CUSTOM_VCF="${2:-}"

PCGR_BUNDLE="${PCGR_BUNDLE:-}"
VEP_CACHE="${VEP_CACHE:-}"
TUMOR_TYPE="${TUMOR_TYPE:-LUAD}"
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/${SAMPLE}}"
PCGR_VERSION="${PCGR_VERSION:-2.2.5}"
PCGR_ASSAY="${PCGR_ASSAY:-WGS}"
TUMOR_PURITY="${TUMOR_PURITY:-0.6}"

DS_OUTDIR="${OUTDIR}/deepsomatic"
INPUT_VCF="${CUSTOM_VCF:-${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.vcf.gz}"
PCGR_OUTDIR="${DS_OUTDIR}/pcgr"

# ============================================================
# VALIDATE PREREQUISITES
# ============================================================
echo "========================================================"
echo "[$(date +%H:%M:%S)] PCGR second-opinion report"
echo "========================================================"

PREREQ_OK=true
if [[ -z "${PCGR_BUNDLE}" ]]; then
  echo "ERROR: PCGR_BUNDLE env var not set."
  echo "       See script header for one-time download instructions (~5 GB)."
  PREREQ_OK=false
elif [[ ! -d "${PCGR_BUNDLE}" ]]; then
  echo "ERROR: PCGR_BUNDLE directory not found: ${PCGR_BUNDLE}"
  PREREQ_OK=false
fi
if [[ -z "${VEP_CACHE}" ]]; then
  echo "ERROR: VEP_CACHE env var not set."
  echo "       See script header for one-time VEP cache download instructions (~30 GB)."
  PREREQ_OK=false
elif [[ ! -d "${VEP_CACHE}" ]]; then
  echo "ERROR: VEP_CACHE directory not found: ${VEP_CACHE}"
  PREREQ_OK=false
fi
if [[ ! -f "${INPUT_VCF}" ]]; then
  echo "ERROR: Input VCF not found: ${INPUT_VCF}"
  echo "       Run run_deepsomatic_pipeline.sh first (Phases 1-2 produce the PASS VCF)."
  PREREQ_OK=false
fi
[[ "${PREREQ_OK}" == false ]] && exit 1

# ============================================================
# ONCOTREE → PCGR TUMOR SITE MAPPING
# ============================================================
declare -A _SITE_MAP
_SITE_MAP=(
  # Lung
  [LUAD]=15 [LUSC]=15 [SCLC]=15 [LUCA]=15 [NSCLC]=15
  # Skin / Melanoma
  [SKCM]=24 [MCC]=24 [MEL]=24
  # Breast
  [BRCA]=6 [IDC]=6 [ILC]=6
  # Colon / Rectum
  [COAD]=9 [READ]=9 [COADREAD]=9
  # Pancreas
  [PAAD]=19
  # CNS / Brain
  [GBM]=8 [LGG]=8 [MBL]=8 [BRAIN]=8
  # Myeloid
  [AML]=17 [MDS]=17 [CML]=17 [MPN]=17
  # Lymphoid
  [DLBCL]=16 [FL]=16 [CLL]=16 [HL]=16 [LYMPH]=16
  # Prostate
  [PRAD]=23
  # Ovary
  [OV]=18 [HGSOC]=18
  # Uterus / Endometrium
  [UCEC]=29 [UCS]=29
  # Liver / Biliary
  [HCC]=14 [CHOL]=14 [IHCH]=14
  # Head and Neck
  [HNSC]=12 [SCCHN]=12
  # Kidney
  [KIRC]=13 [KIRP]=13 [CCRCC]=13
  # Esophagus / Stomach
  [ESCA]=10 [STAD]=10
  # Bladder
  [BLCA]=4 [UTUC]=4
  # Thyroid
  [THCA]=28 [THPA]=28
  # Soft Tissue
  [SARC]=25 [LPS]=25 [RMS]=25
)
PCGR_SITE="${_SITE_MAP[${TUMOR_TYPE}]:-0}"
if [[ "${PCGR_SITE}" == "0" ]]; then
  echo "  NOTE: No PCGR site mapping for OncoTree code '${TUMOR_TYPE}' — using site=0 (Any/Other)"
fi

# ============================================================
# RUN PCGR
# ============================================================
# Skip if HTML output already exists
PCGR_HTML="${PCGR_OUTDIR}/${SAMPLE}.pcgr.grch38.html"
if [[ -f "${PCGR_HTML}" ]]; then
  echo "[$(date +%H:%M:%S)] Skipping PCGR — output exists: ${PCGR_HTML}"
else
  mkdir -p "${PCGR_OUTDIR}"

  # PCGR v2.2.5 requires depth/AF as VCF INFO fields (not FORMAT).
  # DeepSomatic stores these in FORMAT/DP and FORMAT/VAF, so we lift
  # them into INFO/TDP and INFO/TVAF using bcftools/bgzip inside the PCGR
  # container (those tools are not assumed to be on the host).
  PREP_VCF="${PCGR_OUTDIR}/${SAMPLE}.pcgr_input.vcf.gz"
  PREP_AWK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/vcf_add_info_dp_vaf.awk"
  INPUT_ABS="$(cd "$(dirname "${INPUT_VCF}")" && pwd)/$(basename "${INPUT_VCF}")"
  INPUT_SRC_DIR="${INPUT_ABS%/*}"
  INPUT_SRC_BASE="$(basename "${INPUT_ABS}")"
  PREP_BASE="$(basename "${PREP_VCF}")"

  if [[ ! -f "${PREP_VCF}" ]]; then
    echo "  Preprocessing VCF: lifting FORMAT/DP → INFO/TDP, FORMAT/VAF → INFO/TVAF..."
    docker container run --rm \
      -v "${INPUT_SRC_DIR}":/mnt/src_vcf:ro \
      -v "${PCGR_OUTDIR}":/mnt/output \
      -v "${PREP_AWK}":/mnt/scripts/vcf_add_info_dp_vaf.awk:ro \
      "sigven/pcgr:${PCGR_VERSION}" \
      bash -c "set -euo pipefail && \
        bcftools view /mnt/src_vcf/${INPUT_SRC_BASE} \
          | awk -f /mnt/scripts/vcf_add_info_dp_vaf.awk \
          | bgzip > /mnt/output/${PREP_BASE} && \
        bcftools index -t /mnt/output/${PREP_BASE}"
    echo "  Preprocessed VCF: ${PREP_VCF}"
  fi

  INPUT_DIR="${PCGR_OUTDIR}"
  INPUT_BASENAME="${PREP_BASE}"

  echo "  Sample:       ${SAMPLE}"
  echo "  Input VCF:    ${PREP_VCF}"
  echo "  Tumor type:   ${TUMOR_TYPE} → PCGR site ${PCGR_SITE}"
  echo "  Assay:        ${PCGR_ASSAY}"
  echo "  Tumor purity: ${TUMOR_PURITY}"
  echo "  PCGR image:   sigven/pcgr:${PCGR_VERSION}"
  echo ""

  docker container run --rm \
    -v "${VEP_CACHE}":/mnt/.vep \
    -v "${PCGR_BUNDLE}":/mnt/bundle \
    -v "${INPUT_DIR}":/mnt/input \
    -v "${PCGR_OUTDIR}":/mnt/output \
    "sigven/pcgr:${PCGR_VERSION}" \
    pcgr \
      --input_vcf "/mnt/input/${INPUT_BASENAME}" \
      --vep_dir "/mnt/.vep" \
      --refdata_dir "/mnt/bundle" \
      --output_dir "/mnt/output" \
      --genome_assembly "grch38" \
      --sample_id "${SAMPLE}" \
      --tumor_site "${PCGR_SITE}" \
      --assay "${PCGR_ASSAY}" \
      --tumor_dp_tag "TDP" \
      --tumor_af_tag "TVAF" \
      --tumor_purity "${TUMOR_PURITY}" \
      --estimate_tmb \
      --estimate_msi \
      --force_overwrite

  echo "[$(date +%H:%M:%S)] PCGR complete."
fi

echo ""
echo "  HTML report: ${PCGR_OUTDIR}/${SAMPLE}.pcgr.grch38.html"
echo "  TSV report:  ${PCGR_OUTDIR}/${SAMPLE}.pcgr.grch38.snvs_indels.tiers.tsv"
echo ""
echo "Diff against pipeline AMP report:"
echo "  Pipeline MAF:  ${DS_OUTDIR}/${SAMPLE}.deepsomatic.maf"
echo "  PCGR TSV:      ${PCGR_OUTDIR}/${SAMPLE}.pcgr.grch38.snvs_indels.tiers.tsv"
echo ""
echo "Quick comparison (bash):"
echo "  join -t\$'\\t' -1 1 -2 1 \\"
echo "    <(awk -F'\\t' 'NR>1 && \$9<=\"2\"{print \$1}' ${PCGR_OUTDIR}/${SAMPLE}.pcgr.grch38.snvs_indels.tiers.tsv | sort) \\"
echo "    <(awk -F'\\t' 'NR>1 && \$15!=\"\"{print \$1}' ${DS_OUTDIR}/${SAMPLE}.deepsomatic.maf | sort)"
