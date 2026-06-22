#!/usr/bin/env bash
#
# Precision Oncology Pipeline — Unified Orchestrator
#
# Chains four stages end-to-end:
#   FQ2BAM → DeepSomatic → Open CRAVAT → PCGR
#
# Each stage is idempotent: it skips if its sentinel output already exists.
# Stages can be entered mid-pipeline by providing BAMs or a VCF directly.
#
# Usage (single-sample):
#   bash run_pipeline.sh SAMPLE_ID [OPTIONS]
#
# Usage (batch):
#   bash run_pipeline.sh --manifest samples.tsv [OPTIONS]
#
# Examples:
#   # Full pipeline from FASTQs
#   bash run_pipeline.sh PATIENT001 \
#     --tumor-r1 data/tumor_R1.fastq.gz --tumor-r2 data/tumor_R2.fastq.gz \
#     --normal-r1 data/normal_R1.fastq.gz --normal-r2 data/normal_R2.fastq.gz \
#     --tumor-type LUAD
#
#   # Start from existing BAMs (skip FQ2BAM)
#   bash run_pipeline.sh HCC1395 \
#     --tumor-bam /path/to/HCC1395_tumor.bam \
#     --normal-bam /path/to/HCC1395_normal.bam \
#     --tumor-type BRCA
#
#   # Annotate a pre-called VCF (skip FQ2BAM + DeepSomatic)
#   bash run_pipeline.sh PATIENT001 \
#     --vcf /path/to/patient001.pass.vcf.gz \
#     --tumor-type PAAD
#
#   # Batch mode
#   bash run_pipeline.sh --manifest samples.tsv --outdir /results

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================
# LOGGING
# ============================================================
_log_file=""
log()  { printf '[%s] %s\n'   "$(date +%H:%M:%S)" "$*" | tee -a "${_log_file}"; }
step() { printf '\n[%s] ===== %s =====\n' "$(date +%H:%M:%S)" "$*" | tee -a "${_log_file}"; }
err()  { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }
warn() { printf '[%s] WARN:  %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "${_log_file}"; }

# ============================================================
# DEFAULTS (overridden by pipeline.conf, then CLI)
# ============================================================
SAMPLE=""
TUMOR_R1=""; TUMOR_R2=""
NORMAL_R1=""; NORMAL_R2=""
TUMOR_BAM=""; NORMAL_BAM=""
INPUT_VCF=""
CONFIG_FILE="${SCRIPT_DIR}/pipeline.conf"
MANIFEST=""
START_FROM=""
STOP_AFTER="pcgr"
SKIP_OC=false
SKIP_PCGR=false
DRY_RUN=false

# Values populated from config then overridable by CLI
OUTDIR=""
TUMOR_TYPE=""
NUM_GPUS=""
GPU_DEVICE=""
MODEL_TYPE=""
TUMOR_PURITY=""
PCGR_ASSAY=""
WORKERS=""
OC_ANNOTATORS=""
PCGR_BUNDLE=""
VEP_CACHE=""
OC_BIN=""
ONCOKB_TOKEN=""
REF_DIR=""
REF_FASTA=""
KNOWN_SITES=""
GNOMAD_VCF=""
PON_VCF=""
PARABRICKS_IMAGE=""
PCGR_IMAGE=""

# ============================================================
# CLI PARSING
# ============================================================
usage() {
  cat <<EOF
Usage: run_pipeline.sh SAMPLE_ID [OPTIONS]
       run_pipeline.sh --manifest samples.tsv [OPTIONS]

FASTQ inputs (required for full pipeline):
  --tumor-r1 FILE       Tumor R1 FASTQ.gz
  --tumor-r2 FILE       Tumor R2 FASTQ.gz
  --normal-r1 FILE      Normal R1 FASTQ.gz
  --normal-r2 FILE      Normal R2 FASTQ.gz

Intermediate inputs (skip earlier stages):
  --tumor-bam FILE      Use existing tumor BAM — skip FQ2BAM
  --normal-bam FILE     Use existing normal BAM — skip FQ2BAM
  --vcf FILE            Use pre-called PASS VCF — skip FQ2BAM + DeepSomatic

Pipeline control:
  --config FILE         Config file (default: pipeline.conf alongside this script)
  --outdir DIR          Base output directory (overrides config OUTDIR)
  --tumor-type CODE     OncoTree code, e.g. BRCA, LUAD (overrides config TUMOR_TYPE)
  --num-gpus N          Number of GPUs (overrides config NUM_GPUS)
  --gpu-device D        GPU device index, e.g. 0 (overrides config GPU_DEVICE)
  --start-from STAGE    fq2bam | deepsomatic | opencravat | pcgr
  --stop-after STAGE    fq2bam | deepsomatic | opencravat | pcgr  (default: pcgr)
  --skip-oc             Skip Open CRAVAT annotation stage
  --skip-pcgr           Skip PCGR report stage
  --manifest FILE       Batch mode: TSV with columns sample_id tumor_r1 tumor_r2
                        normal_r1 normal_r2 tumor_type tumor_bam normal_bam vcf
  --dry-run             Print what would run without executing
  -h, --help            Show this help

EOF
  exit 0
}

_cli_outdir=""; _cli_tumor_type=""; _cli_num_gpus=""; _cli_gpu_device=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tumor-r1)    TUMOR_R1="$2";    shift 2 ;;
    --tumor-r2)    TUMOR_R2="$2";    shift 2 ;;
    --normal-r1)   NORMAL_R1="$2";   shift 2 ;;
    --normal-r2)   NORMAL_R2="$2";   shift 2 ;;
    --tumor-bam)   TUMOR_BAM="$2";   shift 2 ;;
    --normal-bam)  NORMAL_BAM="$2";  shift 2 ;;
    --vcf)         INPUT_VCF="$2";   shift 2 ;;
    --config)      CONFIG_FILE="$2"; shift 2 ;;
    --outdir)      _cli_outdir="$2"; shift 2 ;;
    --tumor-type)  _cli_tumor_type="$2"; shift 2 ;;
    --num-gpus)    _cli_num_gpus="$2";   shift 2 ;;
    --gpu-device)  _cli_gpu_device="$2"; shift 2 ;;
    --start-from)  START_FROM="$2";  shift 2 ;;
    --stop-after)  STOP_AFTER="$2";  shift 2 ;;
    --skip-oc)     SKIP_OC=true;     shift ;;
    --skip-pcgr)   SKIP_PCGR=true;   shift ;;
    --manifest)    MANIFEST="$2";    shift 2 ;;
    --dry-run)     DRY_RUN=true;     shift ;;
    -h|--help)     usage ;;
    -*)            err "Unknown option: $1. Use -h for help." ;;
    *)
      [[ -z "${SAMPLE}" ]] || err "Unexpected positional argument: $1"
      SAMPLE="$1"; shift ;;
  esac
done

# ============================================================
# LOAD CONFIG
# ============================================================
if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck source=/dev/null
  source "${CONFIG_FILE}"
else
  warn "Config file not found: ${CONFIG_FILE} — using built-in defaults"
fi

# Apply CLI overrides (take precedence over config)
[[ -n "${_cli_outdir}"      ]] && OUTDIR="${_cli_outdir}"
[[ -n "${_cli_tumor_type}"  ]] && TUMOR_TYPE="${_cli_tumor_type}"
[[ -n "${_cli_num_gpus}"    ]] && NUM_GPUS="${_cli_num_gpus}"
[[ -n "${_cli_gpu_device}"  ]] && GPU_DEVICE="${_cli_gpu_device}"

# Apply final defaults for anything still unset
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output}"
TUMOR_TYPE="${TUMOR_TYPE:-LUAD}"
NUM_GPUS="${NUM_GPUS:-1}"
GPU_DEVICE="${GPU_DEVICE:-0}"
MODEL_TYPE="${MODEL_TYPE:-WGS}"
TUMOR_PURITY="${TUMOR_PURITY:-0.6}"
PCGR_ASSAY="${PCGR_ASSAY:-WGS}"
WORKERS="${WORKERS:-5}"
OC_ANNOTATORS="${OC_ANNOTATORS:-civic oncokb}"
PCGR_BUNDLE="${PCGR_BUNDLE:-}"
VEP_CACHE="${VEP_CACHE:-}"
OC_BIN="${OC_BIN:-/mnt/storage/open-cravat-reports/oc-venv/bin/oc}"
ONCOKB_TOKEN="${ONCOKB_TOKEN:-}"
REF_DIR="${REF_DIR:-/mnt/storage/parabricks_test/ref}"

# ============================================================
# BATCH MODE
# ============================================================
if [[ -n "${MANIFEST}" ]]; then
  [[ -f "${MANIFEST}" ]] || err "Manifest not found: ${MANIFEST}"
  log "Batch mode: reading manifest ${MANIFEST}"
  tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r _sample _t_r1 _t_r2 _n_r1 _n_r2 _ttype _t_bam _n_bam _vcf; do
    [[ -z "${_sample}" || "${_sample}" == "#"* ]] && continue
    log "--- Batch sample: ${_sample} ---"
    CMD=(bash "${BASH_SOURCE[0]}" "${_sample}")
    [[ "${_t_r1}"  != "-" && -n "${_t_r1}"  ]] && CMD+=(--tumor-r1  "${_t_r1}")
    [[ "${_t_r2}"  != "-" && -n "${_t_r2}"  ]] && CMD+=(--tumor-r2  "${_t_r2}")
    [[ "${_n_r1}"  != "-" && -n "${_n_r1}"  ]] && CMD+=(--normal-r1 "${_n_r1}")
    [[ "${_n_r2}"  != "-" && -n "${_n_r2}"  ]] && CMD+=(--normal-r2 "${_n_r2}")
    [[ "${_t_bam}" != "-" && -n "${_t_bam}" ]] && CMD+=(--tumor-bam  "${_t_bam}")
    [[ "${_n_bam}" != "-" && -n "${_n_bam}" ]] && CMD+=(--normal-bam "${_n_bam}")
    [[ "${_vcf}"   != "-" && -n "${_vcf}"   ]] && CMD+=(--vcf        "${_vcf}")
    [[ "${_ttype}" != "-" && -n "${_ttype}" ]] && CMD+=(--tumor-type "${_ttype}")
    CMD+=(--config "${CONFIG_FILE}" --outdir "${OUTDIR}")
    [[ "${DRY_RUN}" == true ]] && CMD+=(--dry-run)
    [[ "${SKIP_OC}" == true ]] && CMD+=(--skip-oc)
    [[ "${SKIP_PCGR}" == true ]] && CMD+=(--skip-pcgr)
    "${CMD[@]}"
  done
  exit 0
fi

# ============================================================
# SINGLE-SAMPLE VALIDATION
# ============================================================
[[ -z "${SAMPLE}" ]] && err "SAMPLE_ID is required. Use --manifest for batch mode."

SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE}"
_log_file="${SAMPLE_OUTDIR}/pipeline_run.log"
mkdir -p "${SAMPLE_OUTDIR}"

step "Precision Oncology Pipeline — ${SAMPLE}"
log "Config:       ${CONFIG_FILE}"
log "Output dir:   ${SAMPLE_OUTDIR}"
log "Tumor type:   ${TUMOR_TYPE}"
log "GPUs:         ${NUM_GPUS} (device ${GPU_DEVICE})"
log "Dry run:      ${DRY_RUN}"

# Derived paths used across stages
TUMOR_BAM="${TUMOR_BAM:-${SAMPLE_OUTDIR}/${SAMPLE}_tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:-${SAMPLE_OUTDIR}/${SAMPLE}_normal.bam}"
DS_OUTDIR="${SAMPLE_OUTDIR}/deepsomatic"
DS_PASS_VCF="${DS_OUTDIR}/${SAMPLE}.deepsomatic.pass.vcf.gz"
OC_OUTDIR="${SAMPLE_OUTDIR}/opencravat"
OC_REPORT_NAME="${SAMPLE}-oc-report"
OC_TSV="${OC_OUTDIR}/${OC_REPORT_NAME}.tsv"
PCGR_HTML="${DS_OUTDIR}/pcgr/${SAMPLE}.pcgr.grch38.html"

# If a custom VCF was provided, it is the DeepSomatic substitute
[[ -n "${INPUT_VCF}" ]] && DS_PASS_VCF="${INPUT_VCF}"

# ============================================================
# AUTO-DETECT START STAGE
# ============================================================
if [[ -z "${START_FROM}" ]]; then
  if   [[ -n "${INPUT_VCF}" ]];                            then START_FROM="opencravat"
  elif [[ -n "${TUMOR_R1}" ]];                             then START_FROM="fq2bam"
  elif [[ -n "${TUMOR_BAM_ARG:-}" || -f "${TUMOR_BAM}" ]]; then START_FROM="deepsomatic"
  else
    err "Cannot determine start stage. Provide --tumor-r1/r2 (FASTQs), --tumor-bam (BAMs), or --vcf."
  fi
fi

# Normalize stage ordering
declare -A _STAGE_ORDER=([fq2bam]=1 [deepsomatic]=2 [opencravat]=3 [pcgr]=4)
[[ -z "${_STAGE_ORDER[${START_FROM}]+x}" ]] && err "Invalid --start-from: ${START_FROM}"
[[ -z "${_STAGE_ORDER[${STOP_AFTER}]+x}" ]] && err "Invalid --stop-after: ${STOP_AFTER}"

_should_run() {
  local stage="$1"
  local ord="${_STAGE_ORDER[${stage}]}"
  local start="${_STAGE_ORDER[${START_FROM}]}"
  local stop="${_STAGE_ORDER[${STOP_AFTER}]}"
  [[ "${ord}" -ge "${start}" && "${ord}" -le "${stop}" ]]
}

log "Start stage:  ${START_FROM}"
log "Stop after:   ${STOP_AFTER}"
log "Skip OC:      ${SKIP_OC}"
log "Skip PCGR:    ${SKIP_PCGR}"

# ============================================================
# DRY-RUN HELPER
# ============================================================
_run() {
  if [[ "${DRY_RUN}" == true ]]; then
    log "[DRY-RUN] Would run: $*"
  else
    "$@"
  fi
}

# ============================================================
# STAGE 1: FQ2BAM
# ============================================================
run_fq2bam() {
  step "Stage 1/4: FQ2BAM (Parabricks)"

  if [[ -f "${TUMOR_BAM}" && -f "${NORMAL_BAM}" ]]; then
    log "Skipping FQ2BAM — BAMs exist:"
    log "  ${TUMOR_BAM}"
    log "  ${NORMAL_BAM}"
    return 0
  fi

  [[ -z "${TUMOR_R1}"  ]] && err "FQ2BAM requires --tumor-r1"
  [[ -z "${TUMOR_R2}"  ]] && err "FQ2BAM requires --tumor-r2"
  [[ -z "${NORMAL_R1}" ]] && err "FQ2BAM requires --normal-r1"
  [[ -z "${NORMAL_R2}" ]] && err "FQ2BAM requires --normal-r2"
  [[ ! -f "${TUMOR_R1}" ]]  && err "Tumor R1 not found: ${TUMOR_R1}"
  [[ ! -f "${TUMOR_R2}" ]]  && err "Tumor R2 not found: ${TUMOR_R2}"
  [[ ! -f "${NORMAL_R1}" ]] && err "Normal R1 not found: ${NORMAL_R1}"
  [[ ! -f "${NORMAL_R2}" ]] && err "Normal R2 not found: ${NORMAL_R2}"

  log "Tumor:  ${TUMOR_R1} ${TUMOR_R2}"
  log "Normal: ${NORMAL_R1} ${NORMAL_R2}"

  _run env \
    SAMPLE_ID="${SAMPLE}" \
    TUMOR_R1="${TUMOR_R1}" \
    TUMOR_R2="${TUMOR_R2}" \
    NORMAL_R1="${NORMAL_R1}" \
    NORMAL_R2="${NORMAL_R2}" \
    NUM_GPUS="${NUM_GPUS}" \
    GPU_DEVICE="${GPU_DEVICE}" \
    REF_DIR="${REF_DIR}" \
    OUTDIR="${SAMPLE_OUTDIR}" \
    bash "${SCRIPT_DIR}/scripts/run_parabricks_pipeline.sh" "${SAMPLE}"
}

# ============================================================
# STAGE 2: DeepSomatic
# ============================================================
run_deepsomatic() {
  step "Stage 2/4: DeepSomatic"

  if [[ -f "${DS_PASS_VCF}" ]]; then
    log "Skipping DeepSomatic — VCF exists: ${DS_PASS_VCF}"
    return 0
  fi

  [[ ! -f "${TUMOR_BAM}" ]]  && err "Tumor BAM not found: ${TUMOR_BAM}"
  [[ ! -f "${NORMAL_BAM}" ]] && err "Normal BAM not found: ${NORMAL_BAM}"

  log "Tumor BAM:  ${TUMOR_BAM}"
  log "Normal BAM: ${NORMAL_BAM}"

  _run env \
    SAMPLE="${SAMPLE}" \
    TUMOR_BAM="${TUMOR_BAM}" \
    NORMAL_BAM="${NORMAL_BAM}" \
    NUM_GPUS="${NUM_GPUS}" \
    GPU_DEVICE="${GPU_DEVICE}" \
    MODEL_TYPE="${MODEL_TYPE}" \
    TUMOR_TYPE="${TUMOR_TYPE}" \
    REF_DIR="${REF_DIR}" \
    OUTDIR="${SAMPLE_OUTDIR}" \
    ONCOKB_TOKEN="${ONCOKB_TOKEN}" \
    WORKERS="${WORKERS}" \
    bash "${SCRIPT_DIR}/scripts/run_deepsomatic_pipeline.sh" "${SAMPLE}"
}

# ============================================================
# STAGE 3: Open CRAVAT
# ============================================================
run_opencravat() {
  step "Stage 3/4: Open CRAVAT"

  if [[ -f "${OC_TSV}" ]]; then
    log "Skipping Open CRAVAT — TSV exists: ${OC_TSV}"
    return 0
  fi

  if [[ ! -x "${OC_BIN}" ]]; then
    warn "Open CRAVAT binary not found or not executable: ${OC_BIN}"
    warn "Skipping Open CRAVAT stage."
    return 0
  fi

  [[ ! -f "${DS_PASS_VCF}" ]] && err "PASS VCF not found: ${DS_PASS_VCF}"

  mkdir -p "${OC_OUTDIR}"
  local plain_vcf="${OC_OUTDIR}/${SAMPLE}.pass.vcf"

  log "Input VCF:  ${DS_PASS_VCF}"
  log "Report dir: ${OC_OUTDIR}"
  log "Annotators: ${OC_ANNOTATORS}"

  if [[ "${DRY_RUN}" == true ]]; then
    log "[DRY-RUN] Would run: gunzip -c ${DS_PASS_VCF} | oc run ... -a ${OC_ANNOTATORS}"
    return 0
  fi

  gunzip -c "${DS_PASS_VCF}" > "${plain_vcf}"
  # shellcheck disable=SC2086
  "${OC_BIN}" run "${plain_vcf}" \
    -l hg38 \
    -d "${OC_OUTDIR}" \
    -n "${OC_REPORT_NAME}" \
    -t text \
    -a ${OC_ANNOTATORS}
  rm -f "${plain_vcf}"

  log "Open CRAVAT complete: ${OC_TSV}"
}

# ============================================================
# STAGE 4: PCGR
# ============================================================
run_pcgr() {
  step "Stage 4/4: PCGR"

  if [[ -f "${PCGR_HTML}" ]]; then
    log "Skipping PCGR — HTML exists: ${PCGR_HTML}"
    return 0
  fi

  if [[ -z "${PCGR_BUNDLE}" || ! -d "${PCGR_BUNDLE}" ]]; then
    warn "PCGR_BUNDLE not set or not found: ${PCGR_BUNDLE:-<unset>}"
    warn "Skipping PCGR stage. Set PCGR_BUNDLE in pipeline.conf to enable."
    warn "One-time bundle download (~5 GB):"
    warn "  mkdir -p /mnt/storage/pcgr_bundle"
    warn "  curl -L https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz \\"
    warn "    | tar -xz -C /mnt/storage/pcgr_bundle"
    return 0
  fi

  if [[ -z "${VEP_CACHE}" || ! -d "${VEP_CACHE}" ]]; then
    warn "VEP_CACHE not set or not found: ${VEP_CACHE:-<unset>}"
    warn "Skipping PCGR stage. Set VEP_CACHE in pipeline.conf to enable."
    return 0
  fi

  [[ ! -f "${DS_PASS_VCF}" ]] && err "PASS VCF not found: ${DS_PASS_VCF}"

  log "Input VCF:   ${DS_PASS_VCF}"
  log "Tumor type:  ${TUMOR_TYPE}"
  log "PCGR bundle: ${PCGR_BUNDLE}"

  _run env \
    PCGR_BUNDLE="${PCGR_BUNDLE}" \
    VEP_CACHE="${VEP_CACHE}" \
    TUMOR_TYPE="${TUMOR_TYPE}" \
    OUTDIR="${SAMPLE_OUTDIR}" \
    PCGR_ASSAY="${PCGR_ASSAY}" \
    TUMOR_PURITY="${TUMOR_PURITY}" \
    bash "${SCRIPT_DIR}/scripts/run_pcgr.sh" "${SAMPLE}" "${DS_PASS_VCF}"
}

# ============================================================
# MAIN — run stages in order
# ============================================================
_should_run fq2bam    && run_fq2bam
_should_run deepsomatic && run_deepsomatic
if [[ "${SKIP_OC}" == false ]]; then
  _should_run opencravat && run_opencravat
fi
if [[ "${SKIP_PCGR}" == false ]]; then
  _should_run pcgr && run_pcgr
fi

# ============================================================
# SUMMARY
# ============================================================
step "Pipeline Complete — ${SAMPLE}"
log "Output directory: ${SAMPLE_OUTDIR}"
echo ""

_report_sentinel() {
  local label="$1" path="$2"
  if [[ "${DRY_RUN}" == true ]]; then
    printf "  %-20s %s\n" "${label}:" "(dry-run)"
  elif [[ -f "${path}" ]]; then
    printf "  %-20s %s\n" "${label}:" "${path}"
  else
    printf "  %-20s %s\n" "${label}:" "(not produced)"
  fi
}

_report_sentinel "DeepSomatic VCF"  "${DS_PASS_VCF}"
_report_sentinel "Clinical report"  "${DS_OUTDIR}/${SAMPLE}.deepsomatic_clinical_report.txt"
_report_sentinel "Open CRAVAT TSV"  "${OC_TSV}"
_report_sentinel "PCGR HTML"        "${PCGR_HTML}"
_report_sentinel "PCGR TMB"         "${DS_OUTDIR}/pcgr/${SAMPLE}.pcgr.grch38.tmb.tsv"
echo ""
