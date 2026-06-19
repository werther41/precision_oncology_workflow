#!/usr/bin/env bash
#
# Download all missing reference files for the full production pipeline.
# All files served via public HTTPS — no GCS auth or AWS credentials needed.
#
# Estimated total download: ~4.7 GB
# Destination: /mnt/storage/parabricks_test/ref/
#
# Run once; idempotent (skips files that already exist).

set -euo pipefail

REF_DIR="${REF_DIR:-/mnt/storage/parabricks_test/ref}"
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

BROAD="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
GATK="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"

download() {
  local url="$1" file
  file="$(basename "$url")"
  if [[ -f "${file}" ]]; then
    echo "[skip] ${file} already exists"
    return
  fi
  echo "[download] ${file} ..."
  wget -q --show-progress -O "${file}.tmp" "${url}"
  mv "${file}.tmp" "${file}"
  echo "[done] ${file}"
}

echo "=== Downloading BQSR known sites (Broad bundle) ==="
download "${BROAD}/Homo_sapiens_assembly38.known_indels.vcf.gz"
download "${BROAD}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
download "${BROAD}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
download "${BROAD}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
download "${BROAD}/wgs_calling_regions.hg38.interval_list"

echo ""
echo "=== Downloading Mutect2 somatic resources (GATK best practices) ==="
download "${GATK}/af-only-gnomad.hg38.vcf.gz"
download "${GATK}/af-only-gnomad.hg38.vcf.gz.tbi"
download "${GATK}/1000g_pon.hg38.vcf.gz"
download "${GATK}/1000g_pon.hg38.vcf.gz.tbi"

echo ""
echo "=== Done. Ref dir contents ==="
ls -lh "${REF_DIR}"
