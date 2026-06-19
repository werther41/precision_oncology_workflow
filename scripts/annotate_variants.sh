#!/usr/bin/env bash
#
# Tertiary analysis: Annotate filtered somatic VCF with clinical context
# Tools: bcftools, Ensembl VEP, OncoKB annotator
#
# Output: an annotated VCF + a MAF file ready for clinical review

set -euo pipefail

SAMPLE_ID="${1:-PATIENT001}"
VCF_IN="${2:-/data/output/${SAMPLE_ID}/${SAMPLE_ID}.somatic.filtered.vcf.gz}"
OUTDIR="/data/output/${SAMPLE_ID}/annotated"
mkdir -p "${OUTDIR}"

VEP_CACHE="/data/reference/vep_cache"          # downloaded once from Ensembl
ONCOKB_TOKEN="${ONCOKB_TOKEN:?Set ONCOKB_TOKEN env var with your API key}"

# ============================================================
# 1. Keep only PASS variants
# ============================================================
bcftools view -f PASS "${VCF_IN}" -Oz -o "${OUTDIR}/${SAMPLE_ID}.pass.vcf.gz"
bcftools index -t "${OUTDIR}/${SAMPLE_ID}.pass.vcf.gz"

# ============================================================
# 2. Normalize: left-align indels, split multi-allelic
# ============================================================
bcftools norm -m -any -f /data/reference/Homo_sapiens_assembly38.fasta \
  "${OUTDIR}/${SAMPLE_ID}.pass.vcf.gz" \
  -Oz -o "${OUTDIR}/${SAMPLE_ID}.norm.vcf.gz"

# ============================================================
# 3. Run VEP — functional consequence + gnomAD + ClinVar + COSMIC
# ============================================================
docker run --rm -v "${OUTDIR}:${OUTDIR}" -v "${VEP_CACHE}:/opt/vep/.vep" \
  ensemblorg/ensembl-vep:release_111.0 \
  vep \
    --input_file "${OUTDIR}/${SAMPLE_ID}.norm.vcf.gz" \
    --output_file "${OUTDIR}/${SAMPLE_ID}.vep.vcf" \
    --vcf --everything --cache --offline \
    --assembly GRCh38 \
    --fasta /opt/vep/.vep/Homo_sapiens_assembly38.fasta \
    --plugin CADD,/opt/vep/.vep/whole_genome_SNVs.tsv.gz \
    --plugin REVEL,/opt/vep/.vep/revel_grch38.tsv.gz \
    --custom /opt/vep/.vep/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN \
    --custom /opt/vep/.vep/gnomad.exomes.r2.1.1.sites.vcf.gz,gnomADe,vcf,exact,0,AF \
    --custom /opt/vep/.vep/CosmicCodingMuts.vcf.gz,COSMIC,vcf,exact,0,CNT

# ============================================================
# 4. Convert VCF -> MAF (Mutation Annotation Format)
# vcf2maf is the standard intermediate for clinical pipelines
# ============================================================
docker run --rm -v "${OUTDIR}:${OUTDIR}" \
  vanallenlab/vcf2maf:v1.6.21 \
  vcf2maf.pl \
    --input-vcf "${OUTDIR}/${SAMPLE_ID}.vep.vcf" \
    --output-maf "${OUTDIR}/${SAMPLE_ID}.maf" \
    --tumor-id "${SAMPLE_ID}_T" \
    --normal-id "${SAMPLE_ID}_N" \
    --ncbi-build GRCh38 \
    --inhibit-vep

# ============================================================
# 5. OncoKB annotation — therapy-level evidence (Tier I–IV)
# ============================================================
docker run --rm -v "${OUTDIR}:${OUTDIR}" \
  -e ONCOKB_TOKEN="${ONCOKB_TOKEN}" \
  oncokb/oncokb-annotator:v3.4.1 \
  python MafAnnotator.py \
    -i "${OUTDIR}/${SAMPLE_ID}.maf" \
    -o "${OUTDIR}/${SAMPLE_ID}.oncokb.maf" \
    -b "${ONCOKB_TOKEN}" \
    -t "LUAD"   # tumor type code (OncoTree); change per patient

echo "Annotation complete. Final MAF: ${OUTDIR}/${SAMPLE_ID}.oncokb.maf"
