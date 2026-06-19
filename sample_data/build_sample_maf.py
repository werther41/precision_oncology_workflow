#!/usr/bin/env python3
"""
Build a realistic synthetic OncoKB-annotated MAF for a lung adenocarcinoma case.
Used to demonstrate the clinical report generator without needing a real GPU run.

Real-world coordinates and protein changes are used (GRCh38).
"""
import pandas as pd
from pathlib import Path

variants = [
    # ---- Driver: EGFR L858R - Tier IA (osimertinib FDA-approved for NSCLC) ----
    {
        "Hugo_Symbol": "EGFR", "Chromosome": "chr7", "Start_Position": 55191822,
        "Reference_Allele": "T", "Tumor_Seq_Allele2": "G",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.L858R",
        "t_ref_count": 45, "t_alt_count": 78,    # VAF ~63%
        "n_ref_count": 120, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "pathogenic",
        "HIGHEST_LEVEL": "LEVEL_1",
        "HIGHEST_LEVEL_SUMMARY": "Osimertinib, Erlotinib, Gefitinib, Afatinib (FDA-approved)",
    },
    # ---- Co-driver: TP53 R248Q - Tier IIC ----
    {
        "Hugo_Symbol": "TP53", "Chromosome": "chr17", "Start_Position": 7674220,
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "T",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.R248Q",
        "t_ref_count": 30, "t_alt_count": 55,    # VAF ~65%
        "n_ref_count": 110, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "pathogenic",
        "HIGHEST_LEVEL": "LEVEL_3B",
        "HIGHEST_LEVEL_SUMMARY": "APR-246 (clinical trials)",
    },
    # ---- KRAS G12C - Tier IA (sotorasib FDA-approved 2021) ----
    # Mutually exclusive with EGFR in real biology; included here to show
    # what the report would look like if both were called - good test case.
    {
        "Hugo_Symbol": "KRAS", "Chromosome": "chr12", "Start_Position": 25245350,
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "A",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.G12C",
        "t_ref_count": 60, "t_alt_count": 8,     # VAF ~12% — likely subclonal
        "n_ref_count": 130, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "pathogenic",
        "HIGHEST_LEVEL": "LEVEL_1",
        "HIGHEST_LEVEL_SUMMARY": "Sotorasib, Adagrasib (FDA-approved for KRAS G12C NSCLC)",
    },
    # ---- STK11 loss-of-function - Tier IIC (resistance marker) ----
    {
        "Hugo_Symbol": "STK11", "Chromosome": "chr19", "Start_Position": 1207021,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "A",
        "Variant_Classification": "Nonsense_Mutation",
        "Consequence": "stop_gained",
        "HGVSp_Short": "p.W239*",
        "t_ref_count": 40, "t_alt_count": 42,    # VAF ~51%
        "n_ref_count": 115, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "pathogenic",
        "HIGHEST_LEVEL": "LEVEL_4",
        "HIGHEST_LEVEL_SUMMARY": "Associated with reduced response to immunotherapy",
    },
    # ---- VUS: PIK3CA novel missense ----
    {
        "Hugo_Symbol": "PIK3CA", "Chromosome": "chr3", "Start_Position": 179234297,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "C",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.R412P",                 # not a hotspot - VUS
        "t_ref_count": 55, "t_alt_count": 22,
        "n_ref_count": 125, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "uncertain_significance",
        "HIGHEST_LEVEL": "",
        "HIGHEST_LEVEL_SUMMARY": "",
    },
    # ---- VUS: KEAP1 missense ----
    {
        "Hugo_Symbol": "KEAP1", "Chromosome": "chr19", "Start_Position": 10491676,
        "Reference_Allele": "G", "Tumor_Seq_Allele2": "T",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.G364C",
        "t_ref_count": 50, "t_alt_count": 18,
        "n_ref_count": 118, "n_alt_count": 0,
        "gnomAD_AF": 0.0, "CLIN_SIG": "uncertain_significance",
        "HIGHEST_LEVEL": "",
        "HIGHEST_LEVEL_SUMMARY": "",
    },
    # ---- Polymorphism that should be filtered out (gnomAD > 1%) ----
    {
        "Hugo_Symbol": "MUC16", "Chromosome": "chr19", "Start_Position": 8959520,
        "Reference_Allele": "C", "Tumor_Seq_Allele2": "T",
        "Variant_Classification": "Missense_Mutation",
        "Consequence": "missense_variant",
        "HGVSp_Short": "p.S5824L",
        "t_ref_count": 38, "t_alt_count": 35,
        "n_ref_count": 60, "n_alt_count": 55,    # also in normal -> germline
        "gnomAD_AF": 0.18, "CLIN_SIG": "benign",  # common variant
        "HIGHEST_LEVEL": "",
        "HIGHEST_LEVEL_SUMMARY": "",
    },
]

df = pd.DataFrame(variants)
out = Path("/home/claude/precision_oncology_workflow/sample_data/PATIENT001.oncokb.maf")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, sep="\t", index=False)
print(f"Wrote {len(df)} variants to {out}")

# Patient metadata
import json
meta = {
    "patient_id": "PATIENT001",
    "specimen": "FFPE tumor (lung lobectomy) + matched peripheral blood",
    "tumor_cov": 285,
    "normal_cov": 64,
    "purity": "0.72 (estimated by FACETS)",
    "tmb": "8.4",
    "msi": "MSS (stable)",
}
meta_path = out.parent / "PATIENT001_meta.json"
meta_path.write_text(json.dumps(meta, indent=2))
print(f"Wrote patient metadata to {meta_path}")
