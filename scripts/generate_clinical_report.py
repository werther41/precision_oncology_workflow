#!/usr/bin/env python3
"""
Clinical Report Generator for Precision Oncology
Reads an OncoKB-annotated MAF and produces an AMP/ASCO/CAP tiered clinical report.

Tiers (AMP 2017 guidelines):
  Tier IA  - FDA-approved therapy for this tumor type
  Tier IB  - Well-powered evidence in this tumor type
  Tier IIC - FDA-approved therapy for different tumor type
  Tier IID - Clinical trial evidence
  Tier III - Variant of unknown significance
  Tier IV  - Benign / likely benign
"""

import pandas as pd
import json
import argparse
from datetime import datetime
from pathlib import Path

# OncoKB level -> AMP tier mapping
ONCOKB_TO_TIER = {
    "LEVEL_1":  "Tier IA",   # FDA-recognized biomarker, this tumor type
    "LEVEL_2":  "Tier IB",   # Standard-of-care biomarker, this tumor type
    "LEVEL_3A": "Tier IIC",  # Clinical evidence, this tumor type
    "LEVEL_3B": "Tier IIC",  # Clinical evidence, other tumor type
    "LEVEL_4":  "Tier IID",  # Compelling biological evidence
    "LEVEL_R1": "Tier IA",   # Resistance, standard of care
    "LEVEL_R2": "Tier IID",  # Resistance, clinical evidence
}

CONSEQUENCE_FILTER = {
    "missense_variant", "stop_gained", "stop_lost", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "start_lost",
    "inframe_insertion", "inframe_deletion", "protein_altering_variant",
}


def load_maf(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    # standardize column names we'll touch
    return df


def tier_variant(row) -> str:
    """Assign AMP tier from OncoKB highest level + functional impact."""
    level = str(row.get("HIGHEST_LEVEL", "")).strip()
    if level in ONCOKB_TO_TIER:
        return ONCOKB_TO_TIER[level]
    # No actionable OncoKB level
    consequence = str(row.get("Consequence", ""))
    if any(c in consequence for c in CONSEQUENCE_FILTER):
        clnsig = str(row.get("CLIN_SIG", "")).lower()
        if "pathogenic" in clnsig:
            return "Tier IIC"
        return "Tier III"
    return "Tier IV"


def filter_clinically_relevant(df: pd.DataFrame, min_alt_reads: int = 5) -> pd.DataFrame:
    """Keep variants worth reviewing: protein-altering, VAF >= 5%, not common pop variant."""
    df = df.copy()

    # Variant allele frequency
    df["VAF"] = df["t_alt_count"] / (df["t_alt_count"] + df["t_ref_count"]).replace(0, 1)

    # Filters
    keep = (
        df["Consequence"].isin(CONSEQUENCE_FILTER) &
        (df["VAF"] >= 0.05) &
        (df["t_alt_count"] >= min_alt_reads) &
        (df["gnomAD_AF"].fillna(0).astype(float) < 0.01)  # rare in population
    )
    return df[keep].copy()


def generate_report(maf_path: str, sample_id: str, tumor_type: str,
                    patient_meta: dict, out_path: str, min_alt_reads: int = 5):
    df = load_maf(maf_path)
    df = filter_clinically_relevant(df, min_alt_reads=min_alt_reads)
    df["Tier"] = df.apply(tier_variant, axis=1)
    df = df.sort_values("Tier")  # IA first

    # Build report
    lines = []
    lines.append("=" * 72)
    lines.append("PRECISION ONCOLOGY MOLECULAR PROFILING REPORT")
    lines.append("=" * 72)
    lines.append(f"Patient ID:        {patient_meta.get('patient_id', sample_id)}")
    lines.append(f"Sample ID:         {sample_id}")
    lines.append(f"Tumor type:        {tumor_type}")
    lines.append(f"Specimen:          {patient_meta.get('specimen', 'FFPE tumor + matched blood')}")
    lines.append(f"Report date:       {datetime.now().strftime('%Y-%m-%d')}")
    lines.append(f"Pipeline:          Parabricks 4.7 (GPU) -> Mutect2 -> VEP REST -> OncoKB")
    lines.append(f"Reference:         GRCh38 (Homo_sapiens_assembly38)")
    lines.append("")

    # Summary metrics
    lines.append("-" * 72)
    lines.append("QC SUMMARY")
    lines.append("-" * 72)
    lines.append(f"  Tumor mean coverage:   {patient_meta.get('tumor_cov', 'N/A')}x")
    lines.append(f"  Normal mean coverage:  {patient_meta.get('normal_cov', 'N/A')}x")
    lines.append(f"  Tumor purity (est.):   {patient_meta.get('purity', 'N/A')}")
    lines.append(f"  TMB (mut/Mb):          {patient_meta.get('tmb', 'N/A')}")
    lines.append(f"  MSI status:            {patient_meta.get('msi', 'N/A')}")
    lines.append("")

    # Tier I/II - actionable
    actionable = df[df["Tier"].isin(["Tier IA", "Tier IB", "Tier IIC", "Tier IID"])]
    lines.append("-" * 72)
    lines.append(f"ACTIONABLE VARIANTS  (n={len(actionable)})")
    lines.append("-" * 72)
    if len(actionable) == 0:
        lines.append("  No actionable variants detected.")
    for _, v in actionable.iterrows():
        raw = v.get('HGVSp_Short')
        hgvsp = '' if (raw is None or pd.isna(raw)) else str(raw)
        lines.append(f"\n  [{v['Tier']}] {v['Hugo_Symbol']} {hgvsp or '?'}")
        lines.append(f"      Consequence:   {v['Consequence']}")
        lines.append(f"      VAF:           {v['VAF']:.1%}  (depth: {int(v['t_alt_count'] + v['t_ref_count'])})")
        lines.append(f"      Position:      {v['Chromosome']}:{v['Start_Position']} {v['Reference_Allele']}>{v['Tumor_Seq_Allele2']}")
        if pd.notna(v.get("HIGHEST_LEVEL")) and str(v.get("HIGHEST_LEVEL", "")):
            lines.append(f"      OncoKB level:  {v['HIGHEST_LEVEL']}")
        if pd.notna(v.get("HIGHEST_LEVEL_SUMMARY")) and str(v.get("HIGHEST_LEVEL_SUMMARY", "")):
            lines.append(f"      Therapy:       {v['HIGHEST_LEVEL_SUMMARY']}")

    # Tier III - VUS
    vus = df[df["Tier"] == "Tier III"]
    lines.append("")
    lines.append("-" * 72)
    lines.append(f"VARIANTS OF UNKNOWN SIGNIFICANCE  (n={len(vus)})")
    lines.append("-" * 72)
    for _, v in vus.head(20).iterrows():
        raw = v.get('HGVSp_Short')
        hgvsp = '' if (raw is None or pd.isna(raw)) else str(raw)
        lines.append(f"  {str(v['Hugo_Symbol']):10s} {hgvsp:20s} "
                     f"VAF={v['VAF']:.1%}  {v['Consequence']}")

    lines.append("")
    lines.append("=" * 72)
    lines.append("DISCLAIMER: This report is intended for use by qualified medical")
    lines.append("professionals. Variant interpretation should be performed in the")
    lines.append("context of clinical presentation and additional testing.")
    lines.append("=" * 72)

    report_text = "\n".join(lines)
    Path(out_path).write_text(report_text)

    # Also emit machine-readable JSON for downstream EHR integration
    json_path = Path(out_path).with_suffix(".json")
    json_path.write_text(json.dumps({
        "sample_id": sample_id,
        "tumor_type": tumor_type,
        "report_date": datetime.now().isoformat(),
        "patient_meta": patient_meta,
        "actionable_variants": actionable.to_dict(orient="records"),
        "vus": vus.to_dict(orient="records"),
    }, default=str, indent=2))

    print(f"Report:  {out_path}")
    print(f"JSON:    {json_path}")
    return report_text


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--maf", required=True, help="OncoKB-annotated MAF file")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--tumor-type", default="LUAD", help="OncoTree code")
    ap.add_argument("--out", required=True, help="Output report path (.txt)")
    ap.add_argument("--meta-json", help="Patient metadata JSON file")
    ap.add_argument("--min-alt-reads", type=int, default=5,
                    help="Minimum alt read depth to report a variant (default 5; use 2 for low-depth amplicon smoke tests)")
    args = ap.parse_args()

    meta = {}
    if args.meta_json:
        meta = json.loads(Path(args.meta_json).read_text())

    generate_report(args.maf, args.sample_id, args.tumor_type, meta, args.out,
                    min_alt_reads=args.min_alt_reads)
