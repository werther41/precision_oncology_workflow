#!/usr/bin/env python3
"""
VCF -> MAF annotator using Ensembl VEP REST API + optional OncoKB REST.

Supports both Mutect2 paired-sample VCF and DeepSomatic single-sample VCF.
No local VEP install or cache required — uses the public REST endpoint.
OncoKB annotation runs if ONCOKB_TOKEN env var is set; otherwise falls back
to a built-in hotspot table covering the most common actionable variants.

Usage (Mutect2 paired):
  python3 vcf_to_maf_vep_rest.py \
      --vcf sample.somatic.vcf \
      --tumor-id SAMPLE_T --normal-id SAMPLE_N \
      --tumor-type LUAD \
      --out sample.maf

Usage (DeepSomatic single-sample):
  python3 vcf_to_maf_vep_rest.py \
      --vcf somatic.vcf \
      --tumor-id colo829_tumor \
      --tumor-type SKCM \
      --out colo829.maf
"""

import argparse, os, sys, time, gzip, threading
import requests
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

VEP_REST    = "https://rest.ensembl.org/vep/human/region"
ONCOKB_REST = "https://www.oncokb.org/api/v1"

VEP_BATCH_SIZE = 200   # max variants per VEP REST POST
VEP_RATE_DELAY = 0.1   # seconds between batch submissions (rate-limiter token bucket)

# Thread-local storage for per-thread HTTP sessions (avoids reconnecting per batch)
_thread_local = threading.local()

CONSEQUENCE_RANK = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation",
    "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "feature_elongation", "regulatory_region_variant", "feature_truncation",
    "intergenic_variant",
]

CONSEQUENCE_TO_CLASS = {
    "missense_variant":           "Missense_Mutation",
    "stop_gained":                "Nonsense_Mutation",
    "stop_lost":                  "Nonstop_Mutation",
    "frameshift_variant":         "Frame_Shift_Del",
    "inframe_insertion":          "In_Frame_Ins",
    "inframe_deletion":           "In_Frame_Del",
    "splice_acceptor_variant":    "Splice_Site",
    "splice_donor_variant":       "Splice_Site",
    "splice_region_variant":      "Splice_Region",
    "synonymous_variant":         "Silent",
    "start_lost":                 "Translation_Start_Site",
    "protein_altering_variant":   "Missense_Mutation",
    "5_prime_UTR_variant":        "5'UTR",
    "3_prime_UTR_variant":        "3'UTR",
    "intron_variant":             "Intron",
    "upstream_gene_variant":      "5'Flank",
    "downstream_gene_variant":    "3'Flank",
    "intergenic_variant":         "IGR",
}

HOTSPOT_DB = {
    ("EGFR",  "p.L858R"):   ("LEVEL_1",  "Osimertinib, Erlotinib, Gefitinib (FDA-approved)"),
    ("EGFR",  "p.T790M"):   ("LEVEL_1",  "Osimertinib (FDA-approved, acquired resistance)"),
    ("EGFR",  "p.E746_A750del"): ("LEVEL_1", "Osimertinib, Afatinib (FDA-approved)"),
    ("KRAS",  "p.G12C"):    ("LEVEL_1",  "Sotorasib, Adagrasib (FDA-approved)"),
    ("KRAS",  "p.G12D"):    ("LEVEL_3B", "RAS inhibitors in clinical trials"),
    ("BRAF",  "p.V600E"):   ("LEVEL_1",  "Dabrafenib+Trametinib, Vemurafenib (FDA-approved)"),
    ("BRAF",  "p.V600K"):   ("LEVEL_1",  "Dabrafenib+Trametinib (FDA-approved, melanoma)"),
    ("ALK",   "p.F1174L"):  ("LEVEL_1",  "Alectinib, Crizotinib (FDA-approved)"),
    ("PIK3CA","p.H1047R"):  ("LEVEL_1",  "Alpelisib+Fulvestrant (FDA-approved, HR+HER2-)"),
    ("PIK3CA","p.E545K"):   ("LEVEL_2",  "PIK3CA inhibitors (evidence in HR+ breast)"),
    ("PIK3CA","p.E542K"):   ("LEVEL_2",  "PIK3CA inhibitors (evidence in HR+ breast)"),
    ("ERBB2", "p.S310F"):   ("LEVEL_2",  "Trastuzumab deruxtecan (FDA-approved HER2-mut)"),
    ("MET",   "p.Y1253D"):  ("LEVEL_2",  "Tepotinib, Capmatinib (MET exon14 skip)"),
    ("RET",   "p.M918T"):   ("LEVEL_1",  "Selpercatinib, Pralsetinib (FDA-approved)"),
    ("TP53",  "p.R248Q"):   ("LEVEL_3B", "APR-246 (clinical trials, TP53 reactivation)"),
    ("TP53",  "p.R248W"):   ("LEVEL_3B", "APR-246 (clinical trials, TP53 reactivation)"),
    ("TP53",  "p.R175H"):   ("LEVEL_3B", "APR-246 (clinical trials, TP53 reactivation)"),
    ("BRCA1", "p.M1775R"):  ("LEVEL_1",  "Olaparib, Rucaparib (FDA-approved, BRCA1 mut)"),
    ("BRCA2", "p.T2722R"):  ("LEVEL_1",  "Olaparib, Niraparib (FDA-approved, BRCA2 mut)"),
    ("NRAS",  "p.Q61K"):    ("LEVEL_3B", "MEK inhibitors (SKCM, clinical evidence)"),
    ("NRAS",  "p.Q61R"):    ("LEVEL_3B", "MEK inhibitors (SKCM, clinical evidence)"),
    ("NRAS",  "p.Q61L"):    ("LEVEL_3B", "MEK inhibitors (SKCM, clinical evidence)"),
    ("CDK4",  "p.R24C"):    ("LEVEL_3B", "CDK4/6 inhibitors (clinical evidence, melanoma)"),
}


def parse_vcf(vcf_path, tumor_id, normal_id=None, pass_only=True,
              min_vaf=0.05, min_alt_depth=5):
    """Parse VCF and return list of variant dicts.

    Handles:
    - Mutect2 paired VCF (two sample columns, FORMAT=GT:AD:AF:DP:...)
    - DeepSomatic single-sample VCF (FORMAT=GT:GQ:DP:AD:VAF:PL)
    """
    variants = []
    tumor_col = None
    normal_col = None
    single_sample = (normal_id is None)

    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.strip().split("\t")
                samples = cols[9:]
                for i, s in enumerate(samples):
                    if s == tumor_id:
                        tumor_col = i
                    elif normal_id and s == normal_id:
                        normal_col = i
                if tumor_col is None:
                    if single_sample:
                        tumor_col = 0
                    else:
                        # Mutect2 fallback: last = tumor
                        tumor_col = len(samples) - 1
                        normal_col = 0
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, pos, vid = parts[0], int(parts[1]), parts[2]
            ref, alt = parts[3], parts[4]
            filt = parts[6]

            # Skip multi-allelic — only take first ALT
            if "," in alt:
                alt = alt.split(",")[0]

            if pass_only and filt != "PASS":
                continue

            # Skip indels that are too complex for VEP region format
            # (VEP handles simple SNVs/indels fine)

            fmt = parts[8].split(":")

            def parse_sample(col):
                raw = parts[9 + col].split(":")
                return dict(zip(fmt, raw))

            tumor_vals = parse_sample(tumor_col)

            # VAF: DeepSomatic has a direct VAF field; Mutect2 uses AF
            if "VAF" in tumor_vals:
                try:
                    vaf = float(tumor_vals["VAF"])
                except (ValueError, TypeError):
                    vaf = 0.0
            elif "AF" in tumor_vals:
                try:
                    af_raw = tumor_vals["AF"].split(",")[0]
                    vaf = float(af_raw)
                except (ValueError, TypeError):
                    vaf = 0.0
            else:
                vaf = 0.0

            # Allele depths
            def ad_counts(vals):
                ad = vals.get("AD", "0,0").split(",")
                try:
                    return int(ad[0]), int(ad[1])
                except (ValueError, IndexError):
                    return 0, 0

            t_ref, t_alt = ad_counts(tumor_vals)

            # Normal counts (0 for single-sample)
            n_ref, n_alt = 0, 0
            if not single_sample and normal_col is not None:
                normal_vals = parse_sample(normal_col)
                n_ref, n_alt = ad_counts(normal_vals)

            # Apply quality filters
            if vaf < min_vaf:
                continue
            if t_alt < min_alt_depth:
                continue

            variants.append({
                "Chromosome":       chrom,
                "Start_Position":   pos,
                "Reference_Allele": ref,
                "Tumor_Seq_Allele2": alt,
                "t_ref_count":      t_ref,
                "t_alt_count":      t_alt,
                "n_ref_count":      n_ref,
                "n_alt_count":      n_alt,
                "t_vaf":            vaf,
                "vcf_id":           vid,
            })
    return variants


def _get_session():
    """Return a thread-local requests.Session, creating it on first use."""
    if not hasattr(_thread_local, "session"):
        s = requests.Session()
        s.headers.update({"Content-Type": "application/json", "Accept": "application/json"})
        _thread_local.session = s
    return _thread_local.session


def _vep_post_batch(region_list):
    """POST a single batch of ≤200 variants to VEP REST. Returns list of hits."""
    payload = {"variants": region_list}
    params  = {"hgvs": 1, "canonical": 1, "clinical_significance": 1, "af_gnomade": 1}
    session = _get_session()

    for attempt in range(4):
        try:
            r = session.post(VEP_REST, json=payload, params=params, timeout=90)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 60))
                print(f"[vep]   Rate-limited — sleeping {wait}s...", flush=True)
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except requests.exceptions.RequestException as e:
            if attempt < 3:
                time.sleep(2 ** attempt)
            else:
                raise RuntimeError(f"VEP REST failed after 4 attempts: {e}")
    return []


def _hits_to_maps(hits):
    """Convert a list of VEP hits into (exact_map, by_pos_map) fragments."""
    exact = {}
    by_pos = {}
    for hit in hits:
        chrom  = str(hit.get("seq_region_name", ""))
        pos    = int(hit.get("start", 0))
        allele = hit.get("allele_string", "/")
        parts  = allele.split("/")
        ref    = parts[0] if parts else ""
        alt    = parts[-1] if len(parts) > 1 else ""
        exact[(chrom, pos, ref, alt)] = hit
        by_pos.setdefault((chrom, pos), hit)
    return exact, by_pos


def vep_annotate_all(variants, workers=1):
    """Annotate all variants via VEP REST in batches of VEP_BATCH_SIZE.

    workers > 1 sends multiple batches concurrently — check Ensembl fair-use
    policy before using high values; 3–5 is typically safe for a single client.

    Returns two lookup dicts:
      - exact: (chrom, pos, ref, alt) -> vep_hit
      - by_pos: (chrom, pos) -> vep_hit  (fallback for allele_string mismatches)
    """
    def region_str(v):
        chrom = v["Chromosome"].replace("chr", "")
        ref, alt = v["Reference_Allele"], v["Tumor_Seq_Allele2"]
        pos = v["Start_Position"]
        vid = v["vcf_id"] or "."
        return f"{chrom} {pos} {vid} {ref} {alt} . ."

    total     = len(variants)
    n_batches = (total + VEP_BATCH_SIZE - 1) // VEP_BATCH_SIZE
    exact_map  = {}
    by_pos_map = {}
    lock       = threading.Lock()
    completed  = [0]
    t_start    = time.time()

    def run_batch(batch_i):
        batch   = variants[batch_i * VEP_BATCH_SIZE:(batch_i + 1) * VEP_BATCH_SIZE]
        regions = [region_str(v) for v in batch]
        time.sleep(VEP_RATE_DELAY * (batch_i % max(workers, 1)))  # stagger starts within each round
        try:
            return batch_i, _vep_post_batch(regions)
        except RuntimeError as e:
            print(f"[vep] WARN: {e} — skipping batch {batch_i+1}", flush=True)
            return batch_i, []

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_batch, i): i for i in range(n_batches)}
        for fut in as_completed(futures):
            batch_i, hits = fut.result()
            e_frag, p_frag = _hits_to_maps(hits)
            with lock:
                exact_map.update(e_frag)
                for k, v in p_frag.items():
                    by_pos_map.setdefault(k, v)
                completed[0] += 1
                done = min((batch_i + 1) * VEP_BATCH_SIZE, total)
                elapsed = time.time() - t_start
                rate = completed[0] / elapsed if elapsed > 0 else 0
                eta  = (n_batches - completed[0]) / rate if rate > 0 else 0
                print(f"[vep] Batch {batch_i+1}/{n_batches}  "
                      f"({done}/{total} variants, {100*done//total}%)"
                      f"  ETA {eta/60:.1f}min", flush=True)

    return exact_map, by_pos_map


def _lookup_vep(variant, exact_map, by_pos_map):
    """Find the VEP hit for a variant; fall back to position-only match.

    VEP normalizes indel coordinates: it strips the anchor base and advances
    the position by the anchor length.  For example:
      VCF deletion  43046406 TC→T   →  VEP start=43046407  allele_string=C/-
      VCF insertion 43046406 T→TCA  →  VEP start=43046406  allele_string=-/CA

    So we try both VCF and VEP-normalized keys before falling back to by_pos.
    """
    chrom = variant["Chromosome"].replace("chr", "")
    pos   = variant["Start_Position"]
    ref   = variant["Reference_Allele"]
    alt   = variant["Tumor_Seq_Allele2"]

    hit = exact_map.get((chrom, pos, ref, alt))
    if hit:
        return hit

    if len(ref) > len(alt):
        # Deletion: VEP trims anchor bases, pos advances by len(alt)
        vep_pos = pos + len(alt)
        vep_ref = ref[len(alt):]
        hit = (exact_map.get((chrom, vep_pos, vep_ref, "-"))
               or by_pos_map.get((chrom, vep_pos)))
        if hit:
            return hit
    elif len(alt) > len(ref):
        # Insertion: VEP uses anchor pos, allele "-/inserted_bases"
        vep_pos = pos + len(ref) - 1
        vep_alt = alt[len(ref):]
        hit = (exact_map.get((chrom, vep_pos, "-", vep_alt))
               or by_pos_map.get((chrom, vep_pos)))
        if hit:
            return hit

    return by_pos_map.get((chrom, pos))


AA3_TO_1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V","Ter":"*",
}

def shorten_hgvsp(hgvsp):
    """Convert 'ENSP...:p.Arg248Trp' -> 'p.R248W'."""
    if not hgvsp:
        return ""
    if ":" in hgvsp:
        hgvsp = hgvsp.split(":")[1]
    import re
    def repl(m):
        return AA3_TO_1.get(m.group(0), m.group(0))
    return re.sub(r"[A-Z][a-z]{2}", repl, hgvsp)


def pick_worst_transcript(vep_result):
    """Pick canonical transcript first, then most severe consequence."""
    transcripts = vep_result.get("transcript_consequences", [])
    if not transcripts:
        return {}, []

    def severity(tc):
        consequences = tc.get("consequence_terms", [])
        sev = min((CONSEQUENCE_RANK.index(c) if c in CONSEQUENCE_RANK else 999)
                  for c in consequences) if consequences else 999
        canonical_bonus = 0 if tc.get("canonical") == 1 else 1
        return (canonical_bonus, sev)

    worst = min(transcripts, key=severity)
    return worst, worst.get("consequence_terms", [])


def build_maf_row(variant, vep_hit, worst_tc, consequences):
    """Merge raw variant data with VEP annotation into a MAF row dict."""
    gene = worst_tc.get("gene_symbol") or vep_hit.get("gene_id", "UNKNOWN")

    top_consequence = consequences[0] if consequences else "intergenic_variant"
    variant_class   = CONSEQUENCE_TO_CLASS.get(top_consequence, "Targeted_Region")

    hgvsp_short = shorten_hgvsp(worst_tc.get("hgvsp", "") or "")

    gnomad_af = float(worst_tc.get("gnomade_af", 0) or 0)
    clin_sig_terms = []
    for cv in vep_hit.get("colocated_variants", []):
        freqs = cv.get("frequencies", {})
        gnomad_af = max(gnomad_af,
                        float(freqs.get("gnomade", {}).get("af", 0) or 0),
                        float(freqs.get("af", 0) or 0))
        if cv.get("clin_sig"):
            clin_sig_terms.extend(cv["clin_sig"])

    return {
        "Hugo_Symbol":            gene,
        "Chromosome":             variant["Chromosome"],
        "Start_Position":         variant["Start_Position"],
        "Reference_Allele":       variant["Reference_Allele"],
        "Tumor_Seq_Allele2":      variant["Tumor_Seq_Allele2"],
        "Variant_Classification": variant_class,
        "Consequence":            top_consequence,
        "HGVSp_Short":            hgvsp_short,
        "t_ref_count":            variant["t_ref_count"],
        "t_alt_count":            variant["t_alt_count"],
        "t_vaf":                  variant.get("t_vaf", 0),
        "n_ref_count":            variant["n_ref_count"],
        "n_alt_count":            variant["n_alt_count"],
        "gnomAD_AF":              gnomad_af,
        "CLIN_SIG":               ";".join(sorted(set(clin_sig_terms))),
        "HIGHEST_LEVEL":          "",
        "HIGHEST_LEVEL_SUMMARY":  "",
    }


def oncokb_annotate_rows(rows, tumor_type, token, workers=4):
    """Annotate via OncoKB REST API. Mutates rows in place. Uses concurrent requests."""
    mutations = [
        {"idx": i, "gene": r["Hugo_Symbol"], "alteration": r["HGVSp_Short"],
         "tumorType": tumor_type}
        for i, r in enumerate(rows) if r["HGVSp_Short"]
    ]
    if not mutations:
        return

    url = f"{ONCOKB_REST}/annotate/mutations/byHGVSp"

    def fetch_one(m):
        session = _get_session()
        session.headers["Authorization"] = f"Bearer {token}"
        try:
            resp = session.get(
                url,
                params={"hugoSymbol": m["gene"], "alteration": m["alteration"],
                        "tumorType": m["tumorType"]},
                timeout=15,
            )
            if resp.status_code == 200:
                return m["idx"], resp.json()
        except Exception:
            pass
        return m["idx"], None

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for idx, data in pool.map(fetch_one, mutations):
            if data is None:
                continue
            rows[idx]["HIGHEST_LEVEL"] = (data.get("highestSensitiveLevel")
                                          or data.get("highestLevel") or "")
            treatments = data.get("treatments", [])
            if treatments:
                rows[idx]["HIGHEST_LEVEL_SUMMARY"] = "; ".join(
                    t.get("drugs", [{}])[0].get("drugName", "")
                    for t in treatments[:3] if t.get("drugs")
                )


LOF_CONSEQUENCES = {
    "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant",
}

# Gene-level LoF rules applied when exact HGVSp not in HOTSPOT_DB.
# Clinically, PARP inhibitor eligibility for BRCA1/2 is LoF-agnostic.
GENE_LOF_DB = {
    "BRCA1": ("LEVEL_1", "Olaparib, Rucaparib, Niraparib (FDA-approved, BRCA1 LoF)"),
    "BRCA2": ("LEVEL_1", "Olaparib, Rucaparib, Niraparib (FDA-approved, BRCA2 LoF)"),
    "PALB2": ("LEVEL_2", "Olaparib (clinical evidence, PALB2 LoF, HR+ breast)"),
}


def hotspot_annotate_rows(rows):
    """Fallback: annotate using built-in hotspot DB when OncoKB token is absent."""
    for row in rows:
        key = (row["Hugo_Symbol"], row["HGVSp_Short"])
        if key in HOTSPOT_DB:
            row["HIGHEST_LEVEL"], row["HIGHEST_LEVEL_SUMMARY"] = HOTSPOT_DB[key]
        elif row["Hugo_Symbol"] in GENE_LOF_DB:
            # Gene-level LoF: any frameshift/nonsense/splice in BRCA1/BRCA2/PALB2
            csq = str(row.get("Consequence", ""))
            if any(c in csq for c in LOF_CONSEQUENCES):
                row["HIGHEST_LEVEL"], row["HIGHEST_LEVEL_SUMMARY"] = GENE_LOF_DB[row["Hugo_Symbol"]]


def main():
    ap = argparse.ArgumentParser(description="VCF -> MAF via Ensembl VEP REST + OncoKB")
    ap.add_argument("--vcf",           required=True, help="Input somatic VCF")
    ap.add_argument("--tumor-id",      required=True, help="Tumor sample name in VCF")
    ap.add_argument("--normal-id",     default=None,  help="Normal sample name (omit for single-sample VCF)")
    ap.add_argument("--tumor-type",    default="LUAD", help="OncoTree tumor type code")
    ap.add_argument("--out",           required=True, help="Output MAF path")
    ap.add_argument("--no-pass-filter", action="store_true",
                    help="Include non-PASS variants (default: PASS only)")
    ap.add_argument("--min-vaf",       type=float, default=0.05,
                    help="Minimum tumor VAF (default 0.05)")
    ap.add_argument("--min-alt-depth", type=int,   default=5,
                    help="Minimum tumor alt read depth (default 5)")
    ap.add_argument("--max-variants",  type=int,   default=0,
                    help="Cap variants for testing (0 = no cap)")
    ap.add_argument("--workers",       type=int,   default=5,
                    help="Concurrent VEP batch workers (default 5; use 1 for max politeness)")
    args = ap.parse_args()

    pass_only = not args.no_pass_filter
    print(f"[annotate] Parsing VCF: {args.vcf}")
    print(f"[annotate] Tumor: {args.tumor_id}  Normal: {args.normal_id or '(none — single-sample)'}")
    print(f"[annotate] Filters: PASS={pass_only}, VAF≥{args.min_vaf}, alt_depth≥{args.min_alt_depth}")

    variants = parse_vcf(
        args.vcf, args.tumor_id, args.normal_id,
        pass_only=pass_only,
        min_vaf=args.min_vaf,
        min_alt_depth=args.min_alt_depth,
    )
    print(f"[annotate] {len(variants)} variants after filtering")

    if args.max_variants and len(variants) > args.max_variants:
        print(f"[annotate] Capping at {args.max_variants} variants (--max-variants)")
        variants = variants[:args.max_variants]

    if not variants:
        print("[annotate] No variants found — writing empty MAF.")
        pd.DataFrame(columns=[
            "Hugo_Symbol","Chromosome","Start_Position","Reference_Allele",
            "Tumor_Seq_Allele2","Variant_Classification","Consequence","HGVSp_Short",
            "t_ref_count","t_alt_count","t_vaf","n_ref_count","n_alt_count",
            "gnomAD_AF","CLIN_SIG","HIGHEST_LEVEL","HIGHEST_LEVEL_SUMMARY",
        ]).to_csv(args.out, sep="\t", index=False)
        return

    # VEP annotation (batched, optionally concurrent)
    print(f"[annotate] Querying Ensembl VEP REST in batches of {VEP_BATCH_SIZE} "
          f"({args.workers} worker{'s' if args.workers > 1 else ''})...")
    t0 = time.time()
    exact_map, by_pos_map = vep_annotate_all(variants, workers=args.workers)
    elapsed = time.time() - t0
    print(f"[annotate] VEP done in {elapsed:.0f}s — {len(exact_map)} exact hits")

    rows = []
    unannotated = 0
    for variant in variants:
        vep_hit = _lookup_vep(variant, exact_map, by_pos_map)
        if vep_hit is None:
            unannotated += 1
            vep_hit = {}
        worst_tc, consequences = pick_worst_transcript(vep_hit)
        row = build_maf_row(variant, vep_hit, worst_tc, consequences)
        rows.append(row)

    if unannotated:
        print(f"[annotate] WARN: {unannotated} variants had no VEP hit")

    # OncoKB / hotspot annotation
    token = os.environ.get("ONCOKB_TOKEN", "")
    if token:
        print(f"[annotate] OncoKB annotation (tumor_type={args.tumor_type})...")
        oncokb_annotate_rows(rows, args.tumor_type, token, workers=args.workers)
    else:
        print("[annotate] ONCOKB_TOKEN not set — using built-in hotspot database")
        hotspot_annotate_rows(rows)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[annotate] MAF written: {args.out}  ({len(df)} variants)")

    actionable = df[df["HIGHEST_LEVEL"] != ""].shape[0]
    tier1 = df[df["HIGHEST_LEVEL"].isin(["LEVEL_1","LEVEL_2"])].shape[0]
    print(f"[annotate] Actionable variants: {actionable}  (LEVEL_1/2: {tier1})")

    # Print top hits
    if actionable:
        print("\n[annotate] Top actionable variants:")
        top = df[df["HIGHEST_LEVEL"] != ""].sort_values("HIGHEST_LEVEL").head(10)
        for _, r in top.iterrows():
            vaf_pct = f"{r['t_vaf']*100:.1f}%" if r['t_vaf'] > 0 else ""
            print(f"  {r['Hugo_Symbol']:12s} {r['HGVSp_Short'] or r['Consequence']:20s} "
                  f"VAF={vaf_pct:6s}  {r['HIGHEST_LEVEL']:10s}  {r['HIGHEST_LEVEL_SUMMARY'][:60]}")


if __name__ == "__main__":
    main()
