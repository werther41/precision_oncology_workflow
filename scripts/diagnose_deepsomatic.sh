#!/usr/bin/env bash
# diagnose_deepsomatic.sh
# Root-cause why DeepSomatic shows ~30% precision on HCC1395 (its own training sample).
# Hypothesis under test: the matched normal isn't subtracting → germline floods the PASS set.
#
# Read-only. Runs decisive checks in order and prints a VERDICT per check + a final diagnosis.
# Requires: bcftools, samtools, tabix on PATH (or wrap in your gatk:4.5.0.0 container).

set -uo pipefail

# ---- config (override via env or args) ----
OUTDIR="${OUTDIR:-/mnt/storage/parabricks_test/output/HCC1395}"
RAW_VCF="${1:-${OUTDIR}/deepsomatic/HCC1395.deepsomatic.vcf.gz}"   # the UNFILTERED deepsomatic output
TUMOR_BAM="${2:-${OUTDIR}/HCC1395_tumor.bam}"
NORMAL_BAM="${3:-${OUTDIR}/HCC1395_normal.bam}"
GNOMAD="${GNOMAD:-/mnt/storage/parabricks_test/ref/af-only-gnomad.noalt.vcf.gz}"  # optional check 7

fail=0; warn=0
hr(){ printf '\n========== %s ==========\n' "$1"; }
verdict(){ printf '  >> VERDICT: %s\n' "$1"; }

# ---- Check 0: inputs exist ----
hr "0. Inputs"
for f in "$RAW_VCF" "$TUMOR_BAM" "$NORMAL_BAM"; do
  if [[ -s "$f" ]]; then echo "  ok   $f"; else echo "  MISSING $f"; fail=1; fi
done
[[ -f "${RAW_VCF}.tbi" || -f "${RAW_VCF}.csi" ]] || { echo "  no index for $RAW_VCF — run: tabix -p vcf $RAW_VCF"; }

# ---- Check 1: FILTER distribution (THE decisive diagnostic) ----
hr "1. FILTER distribution (raw DeepSomatic VCF)"
echo "  Expect: a LARGE GERMLINE bucket + a much smaller PASS bucket."
bcftools query -f '%FILTER\n' "$RAW_VCF" 2>/dev/null | sort | uniq -c | sort -rn
pass_n=$(bcftools view -H -f PASS "$RAW_VCF" 2>/dev/null | wc -l)
germ_n=$(bcftools view -H -i 'FILTER="GERMLINE"' "$RAW_VCF" 2>/dev/null | wc -l)
echo "  PASS=$pass_n   GERMLINE=$germ_n"
if [[ "$germ_n" -eq 0 ]]; then
  verdict "FAIL — zero GERMLINE calls. Somatic/germline classification never engaged (normal likely ignored)."; fail=1
elif [[ "$pass_n" -gt 60000 ]]; then
  verdict "FAIL — PASS count ($pass_n) is implausibly high for a somatic callset. Germline is leaking into PASS."; fail=1
else
  verdict "ok — PASS/GERMLINE split looks sane; problem is likely downstream (extraction or eval row)."
fi

# ---- Check 2: sample structure of the VCF ----
hr "2. VCF sample columns"
samples=$(bcftools query -l "$RAW_VCF" 2>/dev/null)
nsamp=$(echo "$samples" | grep -c . )
echo "  $nsamp sample(s): $(echo $samples | tr '\n' ' ')"
echo "  (DeepSomatic T/N is typically SINGLE-sample = tumor. If multi-sample, your eval needs --sample.)"

# ---- Check 3: was the NORMAL actually consumed? ----
hr "3. DeepSomatic command line / inputs (from VCF header)"
hdr=$(bcftools view -h "$RAW_VCF" 2>/dev/null)
echo "$hdr" | grep -iE 'deepsomatic|commandline|in-normal-bam|reads_normal|normal' | head -20
if echo "$hdr" | grep -qiE 'in-normal-bam|reads_normal|normal\.bam|HCC1395_normal'; then
  verdict "ok — normal BAM appears referenced in the run."
else
  verdict "WARN — no normal BAM reference found in header. Confirm the run wasn't tumor-only."; warn=1
fi

# ---- Check 4: BAM read-group SM tags distinct + not swapped/identical ----
hr "4. BAM SM tags + integrity"
sm_t=$(samtools view -H "$TUMOR_BAM"  2>/dev/null | sed -n 's/.*\tSM:\([^\t]*\).*/\1/p' | sort -u | tr '\n' ',')
sm_n=$(samtools view -H "$NORMAL_BAM" 2>/dev/null | sed -n 's/.*\tSM:\([^\t]*\).*/\1/p' | sort -u | tr '\n' ',')
echo "  tumor  SM: ${sm_t:-<none>}"
echo "  normal SM: ${sm_n:-<none>}"
samtools quickcheck "$TUMOR_BAM"  && echo "  tumor  quickcheck ok"  || { echo "  tumor  quickcheck FAILED";  fail=1; }
samtools quickcheck "$NORMAL_BAM" && echo "  normal quickcheck ok"  || { echo "  normal quickcheck FAILED"; fail=1; }
if [[ -n "$sm_t" && "$sm_t" == "$sm_n" ]]; then
  verdict "FAIL — tumor and normal share the SAME SM tag. DeepSomatic can't pair them correctly."; fail=1
elif [[ -z "$sm_t" || -z "$sm_n" ]]; then
  verdict "WARN — missing SM tag on one arm."; warn=1
else
  verdict "ok — distinct SM tags."
fi

# ---- Check 5: did the EVAL input wrongly include GERMLINE? ----
hr "5. What went into vcfeval (PASS-extraction sanity)"
PASS_VCF="${OUTDIR}/deepsomatic/HCC1395.deepsomatic.pass.vcf.gz"
if [[ -s "$PASS_VCF" ]]; then
  eval_n=$(bcftools view -H "$PASS_VCF" 2>/dev/null | wc -l)
  eval_germ=$(bcftools view -H -i 'FILTER="GERMLINE"' "$PASS_VCF" 2>/dev/null | wc -l)
  echo "  eval input records: $eval_n   of which FILTER=GERMLINE: $eval_germ"
  if [[ "$eval_germ" -gt 0 ]]; then
    verdict "FAIL — GERMLINE records present in the eval VCF. Extraction did not use 'bcftools view -f PASS'."; fail=1
  else
    verdict "ok — eval VCF is PASS-only."
  fi
else
  echo "  (no pass.vcf.gz at $PASS_VCF — skip)"
fi

# ---- Check 6 (optional): germline leakage signature via common gnomAD overlap ----
hr "6. Germline-leakage signature (PASS calls that are common population variants)"
if [[ -s "$GNOMAD" ]]; then
  tmp=$(mktemp -d)
  bcftools view -f PASS -v snps "$RAW_VCF" -Oz -o "$tmp/pass_snv.vcf.gz" 2>/dev/null
  tabix -p vcf "$tmp/pass_snv.vcf.gz" 2>/dev/null
  bcftools annotate -a "$GNOMAD" -c 'INFO/GNOMAD_AF:=INFO/AF' "$tmp/pass_snv.vcf.gz" 2>/dev/null \
    | bcftools view -H -i 'INFO/GNOMAD_AF>=0.01' 2>/dev/null | wc -l > "$tmp/common_n"
  total=$(bcftools view -H "$tmp/pass_snv.vcf.gz" 2>/dev/null | wc -l)
  common=$(cat "$tmp/common_n")
  rm -rf "$tmp"
  if [[ "$total" -gt 0 ]]; then
    pct=$(awk -v c="$common" -v t="$total" 'BEGIN{printf "%.1f", 100*c/t}')
    echo "  PASS SNVs that are common gnomAD variants (AF>=1%): $common / $total = ${pct}%"
    if awk -v c="$common" -v t="$total" 'BEGIN{exit !(t>0 && c/t>0.30)}'; then
      verdict "FAIL — >30% of PASS calls are common germline variants. Strong germline-leakage signal."; fail=1
    else
      verdict "ok — low common-variant contamination in PASS."
    fi
  fi
else
  echo "  (gnomAD file not found at $GNOMAD — skip optional check)"
fi

# ---- Final diagnosis ----
hr "DIAGNOSIS"
if [[ "$fail" -ne 0 ]]; then
  cat <<'EOF'
  ROOT CAUSE: germline is in the PASS set — DeepSomatic's normal-subtraction did not take effect.
  Most likely one of:
    (a) normal BAM not passed / run behaved tumor-only   -> see checks 1 & 3
    (b) tumor & normal share SM tag or are swapped/identical -> see check 4
    (c) PASS-extraction kept GERMLINE-tagged records      -> see check 5
  FIX, then re-run Phases 1-3:
    - re-issue `pbrun deepsomatic` with BOTH --in-tumor-bam and --in-normal-bam, distinct SM tags
    - extract calls with exactly: bcftools view -f PASS  (PASS == somatic in DeepSomatic)
  Do NOT bolt the Mutect2 POPAF/MMQ/TLOD cascade on top to "fix" precision — that masks the bug.
EOF
elif [[ "$warn" -ne 0 ]]; then
  echo "  No hard failure, but warnings above. Re-check the normal-arm wiring before trusting numbers."
else
  echo "  Checks pass. If precision is still ~30%, re-read the vcfeval summary: confirm you took the"
  echo "  weighted best-F row, not the Threshold=None baseline, and that GQ (not TLOD) was the score field."
fi
echo
echo "Reminder: even a clean HCC1395 number is training-data memorization. The held-out sample is the real test."