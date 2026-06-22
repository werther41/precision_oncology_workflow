#!/usr/bin/awk -f
# Lifts FORMAT/DP -> INFO/TDP and FORMAT/VAF -> INFO/TVAF.
# Required by PCGR v2.2.5 which expects depth/AF as INFO (not FORMAT) fields.
# Usage: bcftools view input.vcf.gz | awk -f this_script | bgzip > out.vcf.gz
BEGIN { OFS = "\t"; dp_hdr = 0; vaf_hdr = 0 }

/^##FORMAT=<ID=DP,/ {
    if (!dp_hdr) {
        print "##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Tumor sequencing depth lifted from FORMAT/DP\">"
        dp_hdr = 1
    }
    print; next
}

/^##FORMAT=<ID=VAF,/ {
    if (!vaf_hdr) {
        print "##INFO=<ID=TVAF,Number=A,Type=Float,Description=\"Tumor variant allele fraction lifted from FORMAT/VAF\">"
        vaf_hdr = 1
    }
    print; next
}

/^#/ { print; next }

{
    n = split($9, fmt, ":")
    split($10, vals, ":")
    dp = ""; vaf = ""
    for (i = 1; i <= n; i++) {
        if (fmt[i] == "DP")  dp  = vals[i]
        if (fmt[i] == "VAF") vaf = vals[i]
    }
    tag = ""
    if (dp  != "" && dp  != ".") tag = tag (tag == "" ? "" : ";") "TDP="  dp
    if (vaf != "" && vaf != ".") tag = tag (tag == "" ? "" : ";") "TVAF=" vaf
    if ($8 == ".")  $8 = (tag == "" ? "." : tag)
    else            $8 = $8 (tag == "" ? "" : ";" tag)
    print
}
