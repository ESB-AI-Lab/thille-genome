#!/bin/bash
# Reduce one bedtools genomecov -d file (chrom, 1-based pos, depth; gzipped) to 10 kb bins of mean depth, plus the
# chromosome median depth over all positions (zero-depth included, as in the Methods). Streams the file once.
# The per-chromosome genomecov files list every position of the genome (the per-chromosome BAMs keep the full
# header), zero outside their own chromosome, so only rows of <chrom> are used.
# Usage: bin_genomecov.sh <in.coverage.gz> <chrom> <out.tsv> [bin=10000]
set -euo pipefail
IN=$1; CHROM=$2; OUT=$3; BIN=${4:-10000}
# awk exits once its chromosome ends, so zcat gets SIGPIPE (141); only awk's status counts.
set +o pipefail
zcat "$IN" | awk -v B="$BIN" -v C="$CHROM" -v OUT="$OUT.tmp" '
  $1 != C { if (seen) exit; next }
  { seen = 1; d = $3 + 0; h[d]++; n++; b = int(($2 - 1) / B); s[b] += d; c[b]++; if (b > mb) mb = b; chrom = $1 }
  END {
    half = (n + 1) / 2; cum = 0; med = 0
    PROCINFO["sorted_in"] = "@ind_num_asc"
    for (k in h) { cum += h[k]; if (cum >= half) { med = k; break } }
    printf "#chrom=%s positions=%d median_depth=%d bin=%d\n", chrom, n, med, B > OUT
    print "bin_start\tbin_end\tmean_depth" > OUT
    for (i = 0; i <= mb; i++) printf "%d\t%d\t%.2f\n", i * B, i * B + c[i], (c[i] ? s[i] / c[i] : 0) > OUT
  }'
[[ ${PIPESTATUS[1]} -eq 0 ]] || { echo "ERROR: awk failed on $IN" >&2; exit 1; }
mv "$OUT.tmp" "$OUT"
