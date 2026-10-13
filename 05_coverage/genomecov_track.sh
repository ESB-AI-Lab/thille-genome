#!/bin/bash
# Per-position coverage for one alignment track, made the same way as the Thille files
# (BEDTools v2.30.0 genomecov -ibam <bam> -d -split on each per-chromosome BAM). The 12 chromosomes run in parallel.
# Output: covfiles/<track>.<contig>.coverage.gz, contig names from the chromosome map; never overwrites an existing file.
#
# Usage (from the repo root):
#   scripts/genomecov_track.sh <track> <bam prefix> <chromosome map TSV>
#   e.g. scripts/genomecov_track.sh gwen.clr.wi /epool/edwin/pamericana/avoflow/mappingback/gwen.clr.wi wi_cm_chromosomes.tsv
#        (reads <bam prefix>.chr1.bam .. .chr12.bam; the map's column 2 is the contig name, column 3 "Chromosome N")
set -euo pipefail
TRACK="$1"; PREFIX="$2"; MAP="$3"
IMG=containers/aligners_v3.sif
run_one() {
    local bam="$1" out="$2"
    if [ -e "$out" ]; then echo "exists, skipped: $out"; return 0; fi
    # pipefail here: xargs runs this in a fresh bash, which does not inherit the script's set options, and an
    # interrupted bedtools must never leave a finished-looking file behind
    set -o pipefail
    if singularity exec --bind /epool "$IMG" bedtools genomecov -ibam "$bam" -d -split | gzip -c > "${out}.tmp"; then
        mv "${out}.tmp" "$out"
    else
        echo "FAILED: $out (partial output left as ${out}.tmp)" >&2
        return 1
    fi
    echo "done: $out"
}
export -f run_one; export IMG
for n in $(seq 1 12); do
    contig=$(awk -F'\t' -v n="$n" '$0 !~ /^#/ && $3 == "Chromosome " n {print $2}' "$MAP")
    [ -n "$contig" ] || { echo "no contig for chromosome $n in $MAP" >&2; exit 1; }
    [ -s "${PREFIX}.chr${n}.bam" ] || { echo "missing ${PREFIX}.chr${n}.bam" >&2; exit 1; }
    echo "${PREFIX}.chr${n}.bam covfiles/${TRACK}.${contig}.coverage.gz"
done | xargs -P 12 -L 1 bash -c 'run_one "$0" "$1"'
