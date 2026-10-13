#!/bin/bash
# samtools stats split by chromosome (plus a count of unplaced unmapped reads) for the West Indian alignments (Table Z).
# Each primary read is placed on exactly one chromosome, so the per-chromosome SN totals sum to the whole-BAM totals.
# Usage (from the repo root): summaries/mapping_stats/run_stats_wi.sh [parallel jobs, default 16]
set -euo pipefail
JOBS="${1:-16}"
R=/epool/edwin/pamericana/avoflow/mappingback
IMG=containers/aligners_v3.sif
OUT=summaries/mapping_stats/parts
export IMG OUT
one() {
    local name="$1" bam="$2" region="$3" out
    if [ "$region" = "UNMAPPED" ]; then
        out="$OUT/$name.unmapped_unplaced.count"
        [ -e "$out" ] && return 0
        singularity exec --bind /epool "$IMG" samtools view -c -f 4 "$bam" '*' > "$out.tmp" && mv "$out.tmp" "$out"
    else
        out="$OUT/$name.$region.stats.txt"
        [ -e "$out" ] && return 0
        singularity exec --bind /epool "$IMG" samtools stats "$bam" "$region" > "$out.tmp" && mv "$out.tmp" "$out"
    fi
    echo "done $name $region"
}
export -f one
{
for pair in "gwen.clr.wi:$R/gwen.clr.wi.bam" "thille.hifi.wi:$R/minimap2/thille.hifi.wi.bam"; do
    name=${pair%%:*}; bam=${pair#*:}
    for region in $(singularity exec --bind /epool "$IMG" samtools idxstats "$bam" | awk '$1 != "*" {print $1}') UNMAPPED; do
        echo "$name $bam $region"
    done
done
} | xargs -P "$JOBS" -L 1 bash -c 'one "$0" "$1" "$2"'
echo ALL_DONE
