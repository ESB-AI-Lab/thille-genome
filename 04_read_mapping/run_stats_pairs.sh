#!/bin/bash
# samtools stats split by chromosome (plus a count of unplaced unmapped reads) for the alignments given as name:bam arguments (Table Z).
# Each primary read is placed on exactly one chromosome, so the per-chromosome SN totals sum to the whole-BAM totals.
# Usage (from the repo root): JOBS=26 summaries/mapping_stats/run_stats_pairs.sh westindian.hifi.thille:/path/westindian.hifi.thille.bam [...]
set -euo pipefail
JOBS="${JOBS:-16}"
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
for pair in "$@"; do
    name=${pair%%:*}; bam=${pair#*:}
    for region in $(singularity exec --bind /epool "$IMG" samtools idxstats "$bam" | awk '$1 != "*" {print $1}') UNMAPPED; do
        echo "$name $bam $region"
    done
done
} | xargs -P "$JOBS" -L 1 bash -c 'one "$0" "$1" "$2"'
echo ALL_DONE
