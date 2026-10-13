#!/bin/bash
# samtools stats split by chromosome (plus a count of unplaced unmapped reads) for the four full BAMs (Table Z).
# Each primary read is placed on exactly one chromosome, so the per-chromosome SN totals sum to the whole-BAM totals.
# Usage (from the repo root): summaries/mapping_stats/run_stats_parallel.sh [parallel jobs, default 16]
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
for pair in "gwen.clr.thille:$R/gwen.clr.thille.bam" "thille.hifi.thille:$R/minimap2/thille.hifi.thille.bam" \
            "gwen.clr.hass_uoq:$R/gwen.clr.hass_uoq.bam" "thille.hifi.hass_uoq:$R/minimap2/thille.hifi.hass_uoq.bam"; do
    name=${pair%%:*}; bam=${pair#*:}
    for region in $(singularity exec --bind /epool "$IMG" samtools idxstats "$bam" | awk '$1 != "*" {print $1}') UNMAPPED; do
        echo "$name $bam $region"
    done
done
} | xargs -P "$JOBS" -L 1 bash -c 'one "$0" "$1" "$2"'
echo ALL_DONE
