#!/bin/bash
# Per-position coverage for the UoQ Hass alignments, made the same way as the Thille files
# (BEDTools v2.30.0 genomecov -ibam <bam> -d -split on each per-chromosome BAM, the command used for the Thille files).
# Output: covfiles/<track>.<CM accession>.coverage.gz ; never overwrites an existing file.
#
# Usage (from the repo root): scripts/uoq_genomecov.sh [parallel jobs, default 12]
set -euo pipefail
JOBS="${1:-12}"
R=/epool/edwin/pamericana/avoflow/mappingback
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
{
for n in $(seq 1 12); do
    cm=$(awk -F'\t' -v n="$n" '$0 !~ /^#/ && $3 == "Chromosome " n {print $2}' hass_cm_chromosomes.tsv)
    echo "$R/gwen.clr.hass_uoq.chr${n}.bam covfiles/gwen.clr.hass_uoq.${cm}.coverage.gz"
    echo "$R/minimap2/thille.hifi.hass_uoq.chr${n}.bam covfiles/thille.hifi.hass_uoq.${cm}.coverage.gz"
done
} | xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1"'
