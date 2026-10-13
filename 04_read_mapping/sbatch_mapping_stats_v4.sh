#!/bin/bash
#SBATCH --job-name=mapping_stats.v4
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/mapping_stats_v4.%j.log
# Table 3 statistics for the two read set and reference pairs not yet computed: Puebla HiFi on Thille v4 (SLURM 13173)
# and Unnamed Lowland HiFi on Hass (UoQ, 12 chromosomes). Same three readers as the other rows: samtools stats (check),
# scripts/mapping_stats_full.py (mapping rate, MAPQ, NM error, QV) and scripts/identity_stats.py (BLAST identity).
# Writes summaries/v4/mapping_stats/<name>.{samtools_stats.txt,full_stats.tsv,identity.tsv}; never overwrites.
# Submit from the repo root: sbatch scripts/sbatch_mapping_stats_v4.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node $(hostname)  cpus ${SLURM_CPUS_PER_TASK:-?}  start $(date)"
cd /gpool/edwin/avocados/thille_paper
export PATH="/gpool/bin/samtools-1.23.1/bin:$PATH"
O=summaries/v4/mapping_stats
run() {  # name bam
    local n=$1 b=$2
    for f in $O/$n.samtools_stats.txt $O/$n.full_stats.tsv $O/$n.identity.tsv; do [ ! -e "$f" ] || { echo "exists: $f"; return 1; }; done
    samtools stats -@ 2 "$b" > $O/$n.samtools_stats.txt &
    samtools view -@ 1 "$b" | python3 scripts/mapping_stats_full.py "$n" > $O/$n.full_stats.tsv &
    samtools view -@ 1 -F 0x904 "$b" | python3 scripts/identity_stats.py "$n" > $O/$n.identity.tsv &
    wait
    echo "done $n $(date)"
}
run puebla.hifi.thille_v4 /epool/edwin/pamericana/avoflow/mappingback/puebla/minimap2/puebla.hifi.thille_v4.bam &
run westindian.hifi.hass_uoq /epool/edwin/pamericana/avoflow/mappingback/minimap2/westindian.hifi.hass_uoq.bam &
wait
grep -h -v "^label" $O/puebla.hifi.thille_v4.full_stats.tsv $O/westindian.hifi.hass_uoq.full_stats.tsv $O/puebla.hifi.thille_v4.identity.tsv $O/westindian.hifi.hass_uoq.identity.tsv
echo "end $(date)"
