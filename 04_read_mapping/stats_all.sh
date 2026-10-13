#!/bin/bash
# Whole-genome stats (samtools stats as the independent check, and mapping_stats_full.py) for the four Thille HiFi BAMs.
set -uo pipefail
export PATH="/gpool/bin/samtools-1.23.1/bin:$PATH"
O=summaries/mapping_stats/rerun_20261002; A=/epool/edwin/pamericana/avoflow/mappingback/minimap2; B=/epool/edwin/pamericana/avoflow/mappingback/rerun_20261002
run() { n=$1; bam=$2
  [ -s $O/$n.samtools_stats.txt ] || samtools stats -@ 2 $bam > $O/$n.samtools_stats.txt.tmp && mv $O/$n.samtools_stats.txt.tmp $O/$n.samtools_stats.txt 2>/dev/null
  samtools view -@ 2 $bam | python3 scripts/mapping_stats_full.py $n > $O/$n.full_stats.tsv.tmp && mv $O/$n.full_stats.tsv.tmp $O/$n.full_stats.tsv
  echo "done $n $(date)"; }
[ -s $O/thille.hifi.thille.chr12only.samtools_stats.txt ] && mv $O/thille.hifi.thille.chr12only.samtools_stats.txt $O/thille_12chr.samtools_stats.txt
run thille_12chr $A/thille.hifi.thille.bam &
run hass_uoq_12chr $A/thille.hifi.hass_uoq.bam &
run thille_13seq $B/thille.hifi.thille_chr_scf13.bam &
run hass_uoq_full $B/thille.hifi.hass_uoq_full.bam &
wait; echo ALL_DONE
