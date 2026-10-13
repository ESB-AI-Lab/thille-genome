#!/bin/bash
# Whole-genome stats for the existing 12-chromosome BAMs: samtools stats (independent check) and mapping_stats_full.py.
set -euo pipefail
export PATH="/gpool/bin/samtools-1.23.1/bin:$PATH"
O=summaries/mapping_stats/rerun_20261002; R=/epool/edwin/pamericana/avoflow/mappingback/minimap2
for n in thille.hifi.thille thille.hifi.hass_uoq; do
  samtools stats -@ 2 $R/$n.bam > $O/$n.chr12only.samtools_stats.txt
  samtools view -@ 2 $R/$n.bam | python3 scripts/mapping_stats_full.py ${n}_12chr > $O/$n.chr12only.full_stats.tsv
  echo "done $n $(date)"
done
