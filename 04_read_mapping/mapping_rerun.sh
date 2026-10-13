#!/bin/bash
# Re-map Thille HiFi reads (Edwin, 2026-10-02) with the same minimap2 2.30 / samtools 1.23.1 binaries and options as
# avoflow/mappingback/sbatch_align_hifi.sh: minimap2 -ax map-hifi --secondary=no, sorted and indexed BAM.
# Usage: scripts/mapping_rerun.sh <name> <reference.fasta> [map threads]
set -euo pipefail
export PATH="/gpool/bin/minimap2-2.30_x64-linux:/gpool/bin/samtools-1.23.1/bin:/gpool/bin/samtools-1.23.1/htslib-1.23.1/bin:$PATH"
NAME="$1"; REF="$2"; T="${3:-17}"
READS=/epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz
OUT=/epool/edwin/pamericana/avoflow/mappingback/rerun_20261002
mkdir -p "$OUT"; BAM="$OUT/$NAME.bam"
[ -e "$BAM" ] && { echo "exists: $BAM"; exit 1; }
echo "start $(date) host $(hostname) ref $REF threads $T"; minimap2 --version; samtools --version | head -1
minimap2 -ax map-hifi -t "$T" --secondary=no "$REF" "$READS" \
  | samtools sort -@ 2 -m 6G -T "$OUT/$NAME.sorttmp" -o "$BAM.tmp" -
mv "$BAM.tmp" "$BAM"; samtools index -@ 4 "$BAM"; touch "$BAM.done"
echo "done $(date)"
