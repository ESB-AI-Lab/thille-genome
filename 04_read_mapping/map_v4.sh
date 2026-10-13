#!/bin/bash
# Align a read set to the Thille v4 assembly (Edwin, 2026-10-04) with the same minimap2 2.30 / samtools 1.23.1 binaries
# and options as avoflow/mappingback/sbatch_align_{hifi,clr}.sh: minimap2 -ax <preset> --secondary=no, sorted+indexed BAM,
# then per-chromosome BAMs that keep the full header (only for references of at most 20 sequences).
# Usage: scripts/map_v4.sh <name> <map-hifi|map-pb> <reads.fastq.gz> <reference.fasta> [threads]
set -euo pipefail
export PATH="/gpool/bin/minimap2-2.30_x64-linux:/gpool/bin/samtools-1.23.1/bin:/gpool/bin/samtools-1.23.1/htslib-1.23.1/bin:$PATH"
NAME="$1"; PRESET="$2"; READS="$3"; REF="$4"; T="${5:-17}"
OUT=/epool/edwin/pamericana/avoflow/mappingback/v4
mkdir -p "$OUT"; BAM="$OUT/$NAME.bam"
[ -e "$BAM" ] && { echo "exists: $BAM"; exit 1; }
echo "start $(date) host $(hostname) preset $PRESET ref $REF reads $READS threads $T"; minimap2 --version; samtools --version | head -1
minimap2 -ax "$PRESET" -t "$T" --secondary=no "$REF" "$READS" | samtools sort -@ 2 -m 6G -T "$OUT/$NAME.sorttmp" -o "$BAM.tmp" -
mv "$BAM.tmp" "$BAM"; samtools index -@ 4 "$BAM"
NSEQ=$(wc -l < "$REF.fai"); [ "$NSEQ" -le 20 ] && for c in $(cut -f1 "$REF.fai"); do samtools view -@ 2 -b -o "$OUT/$NAME.$c.bam" "$BAM" "$c" && samtools index "$OUT/$NAME.$c.bam"; done
touch "$BAM.done"; echo "done $(date)"
