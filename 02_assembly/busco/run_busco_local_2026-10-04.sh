#!/bin/bash
# Local (non-SLURM) run of the two jobfile.txt entries, same image and output layout as sbatch_busco.sh.
set -euo pipefail
cd /epool/edwin/pamericana/avoflow/annotation/tiberius/thille.hic
SIF=/gpool/bin/singularityimages/busco6.1.sif
DL=/dpool/github/HapSolo/busco_downloads
export SINGULARITYENV_TINI_SUBREAPER=1
mkdir -p busco
while read ASM LINEAGE TYPE; do
  OUT="$(basename "$ASM" .fasta)_${LINEAGE}"
  ( echo "$(date) START $OUT ($TYPE)"
    singularity exec -B /epool -B /dpool "$SIF" busco -i "$ASM" -l "$LINEAGE" -o "$OUT" --out_path busco \
        -m "$TYPE" -c 20 --download_path "$DL" --offline
    echo "$(date) DONE $OUT"
    grep -E "C:" "busco/$OUT/short_summary.specific.${LINEAGE}.${OUT}.txt" || true
  ) > "busco/${OUT}.local.log" 2>&1 &
done < jobfile.txt
wait
