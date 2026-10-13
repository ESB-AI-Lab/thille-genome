#!/bin/bash
#SBATCH -J hic.juicer.v4
#SBATCH -o logs/hic.juicer.v4.o%A.%a
#SBATCH -e logs/hic.juicer.v4.e%A.%a
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-1
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Thille Hi-C .hic for Juicebox via juicer_tools 2.20.00 pre, following the YaHS README recipe
# (sort -k2,2d -k6,6d, then pre with a chrom sizes file). Array task 0 = chr12, 1 = all.
# --mem is set because sort (-S 24G) plus the JVM (-Xmx32G) exceed the 2.5 GB/CPU default.
# v4 copy (2026-10-05): v4 pairs and sizes, outputs hic/v4/<scope>/<scope>.hic.
# Submit from assemblies/thille.hic: sbatch hic/sbatch_hic_v4.sh [pairs.txt]
set -euo pipefail

JUICER=/gpool/bin/juicebox/juicer_tools.2.20.00.jar
PAIRS=${1:-hic/v4/thille.hifiv4.asm.hic.p_ctg.fasta_out_jbmap.txt}
SCOPES=(chr12 all)
SCOPE=${SCOPES[$SLURM_ARRAY_TASK_ID]}
OUT=hic/v4/$SCOPE
TMP=$OUT/tmp_hic_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$OUT" "$TMP"

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} scope=$SCOPE"
java -version 2>&1 | head -1
free -g

SORTED=$OUT/$SCOPE.pairs.sorted.txt
if [[ -f $SORTED.done ]]; then
  echo "[$(date)] sorted pairs exist with .done, skipping sort"
else
  echo "[$(date)] rename + sort start"
  awk -v SIZES=hic/v4/sizes/$SCOPE.sizes -v RENAME=hic/v4/sizes/rename.tsv -v STATS=$TMP/pairstats.tsv \
      -f hic/rename_pairs.awk "$PAIRS" \
    | LC_ALL=C sort -k2,2d -k6,6d -T "$TMP" --parallel="${SLURM_CPUS_PER_TASK}" -S 24G > "$SORTED.tmp"
  mv "$SORTED.tmp" "$SORTED"
  date > "$SORTED.done"
fi

echo "[$(date)] juicer_tools pre start"
java -Xmx32G -jar "$JUICER" pre -j "${SLURM_CPUS_PER_TASK}" -t "$TMP" \
    "$SORTED" "$OUT/$SCOPE.hic.tmp" "hic/v4/sizes/$SCOPE.sizes"
mv "$OUT/$SCOPE.hic.tmp" "$OUT/$SCOPE.hic"
echo "[$(date)] done: $OUT/$SCOPE.hic"
