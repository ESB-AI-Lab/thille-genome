#!/bin/bash
#SBATCH -J hic.cool.v4
#SBATCH -o logs/hic.cool.v4.o%A.%a
#SBATCH -e logs/hic.cool.v4.e%A.%a
#SBATCH -c 8
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-1
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Thille Hi-C: juicer-pre pairs (YaHS scaffold coords, 0-based) -> 50 kb .cool -> balanced .mcool -> figures.
# Array task 0 = chr12 (chr1-12 only), 1 = all (all 1,491 scaffolds), 2 = chr12scf13 (chr1-12 + scaffold_13;
# submit with --array=2).
# v4 copy (2026-10-04): pairs from juicer pre with the v4 AGP (already chr-named), outputs in hic/v4/.
# Submit from assemblies/thille.hic: sbatch hic/sbatch_cool_v4.sh [pairs.txt]
# No extra binds: inputs and outputs are under the auto-bound cwd.
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/hic/.mplconfig

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
PAIRS=${1:-hic/v4/thille.hifiv4.asm.hic.p_ctg.fasta_out_jbmap.txt}
SCOPES=(chr12 all)
SCOPE=${SCOPES[$SLURM_ARRAY_TASK_ID]}
BASE=50000
RES=50000,100000,250000,500000,1000000
OUT=hic/v4/$SCOPE
mkdir -p "$OUT" "$MPLCONFIGDIR"

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} scope=$SCOPE pairs=$PAIRS"
free -g
run() { singularity exec --no-home "$SIF" "$@"; }
run cooler --version

if [[ -f $OUT/$SCOPE.${BASE}.cool.done ]]; then
  echo "[$(date)] base cool exists with .done, skipping cload"
else
  echo "[$(date)] cload start"
  awk -v SIZES=hic/v4/sizes/$SCOPE.sizes -v RENAME=hic/v4/sizes/rename.tsv -v STATS=$OUT/$SCOPE.pairstats.tsv \
      -f hic/rename_pairs.awk "$PAIRS" \
    | run cooler cload pairs -0 -c1 2 -p1 3 -c2 6 -p2 7 --temp-dir "$OUT" \
        "hic/v4/sizes/$SCOPE.sizes:$BASE" - "$OUT/$SCOPE.$BASE.cool.tmp"
  mv "$OUT/$SCOPE.$BASE.cool.tmp" "$OUT/$SCOPE.$BASE.cool"
  date > "$OUT/$SCOPE.$BASE.cool.done"
  cat "$OUT/$SCOPE.pairstats.tsv"
fi

echo "[$(date)] zoomify + balance start"
run cooler zoomify --balance --balance-args "--nproc ${SLURM_CPUS_PER_TASK}" -p "${SLURM_CPUS_PER_TASK}" \
    -r "$RES" -o "$OUT/$SCOPE.mcool.tmp" "$OUT/$SCOPE.$BASE.cool"
mv "$OUT/$SCOPE.mcool.tmp" "$OUT/$SCOPE.mcool"

echo "[$(date)] plots"
run python3 hic/plot_contacts.py "$OUT/$SCOPE.mcool" "$SCOPE" hic/v4/figures "$OUT/$SCOPE.pairstats.tsv"
echo "[$(date)] done"
