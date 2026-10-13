#!/bin/bash
#SBATCH -J compare.gwen.v4
#SBATCH -o logs/compare.gwen.v4.o%A.%a
#SBATCH -e logs/compare.gwen.v4.e%A.%a
#SBATCH -c 8
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-23
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Thille chrN vs Gwen pseudo-chrN (--compare-only). Tasks 0-11 = asm5 chr1-12, 12-23 = asm10 chr1-12.
# Submit from assemblies/thille.hic after compare/sbatch_prep_v4.sh:
#   sbatch --dependency=afterok:<prep_jobid> compare/sbatch_cmp_gwen_v4.sh [tag]
# tag selects the placement set: gwen (default, approximate mapping) or gwen_c (base-level, -c).
# No extra binds: all inputs/outputs are under the auto-bound cwd.
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
TAG=${1:-gwen}
PRESETS=(asm5 asm10)
i=${SLURM_ARRAY_TASK_ID}
PRESET=${PRESETS[$((i / 12))]}
CHR=chr$((i % 12 + 1))
REF=compare/v4/thille_chrs/$CHR.fasta
QRY=compare/v4/$TAG/$PRESET/${TAG}_$PRESET.$CHR.fasta
OUT=compare/v4/$TAG/$PRESET/plots

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${i} tag=$TAG preset=$PRESET chr=$CHR"
free -g
for f in "$REF" "$QRY"; do
  [[ -s $f ]] || { echo "ERROR: missing input $f (no Gwen contigs assigned to $CHR?)" >&2; exit 1; }
done
mkdir -p "$OUT"
singularity exec --no-home "$SIF" moddotplot static -f "$REF" "$QRY" --compare-only -o "$OUT"
echo "[$(date)] done"
