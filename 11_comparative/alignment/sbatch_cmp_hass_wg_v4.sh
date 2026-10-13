#!/bin/bash
#SBATCH -J compare.hass.wg.v4
#SBATCH -o logs/compare.hass.wg.v4.o%j
#SBATCH -e logs/compare.hass.wg.v4.e%j
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Whole-genome Thille vs Hass UoQ (chr1..12 concatenated, 100 kb N spacers), --compare-only.
# --mem is set because the default (2.5 GB/CPU) is far below the ~120 GB extrapolated from a
# measured 7.4 GB for chr12-vs-chr12 (109 Mb total).
# Submit from assemblies/thille.hic after compare/sbatch_prep_v4.sh:
#   sbatch --dependency=afterok:<prep_jobid> compare/sbatch_cmp_hass_wg_v4.sh [query.fa ref.fa outdir]
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
# Optional args (defaults = original run): query fasta, reference fasta, output dir.
QRY=${1:-compare/v4/wg/thille_wg.fasta}
REF=${2:-compare/v4/wg/hass_wg.fasta}
OUT=${3:-compare/v4/wg/plots}

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
free -g
mkdir -p "$OUT"
/usr/bin/time -v singularity exec --no-home "$SIF" \
  moddotplot static -f "$QRY" "$REF" --compare-only -o "$OUT"
echo "[$(date)] done"
