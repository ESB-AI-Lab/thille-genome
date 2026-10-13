#!/bin/bash
#SBATCH --job-name=bin_puebla
#SBATCH --partition=gcbigmem,gcluster
#SBATCH --cpus-per-task=2
#SBATCH --array=0-1
#SBATCH --output=compare/v4/regions/cov/bin_puebla.%A_%a.log
# Puebla HiFi on Thille v4 binned to 10 kb for the chr3 and chr10 region track figures (audit A3, 2026-10-07), with
# the same bin_genomecov.sh used for the other read sets. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
CHR=(chr3 chr10); C=${CHR[$SLURM_ARRAY_TASK_ID]}
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} chrom=$C"
OUT=compare/v4/regions/cov/puebla.hifi.thille_v4.$C.10kb.tsv
[[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
bash compare/regions/cov/bin_genomecov.sh /gpool/edwin/avocados/thille_paper/covfiles/puebla.hifi.thille_v4.$C.coverage.gz $C $OUT
head -1 $OUT
