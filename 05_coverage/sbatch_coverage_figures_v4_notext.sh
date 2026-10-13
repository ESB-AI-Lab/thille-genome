#!/bin/bash
#SBATCH --job-name=covfig.notext
#SBATCH --partition=gcluster,gcbigmem
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/covfig_notext.%j.log
# Usage (repo root): sbatch scripts/sbatch_coverage_figures_v4_notext.sh
export PYTHONUNBUFFERED=1
echo "node $(hostname) cpus $SLURM_CPUS_PER_TASK start $(date)"
bash scripts/make_coverage_figures_v4_notext.sh
echo "end $(date)"
