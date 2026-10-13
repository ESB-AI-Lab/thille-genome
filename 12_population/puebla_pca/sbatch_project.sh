#!/bin/bash
#SBATCH -J pca_project
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 8
#SBATCH --mem=64G   # bcftools sort of ~34 M sites per sample and PLINK 2 scoring of the panel
#SBATCH -o puebla_pca/project.%j.log
# Plan v19 Q4 for one panel. Usage (from assemblies/thille.hic):
#   sbatch puebla_pca/sbatch_project.sh <panel_dir> <published.eigenvec> <geno_dir> <out_dir>
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node=$(hostname) job=$SLURM_JOB_ID args=$*"
[[ -e $4/projection.tsv ]] && { echo "exists: $4/projection.tsv"; exit 1; }
python3 puebla_pca/project.py "$1" "$2" "$3" "$4" --threads 8
echo "[$(date)] done"
