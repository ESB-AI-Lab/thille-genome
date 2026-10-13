#!/bin/bash
#SBATCH -J inv.pileup
#SBATCH -o logs/inv.pileup.o%A.%a
#SBATCH -e logs/inv.pileup.e%A.%a
#SBATCH -c 2
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-33
#SBATCH --export=ALL
# Raw-read mpileups and breakpoint support for the 337 inversion windows (compare/v4/inversions/inv_pileup.py),
# 34 tasks of about 10 windows. One singularity exec; binds limited to the samtools build and the BAM directory
# (read-only); the working directory is bound by default.
# Submit from assemblies/thille.hic: sbatch compare/v4/inversions/sbatch_inv_pileup.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
free -g
singularity exec --no-home -B /gpool/bin/samtools-1.24:/gpool/bin/samtools-1.24:ro \
  -B /epool/edwin/pamericana/avoflow/mappingback:/epool/edwin/pamericana/avoflow/mappingback:ro \
  /gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif \
  python3 compare/v4/inversions/inv_pileup.py --task "$SLURM_ARRAY_TASK_ID" --ntasks 34
echo "[$(date)] done"
