#!/bin/bash
#SBATCH -J inv.pileup.pb
#SBATCH -o logs/inv.pileup.pb.o%A.%a
#SBATCH -e logs/inv.pileup.pb.e%A.%a
#SBATCH -c 2
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-33
#SBATCH --export=ALL
# Raw-read mpileups and breakpoint support for the 337 inversion windows (compare/v4/inversions/inv_pileup.py),
# 34 tasks of about 10 windows. One singularity exec; binds limited to the samtools build and the BAM directory
# (read-only); the working directory is bound by default.
# Puebla read sets (2026-10-06): Puebla HiFi on every Thille window and all four read sets on the Puebla windows;
# support rows go to support/task<N>_puebla.tsv. Adds a read-only bind of the Puebla reference directory.
# Submit from assemblies/thille.hic: sbatch compare/v4/inversions/sbatch_inv_pileup_puebla.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
free -g
singularity exec --no-home -B /gpool/bin/samtools-1.24:/gpool/bin/samtools-1.24:ro \
  -B /epool/edwin/pamericana/avoflow/mappingback:/epool/edwin/pamericana/avoflow/mappingback:ro \
  -B /gpool/edwin/avocados/thille_paper/genome/puebla/GCA_060308345.1/ncbi_dataset/data/GCA_060308345.1:/gpool/edwin/avocados/thille_paper/genome/puebla/GCA_060308345.1/ncbi_dataset/data/GCA_060308345.1:ro \
  /gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif \
  python3 compare/v4/inversions/inv_pileup.py --task "$SLURM_ARRAY_TASK_ID" --ntasks 34 --select "puebla:*" thille:puebla_hifi --tag _puebla
echo "[$(date)] done"
