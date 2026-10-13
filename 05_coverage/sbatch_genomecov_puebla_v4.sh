#!/bin/bash
#SBATCH --job-name=genomecov.puebla.thille_v4
#SBATCH --partition=debug
#SBATCH --cpus-per-task=12
#SBATCH --output=logs/genomecov_v4_puebla.hifi.%j.log
# Per-position coverage of Puebla HiFi on the Thille v4 chromosomes (BAMs from SLURM job 13173), same command as the
# other read sets: scripts/genomecov_track.sh (BEDTools v2.30.0 genomecov -ibam -d -split in aligners_v3.sif, 12
# chromosomes in parallel). Writes covfiles/puebla.hifi.thille_v4.chrN.coverage.gz; never overwrites.
# Submit from the repo root: sbatch scripts/sbatch_genomecov_puebla_v4.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node $(hostname)  cpus ${SLURM_CPUS_PER_TASK:-?}  start $(date)"
cd /gpool/edwin/avocados/thille_paper
scripts/genomecov_track.sh puebla.hifi.thille_v4 /epool/edwin/pamericana/avoflow/mappingback/puebla/minimap2/puebla.hifi.thille_v4 thille_chromosomes.tsv
echo "done $(date)"
