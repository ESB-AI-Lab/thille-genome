#!/bin/bash
#SBATCH -J merqury
#SBATCH -o logs/merqury.o%j
#SBATCH -e logs/merqury.e%j
#SBATCH -c 32
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Submit from assemblies/thille.hic (auto-bound as the container cwd):
#   sbatch -J merqury.thille merqury/sbatch_merqury.sh <reads.fastq.gz> <asm.fasta in cwd> <name> [k]
#   sbatch -J merqury.thille merqury/sbatch_merqury.sh /epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz \
#          thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta thille_hifi 21
# Outputs: merqury/<name>_k<k>/ ; only the reads directory is bound (read-only).
set -euo pipefail
export PYTHONUNBUFFERED=1

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
READS=$(readlink -f "$1")
ASM=$2
NAME=$3
K=${4:-21}
READDIR=$(dirname "$READS")

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_JOB_ID}"
echo "reads=$READS asm=$ASM name=$NAME k=$K bind=$READDIR:ro"
free -g

singularity exec --no-home -B "$READDIR:$READDIR:ro" "$SIF" \
  bash merqury/run_merqury.sh "$READS" "$ASM" "merqury/${NAME}_k${K}" "$NAME" "$K" "${SLURM_CPUS_PER_TASK}" 64
