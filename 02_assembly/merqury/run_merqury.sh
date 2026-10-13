#!/bin/bash
# Build a k-mer DB from HiFi reads with meryl, then run Merqury (QV, completeness, spectra-cn).
# Runs INSIDE the asmbenchmarkers container (see sbatch_merqury.sh).
#
# Usage: run_merqury.sh <reads.fastq.gz> <asm.fasta> <outdir> <out_prefix> [k] [threads] [meryl_mem_gb]
#   run_merqury.sh reads.fastq.gz asm.fasta merqury/k21 thille.hic 21 32 64
set -euo pipefail
export PYTHONUNBUFFERED=1

READS=$(readlink -f "$1")
ASM=$(readlink -f "$2")
OUTDIR=$3
PREFIX=$4
K=${5:-21}
THREADS=${6:-${SLURM_CPUS_PER_TASK:-8}}
MEM=${7:-64}

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
DB=$OUTDIR/$(basename "$READS" | sed -E 's/\.(fastq|fq)(\.gz)?$//').k${K}.meryl

echo "[$(date)] host=$(hostname) threads=$THREADS meryl_mem=${MEM}G k=$K"
echo "  reads=$READS"
echo "  asm=$ASM"
echo "  db=$DB"
meryl --version 2>&1 | head -1

# Sentinel lets a rerun skip the expensive count; a partial DB without it is not trusted.
if [[ -f $DB.done ]]; then
  echo "[$(date)] meryl DB exists with .done sentinel, skipping count"
else
  if [[ -e $DB ]]; then
    echo "ERROR: $DB exists without $DB.done (incomplete previous run?). Move it aside and rerun." >&2
    exit 1
  fi
  echo "[$(date)] meryl count start"
  meryl count k="$K" threads="$THREADS" memory="$MEM" "$READS" output "$DB"
  touch "$DB.done"
  echo "[$(date)] meryl count done"
fi

echo "[$(date)] merqury start"
cd "$OUTDIR"
merqury.sh "$DB" "$ASM" "$PREFIX"
echo "[$(date)] merqury done"
echo "--- QV ---";           cat "$PREFIX.qv"
echo "--- completeness ---"; cat "$PREFIX.completeness.stats"
