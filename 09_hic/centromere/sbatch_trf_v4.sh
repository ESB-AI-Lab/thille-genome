#!/bin/bash
#SBATCH -J cen.trf.v4
#SBATCH -o logs/cen.trf.v4.o%A.%a
#SBATCH -e logs/cen.trf.v4.e%A.%a
#SBATCH -c 4
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=8
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Tandem Repeats Finder per Thille chromosome (array task N = chrN) for the centromere satellite track.
# Parameters are TRF's recommended weights (Match 2, Mismatch 7, Delta 7, PM 80, PI 10, Minscore 50) with MaxPeriod 2000
# so satellite monomers of 100-2000 bp are reported. -c 4 only for memory headroom (TRF is single-threaded;
# 5 Mb of chr1 used 1.26 GB). Submit from assemblies/thille.hic: sbatch centromere/sbatch_trf.sh
set -euo pipefail
TRF=/gpool/bin/TRF/bin/trf
CHR=chr${SLURM_ARRAY_TASK_ID}
IN=compare/v4/thille_chrs/$CHR.fasta
OUT=centromere/v4/trf/$CHR.dat
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} chr=$CHR"
[[ -s $IN ]] || { echo "ERROR: missing $IN" >&2; exit 1; }
if [[ -f $OUT.done ]]; then echo "exists with .done, skipping"; exit 0; fi
/usr/bin/time -v $TRF "$IN" 2 7 7 80 10 50 2000 -h -ngs > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
date > "$OUT.done"
echo "arrays: $(grep -vc '^@' "$OUT")"
