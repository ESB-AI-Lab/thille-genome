#!/bin/bash
#SBATCH -J cen.trfwin
#SBATCH -o logs/cen.trfwin.o%A.%a
#SBATCH -e logs/cen.trfwin.e%A.%a
#SBATCH -c 2
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --export=ALL
# TRF on one window of a chromosome (same settings as sbatch_trf.sh). Windows are WIN_BP long and start every
# STEP_BP from OFFSET_BP, so neighbouring windows overlap by WIN_BP - STEP_BP; array task i covers
# [OFFSET + i*STEP, OFFSET + i*STEP + WIN). Windows go to centromere/trf/<chr>_windows<TAG>/.
# Used for chr6, whose whole-chromosome TRF run did not finish in 10+ h. Coordinates are shifted back to the
# chromosome by merge_trf_windows.py.
# Submit from assemblies/thille.hic:
#   sbatch --array=0-8 centromere/sbatch_trf_windows.sh chr6 10000000 8000000            # whole chr6
#   sbatch --array=0-12 centromere/sbatch_trf_windows.sh chr6 2000000 1500000 24000000 _24-34Mb
set -euo pipefail
CHR=$1; WIN_BP=${2:-10000000}; STEP_BP=${3:-8000000}; OFFSET_BP=${4:-0}; TAG=${5:-}
i=${SLURM_ARRAY_TASK_ID}
D=centromere/trf/${CHR}_windows${TAG}
mkdir -p "$D"
S=$(( OFFSET_BP + i * STEP_BP )); E=$(( S + WIN_BP ))
W=$D/${CHR}_w${i}_${S}.fa
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${i} $CHR window $S-$E"
python3 - "$CHR" "$S" "$E" "$W" <<'PY'
import sys; sys.path.insert(0, "compare")
from fastaio import read_fasta, write_record
c, s, e, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
(h, seq), = read_fasta(f"compare/thille_chrs/{c}.fasta")
if s >= len(seq):
    sys.exit(f"window start {s} beyond {c} length {len(seq)}")
write_record(open(out, "w"), f"{c}_{s}", seq[s:e])
PY
/usr/bin/time -v /gpool/bin/TRF/bin/trf "$W" 2 7 7 80 10 50 2000 -h -ngs > "$W.dat.tmp"
mv "$W.dat.tmp" "$W.dat"
date > "$W.dat.done"
