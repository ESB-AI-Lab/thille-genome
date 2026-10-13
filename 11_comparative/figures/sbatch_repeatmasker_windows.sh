#!/bin/bash
#SBATCH -J rm.windows
#SBATCH -o logs/rm.windows.o%A.%a
#SBATCH -e logs/rm.windows.e%A.%a
#SBATCH -c 16
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-3
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# TE annotation of the Hass UoQ and West Indian plot windows around the chr3 (R3) and chr10 (R10) SV loci, with the same
# RepeatMasker command and library used for the Thille genome (annotation/repeatmodeler/thille.hic/sbatch_repeatmasker.sh:
# -s -lib consensi.fa.classified -a -poly -ace -excln -html -gff), so the three TE tracks are comparable.
# Windows are the same-strand FASTAs in compare/regions/trf/win_<region>_<genome>_<start>.fa.
# Submit from assemblies/thille.hic: sbatch compare/regions/sbatch_repeatmasker_windows.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
SIF=/gpool/bin/singularityimages/repeatmodeler/repeatmodeler.sif
LIBDIR=/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic
WINS=(compare/regions/trf/win_R3_hass_9041804.fa compare/regions/trf/win_R3_westindian_7601466.fa
      compare/regions/trf/win_R10_hass_9734504.fa compare/regions/trf/win_R10_westindian_9639957.fa)
FA=${WINS[$SLURM_ARRAY_TASK_ID]}
OUT=compare/regions/repeatmasker/$(basename "$FA" .fa)
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} window=$FA"
free -g
mkdir -p "$OUT"
singularity exec --no-home -B "$LIBDIR/consensi.fa.classified:/lib/consensi.fa.classified:ro" "$SIF" \
  RepeatMasker -pa "${SLURM_CPUS_PER_TASK}" -s -lib /lib/consensi.fa.classified -a -poly -ace -excln -html -gff \
  -dir "$OUT" "$FA"
echo "[$(date)] done"
