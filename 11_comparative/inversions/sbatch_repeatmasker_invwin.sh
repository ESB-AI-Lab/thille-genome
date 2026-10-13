#!/bin/bash
#SBATCH -J rm.invwin
#SBATCH -o logs/rm.invwin.o%A.%a
#SBATCH -e logs/rm.invwin.e%A.%a
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
# Copy (2026-10-05) for the inversion windows of Hass UoQ, Unnamed Lowland and Puebla (compare/v4/inversions/windows.tsv,
# original accession strand; sequence names <genome>__<accession>__<0-based window start>), 4 chunks.
# Submit from assemblies/thille.hic: sbatch compare/v4/inversions/te/sbatch_repeatmasker_invwin.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
SIF=/gpool/bin/singularityimages/repeatmodeler/repeatmodeler.sif
LIBDIR=/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic
WINS=(compare/v4/inversions/te/windows_part0.fa compare/v4/inversions/te/windows_part1.fa
      compare/v4/inversions/te/windows_part2.fa compare/v4/inversions/te/windows_part3.fa)
FA=${WINS[$SLURM_ARRAY_TASK_ID]}
OUT=compare/v4/inversions/te/$(basename "$FA" .fa)
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} window=$FA"
free -g
mkdir -p "$OUT"
singularity exec --no-home -B "$LIBDIR/consensi.fa.classified:/lib/consensi.fa.classified:ro" "$SIF" \
  RepeatMasker -pa "${SLURM_CPUS_PER_TASK}" -s -lib /lib/consensi.fa.classified -a -poly -ace -excln -html -gff \
  -dir "$OUT" "$FA"
echo "[$(date)] done"
