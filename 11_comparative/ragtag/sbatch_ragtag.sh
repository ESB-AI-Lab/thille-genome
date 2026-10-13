#!/bin/bash
#SBATCH -J ragtag.v4
#SBATCH -o logs/ragtag.v4.o%A.%a
#SBATCH -e logs/ragtag.v4.e%A.%a
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster
#SBATCH --array=0-2
#SBATCH --export=ALL
# RagTag v2.1.1 reference-guided scaffolding on the Thille v4 chromosomes (PLAN.md v16). Contigs are only ordered and
# oriented (no ragtag correct: Gwen contigs are never broken, per Edwin). Task 0 = Gwen primary contigs, 1 = Gwen alt
# (HapSolo secondary) contigs, 2 = Hass UoQ split at its N gaps (split_at_gaps.py). minimap2 2.30 -x asm5 (RagTag
# default), -u adds the _RagTag suffix for AGP compatibility, 100 bp N per join (gap sizes not inferred).
# --mem: whole-genome minimap2 asm5 runs here peaked at 55-73 GB.
# Submit from assemblies/ragtag_thille_v4: sbatch sbatch_ragtag.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
SIF=/gpool/bin/singularityimages/ragtag/ragtag_2.1.1.sif
REF=/epool/edwin/pamericana/avoflow/assemblies/thille.hic/thille.hifiv4.asm.hic.scaffold.chr.only.fasta
G=/epool/edwin/pamericana/gwen/scaffolding
P=pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004
NAMES=(gwen_primary gwen_alt hass_uoq_broken)
# Gwen contigs are read through links in inputs/ so RagTag can write its .fai index (the Gwen folder is bound read-only)
QRYS=(inputs/gwen_primary_contigs.fasta inputs/gwen_alt_contigs.fasta hass_uoq_broken/hass_uoq_contigs.fasta)
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}; QRY=${QRYS[$SLURM_ARRAY_TASK_ID]}
OUT=$NAME/ragtag_scaffold
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} query=$QRY"
free -g
[[ -e $OUT/ragtag.scaffold.fasta ]] && { echo "ERROR: $OUT/ragtag.scaffold.fasta exists" >&2; exit 1; }
mkdir -p "$NAME"
singularity exec --no-home -B "$REF:$REF:ro" -B "$REF.fai:$REF.fai:ro" -B "$G:$G:ro" "$SIF" bash -c "
set -euo pipefail
ragtag.py --version; minimap2 --version
ragtag.py scaffold -u -t $SLURM_CPUS_PER_TASK --mm2-params '-x asm5 -t $SLURM_CPUS_PER_TASK' -o $OUT $REF $QRY
"
echo "[$(date)] done"; grep -c '>' "$OUT/ragtag.scaffold.fasta"
