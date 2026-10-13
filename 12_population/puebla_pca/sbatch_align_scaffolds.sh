#!/bin/bash
#SBATCH -J lr2gwenscf
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 32
#SBATCH --array=0-2
#SBATCH -o puebla_pca/align_scf.%A_%a.log
# Plan v19 (update 2026-10-10): Puebla HiFi, Thille HiFi and Gwen CLR on the Gwen scaffolds used by the published
# Gwen-paper PCAs (gwen_scaffolds.fasta, 13 sequences), same minimap2 settings as before. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
MM=/gpool/bin/minimap2-2.30_x64-linux/minimap2; S=/gpool/bin/samtools-1.24/samtools
REF=/epool/edwin/pamericana/gwen/scaffolding/scaffolds/gwen_scaffolds.fasta
NAMES=(puebla.hifi thille.hifi gwen.clr); PRE=(map-hifi map-hifi map-pb)
FQS=(/epool/edwin/pamericana/avoflow/rawdata/hifi/puebla/SRR36605778_subreads.fastq.gz
     /epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz
     /epool/edwin/pamericana/gwen/rawdata/pacbio/gwen/gwen_pacb_subreads.fastq.gz)
i=$SLURM_ARRAY_TASK_ID; OUT=puebla_pca/scf/${NAMES[$i]}.gwen_scaffolds.bam
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_$i reads=${FQS[$i]}"; free -g | head -2
mkdir -p puebla_pca/scf; [[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
T=puebla_pca/scf/tmp_${SLURM_ARRAY_JOB_ID}_$i; mkdir -p $T
$MM -ax ${PRE[$i]} --secondary=no -t $((SLURM_CPUS_PER_TASK - 4)) $REF ${FQS[$i]} | $S sort -@ 4 -m 4G -T $T/s -o $OUT.tmp -
mv $OUT.tmp $OUT; $S index -@ 8 $OUT; $S flagstat -@ 8 $OUT > $OUT.flagstat; rmdir $T
echo "[$(date)] done"
