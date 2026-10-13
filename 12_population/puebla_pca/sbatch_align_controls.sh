#!/bin/bash
#SBATCH -J ctrl2gwen
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 32
#SBATCH --array=0-1
#SBATCH -o puebla_pca/align_ctrl.%A_%a.log
# Plan v19 Q4 controls: Thille HiFi (task 0, -x map-hifi) and Gwen CLR (task 1, -x map-pb) on the Gwen primary
# contigs, aligned and genotyped exactly as Puebla, so that their projections can be compared with their own
# Illumina points in the panel PCA. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
MM=/gpool/bin/minimap2-2.30_x64-linux/minimap2; S=/gpool/bin/samtools-1.24/samtools
REF=/epool/edwin/pamericana/gwen/scaffolding/pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta
NAMES=(thille.hifi gwen.clr); PRE=(map-hifi map-pb)
FQS=(/epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz /epool/edwin/pamericana/gwen/rawdata/pacbio/gwen/gwen_pacb_subreads.fastq.gz)
i=$SLURM_ARRAY_TASK_ID; OUT=puebla_pca/${NAMES[$i]}.gwen_primary.bam
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK job=${SLURM_ARRAY_JOB_ID}_$i reads=${FQS[$i]}"; free -g | head -2
[[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
T=puebla_pca/tmp_${SLURM_ARRAY_JOB_ID}_$i; mkdir -p $T
$MM -ax ${PRE[$i]} --secondary=no -t $((SLURM_CPUS_PER_TASK - 4)) $REF ${FQS[$i]} | $S sort -@ 4 -m 4G -T $T/s -o $OUT.tmp -
mv $OUT.tmp $OUT; $S index -@ 8 $OUT; $S flagstat -@ 8 $OUT > $OUT.flagstat; rmdir $T
echo "[$(date)] done"
