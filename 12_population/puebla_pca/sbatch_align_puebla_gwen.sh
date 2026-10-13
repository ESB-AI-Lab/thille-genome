#!/bin/bash
#SBATCH -J puebla2gwen
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 32
#SBATCH -o puebla_pca/align.%j.log
# Plan v19 Q1: Puebla HiFi (SRR36605778) on the Gwen primary contigs used by the Gwen-paper PCA, with the settings of
# the other Puebla alignments (minimap2 2.30 -ax map-hifi --secondary=no). Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
MM=/gpool/bin/minimap2-2.30_x64-linux/minimap2; S=/gpool/bin/samtools-1.24/samtools
REF=/epool/edwin/pamericana/gwen/scaffolding/pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta
FQ=/epool/edwin/pamericana/avoflow/rawdata/hifi/puebla/SRR36605778_subreads.fastq.gz
OUT=puebla_pca/puebla.hifi.gwen_primary.bam
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK job=$SLURM_JOB_ID"; free -g | head -2
[[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
mkdir -p puebla_pca/tmp_$SLURM_JOB_ID
$MM -ax map-hifi --secondary=no -t $((SLURM_CPUS_PER_TASK - 4)) $REF $FQ \
  | $S sort -@ 4 -m 4G -T puebla_pca/tmp_$SLURM_JOB_ID/s -o $OUT.tmp -
mv $OUT.tmp $OUT
$S index -@ 8 $OUT
$S flagstat -@ 8 $OUT > $OUT.flagstat
rmdir puebla_pca/tmp_$SLURM_JOB_ID
echo "[$(date)] done"
