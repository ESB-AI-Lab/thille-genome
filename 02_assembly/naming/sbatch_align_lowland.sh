#!/bin/bash
#SBATCH -J name.lowland
#SBATCH -o name.lowland.o%j
#SBATCH -e name.lowland.e%j
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster
#SBATCH --export=ALL
# Curated Thille scaffolds (juicer post v2, chr8 joined) -> Lowland accession T2T (GCA_051132755.1, Pa01-Pa12), to name
# the chromosome-scale scaffolds chr1-chr12 (name_by_lowland.py). Approximate mapping (no -c) is enough for naming.
# --mem because the earlier whole-genome Lowland run needed ~100 GB.
# Submit from assemblies/thille.hic/curation/naming: sbatch sbatch_align_lowland.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
CUR=/epool/edwin/pamericana/avoflow/hifiasm/thille/hifi.hic.gcluster/yahs_thille.hifiv3.asm.hic.p_ctg.curation/out_JBAT.v2.FINAL.fa
WI=/epool/edwin/pamericana/avoflow/assemblies/thille.hic/westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna
PAF=curated_v2_vs_lowland.asm5.paf
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK mem=${SLURM_MEM_PER_NODE:-NA}MB job=$SLURM_JOB_ID"
echo "minimap2 $($MINIMAP2 --version)"
free -g
if [[ -s $PAF.done ]]; then echo "PAF exists with .done, skipping"; else
  echo "[$(date)] minimap2 -x asm5 start"
  $MINIMAP2 -x asm5 --secondary=no -t "$SLURM_CPUS_PER_TASK" "$WI" "$CUR" > "$PAF.tmp"
  mv "$PAF.tmp" "$PAF"; date > "$PAF.done"
fi
echo "[$(date)] done: $(wc -l < $PAF) alignments"
