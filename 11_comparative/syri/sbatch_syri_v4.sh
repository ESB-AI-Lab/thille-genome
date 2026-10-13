#!/bin/bash
#SBATCH -J compare.syri.v4
#SBATCH -o logs/compare.syri.v4.o%A.%a
#SBATCH -e logs/compare.syri.v4.e%A.%a
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-1
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Structural rearrangements and large indels between Thille (reference) and a chromosome-level query with SyRI 1.7.0,
# following the SyRI docs: `minimap2 -ax asm5 --eqx ref qry`, BAM input (-F B), homologous chromosomes with identical
# ids on the same strand. Query chromosomes are renamed chr1..chr12 and the six stored on the opposite strand
# (majority minus-strand matches in the dot-plot flip tables: 1,4,5,6,10,11 in both Hass and West Indian) are
# reverse-complemented first. Task 0 = Hass UoQ, 1 = West Indian T2T. --mem: previous whole-genome minimap2 runs
# peaked at 69-73 GB.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_syri_v4.sh
# Outputs: compare/v4/syri/<name>/ (query FASTA, BAM, <name>syri.out, syri.vcf, syri.summary, plotsr PDF/PNG).
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SAMTOOLS=/gpool/bin/samtools-1.24/samtools
SIF=/gpool/bin/singularityimages/graphgenomics.sif
REF=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
RC="chr1 chr4 chr5 chr6 chr10 chr11"
NAMES=(hass westindian)
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
case $NAME in
  hass)       SRC=hass_UoQ_chr_only.fasta; IDS=(CM056809.1 CM056810.1 CM056811.1 CM056812.1 CM056813.1 CM056814.1
                                               CM056815.1 CM056816.1 CM056817.1 CM056818.1 CM056819.1 CM056820.1) ;;
  westindian) SRC=westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna
              IDS=(CP152278.1 CP152279.1 CP152280.1 CP152281.1 CP152282.1 CP152283.1
                   CP152284.1 CP152285.1 CP152286.1 CP152287.1 CP152288.1 CP152289.1) ;;
esac
T=${SLURM_CPUS_PER_TASK}
OUT=compare/v4/syri/$NAME
QRY=$OUT/${NAME}.chr.samestrand.fasta
BAM=$OUT/thille_vs_${NAME}.asm5.eqx.bam

echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} query=$NAME"
echo "minimap2 $($MINIMAP2 --version)  $($SAMTOOLS --version | head -1)"
free -g
mkdir -p "$OUT"

if [[ ! -s $QRY.done ]]; then
  [[ -s $SRC.fai ]] || $SAMTOOLS faidx "$SRC"
  : > "$QRY.tmp"
  for i in "${!IDS[@]}"; do
    chr=chr$((i + 1)); id=${IDS[$i]}
    if [[ " $RC " == *" $chr "* ]]; then flag=-i; else flag=; fi
    $SAMTOOLS faidx $flag --mark-strand no "$SRC" "$id" | sed "1s/^>.*/>$chr/" >> "$QRY.tmp"
    echo "$chr <- $id ${flag:+(reverse-complemented)}"
  done
  mv "$QRY.tmp" "$QRY"; $SAMTOOLS faidx "$QRY"; date > "$QRY.done"
fi

if [[ ! -s $BAM.done ]]; then
  echo "[$(date)] minimap2 -ax asm5 --eqx start"
  $MINIMAP2 -ax asm5 --eqx -t "$T" "$REF" "$QRY" | $SAMTOOLS view -b -o "$BAM.tmp" -
  mv "$BAM.tmp" "$BAM"; date > "$BAM.done"
fi

printf '#file\tname\ttags\n%s\tThille\tlw:1.5\n%s\t%s\tlw:1.5\n' "$PWD/$REF" "$PWD/$QRY" "$NAME" > "$OUT/genomes.txt"

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
cd $OUT
if [[ ! -s ${NAME}syri.out ]]; then
  echo \"[\$(date)] syri start\"
  python -c 'import sys; from syri.scripts.syri import main; sys.argv[0] = \"syri\"; main()' \
    -c $(basename "$BAM") -r ../../../../$REF -q $(basename "$QRY") -k -F B --prefix $NAME --nc 12
fi
echo \"[\$(date)] plotsr\"
plotsr --sr ${NAME}syri.out --genomes genomes.txt -o thille_vs_${NAME}_plotsr.pdf -s 100000 -H 10 -W 7
plotsr --sr ${NAME}syri.out --genomes genomes.txt -o thille_vs_${NAME}_plotsr.png -s 100000 -H 10 -W 7
"
echo "[$(date)] done"
