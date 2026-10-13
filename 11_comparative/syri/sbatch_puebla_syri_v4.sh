#!/bin/bash
#SBATCH -J compare.puebla.syri.v4
#SBATCH -o logs/compare.puebla.syri.v4.o%A.%a
#SBATCH -e logs/compare.puebla.syri.v4.e%A.%a
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-1
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Puebla step 2, SyRI 1.7.0 with the same commands as sbatch_syri.sh (minimap2 -ax asm5 --eqx, -F B), query = the
# same-strand Puebla FASTA from sbatch_puebla_prep.sh:
#   task 0: Thille (reference) -> Puebla; pairwise plotsr
#   task 1: Hass UoQ same-strand FASTA (reference) -> Puebla, then the four-genome plotsr
#           West Indian -> Thille -> Hass -> Puebla (reuses compare/v4/syri/wi_ref and compare/v4/syri/hass)
#           and the three-genome plotsr West Indian -> Thille -> Puebla
# Submit from assemblies/thille.hic after sbatch_puebla_prep.sh: sbatch compare/sbatch_puebla_syri_v4.sh
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1
MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SAMTOOLS=/gpool/bin/samtools-1.24/samtools
SIF=/gpool/bin/singularityimages/graphgenomics.sif
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
HASS=compare/v4/syri/hass/hass.chr.samestrand.fasta
WI=compare/v4/syri/westindian/westindian.chr.samestrand.fasta
QRY=compare/v4/puebla/puebla.chr.samestrand.fasta
NAMES=(puebla hass_puebla)
REFS=("$THILLE" "$HASS")
NAME=${NAMES[$SLURM_ARRAY_TASK_ID]}
REF=${REFS[$SLURM_ARRAY_TASK_ID]}
OUT=compare/v4/syri/$NAME
BAM=$OUT/${NAME}.asm5.eqx.bam
T=${SLURM_CPUS_PER_TASK}
echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ref=$REF"
free -g
mkdir -p "$OUT" compare/v4/syri/fourway
if [[ ! -s $BAM.done ]]; then
  echo "[$(date)] minimap2 -ax asm5 --eqx start"
  $MINIMAP2 -ax asm5 --eqx -t "$T" "$REF" "$QRY" | $SAMTOOLS view -b -o "$BAM.tmp" -
  mv "$BAM.tmp" "$BAM"; date > "$BAM.done"
fi
cp "$QRY" "$OUT/puebla.chr.samestrand.fasta"
printf '#file\tname\ttags\n%s\tThille\tlw:1.5\n%s\tPuebla\tlw:1.5\n' "$PWD/$THILLE" "$PWD/$QRY" > "$OUT/genomes.txt"
seq 1 12 | sed "s/^/chr/" > "$OUT/chrord.txt"
singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
if [[ ! -s $OUT/${NAME}syri.out ]]; then
  cd $OUT
  echo \"[\$(date)] syri start\"
  python -c 'import sys; from syri.scripts.syri import main; sys.argv[0] = \"syri\"; main()' \
    -c $(basename "$BAM") -r $PWD/$REF -q puebla.chr.samestrand.fasta -k -F B --prefix $NAME --nc 12
  cd - > /dev/null
fi
if [[ $NAME == puebla ]]; then
  for ext in pdf png; do
    plotsr --sr $OUT/pueblasyri.out --genomes $OUT/genomes.txt --chrord $OUT/chrord.txt \
      -o $OUT/thille_vs_puebla_plotsr.\$ext -s 100000 -H 10 -W 7
  done
  # three-genome ribbon West Indian -> Thille -> Puebla (reuses the West Indian-reference SyRI)
  mkdir -p compare/v4/syri/wi_thille_puebla
  printf '#file\tname\ttags\n%s\tWest Indian\tlw:1.5\n%s\tThille\tlw:1.5\n%s\tPuebla\tlw:1.5\n' \
    $PWD/$WI $PWD/$THILLE $PWD/$QRY > compare/v4/syri/wi_thille_puebla/genomes.txt
  seq 1 12 | sed 's/^/chr/' > compare/v4/syri/wi_thille_puebla/chrord.txt
  for ext in pdf png; do
    plotsr --sr compare/v4/syri/wi_ref/wi_refsyri.out --sr $OUT/pueblasyri.out --genomes compare/v4/syri/wi_thille_puebla/genomes.txt \
      --chrord compare/v4/syri/wi_thille_puebla/chrord.txt -o compare/v4/syri/wi_thille_puebla/westindian_thille_puebla_plotsr.\$ext \
      -s 100000 -H 14 -W 8
  done
else
  printf '#file\tname\ttags\n%s\tWest Indian\tlw:1.5\n%s\tThille\tlw:1.5\n%s\tHass\tlw:1.5\n%s\tPuebla\tlw:1.5\n' \
    $PWD/$WI $PWD/$THILLE $PWD/$HASS $PWD/$QRY > compare/v4/syri/fourway/genomes.txt
  seq 1 12 | sed 's/^/chr/' > compare/v4/syri/fourway/chrord.txt
  for ext in pdf png; do
    plotsr --sr compare/v4/syri/wi_ref/wi_refsyri.out --sr compare/v4/syri/hass/hasssyri.out --sr $OUT/hass_pueblasyri.out \
      --genomes compare/v4/syri/fourway/genomes.txt --chrord compare/v4/syri/fourway/chrord.txt \
      -o compare/v4/syri/fourway/westindian_thille_hass_puebla_plotsr.\$ext -s 100000 -H 18 -W 8
  done
fi
"
echo "[$(date)] done"
