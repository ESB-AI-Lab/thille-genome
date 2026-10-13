#!/bin/bash
#SBATCH -J compare.syri3.v4
#SBATCH -o logs/compare.syri3.v4.o%j
#SBATCH -e logs/compare.syri3.v4.e%j
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Three-genome synteny plot, West Indian (top) -> Thille -> Hass (bottom), with plotsr.
# plotsr chains SyRI runs genome1->genome2, genome2->genome3 (reference = upper genome), all genomes with the same
# chromosome ids on the same strand. Thille->Hass is reused from sbatch_syri.sh (compare/v4/syri/hass). This job adds
# West Indian (reference, the chr1..chr12 same-strand FASTA built by sbatch_syri.sh) -> Thille (query), with the same
# commands as sbatch_syri.sh (minimap2 -ax asm5 --eqx, SyRI 1.7.0 -F B).
# Submit from assemblies/thille.hic: sbatch compare/sbatch_syri_threeway_v4.sh
# Outputs: compare/v4/syri/wi_ref/ (BAM, wi_refsyri.out) and compare/v4/syri/threeway/ (genomes.txt, plotsr PDF/PNG).
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SAMTOOLS=/gpool/bin/samtools-1.24/samtools
SIF=/gpool/bin/singularityimages/graphgenomics.sif
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
WI=compare/v4/syri/westindian/westindian.chr.samestrand.fasta
HASS=compare/v4/syri/hass/hass.chr.samestrand.fasta
OUT=compare/v4/syri/wi_ref
BAM=$OUT/westindian_vs_thille.asm5.eqx.bam
T=${SLURM_CPUS_PER_TASK}

echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
echo "minimap2 $($MINIMAP2 --version)  $($SAMTOOLS --version | head -1)"
free -g
mkdir -p "$OUT" compare/v4/syri/threeway

if [[ ! -s $BAM.done ]]; then
  echo "[$(date)] minimap2 -ax asm5 --eqx (ref West Indian, query Thille) start"
  $MINIMAP2 -ax asm5 --eqx -t "$T" "$WI" "$THILLE" | $SAMTOOLS view -b -o "$BAM.tmp" -
  mv "$BAM.tmp" "$BAM"; date > "$BAM.done"
fi

seq 1 12 | sed "s/^/chr/" > compare/v4/syri/threeway/chrord.txt
printf '#file\tname\ttags\n%s\tWest Indian\tlw:1.5\n%s\tThille\tlw:1.5\n%s\tHass\tlw:1.5\n' \
  "$PWD/$WI" "$PWD/$THILLE" "$PWD/$HASS" > compare/v4/syri/threeway/genomes.txt

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
if [[ ! -s $OUT/wi_refsyri.out ]]; then
  cd $OUT
  echo \"[\$(date)] syri start\"
  python -c 'import sys; from syri.scripts.syri import main; sys.argv[0] = \"syri\"; main()' \
    -c $(basename "$BAM") -r ../westindian/$(basename "$WI") -q ../../../../$THILLE -k -F B --prefix wi_ref --nc 12
  cd - > /dev/null
fi
echo \"[\$(date)] plotsr\"
for ext in pdf png; do
  plotsr --sr $OUT/wi_refsyri.out --sr compare/v4/syri/hass/hasssyri.out --genomes compare/v4/syri/threeway/genomes.txt \
    --chrord compare/v4/syri/threeway/chrord.txt -o compare/v4/syri/threeway/westindian_thille_hass_plotsr.\$ext -s 100000 -H 14 -W 8
done
"
echo "[$(date)] done"
