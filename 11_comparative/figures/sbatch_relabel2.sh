#!/bin/bash
#SBATCH -J relabel2.v4
#SBATCH -o logs/relabel2.v4.o%j
#SBATCH -e logs/relabel2.v4.e%j
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --export=ALL
# Edwin's display name (2026-10-05, refined): GCA_051132755.1 is "Unnamed Lowland" (no "accession"). Redraws the v4 figures that show
# it from the existing v4 SyRI/PAF results (no realignment) into compare/v4/relabel2/; earlier v4 figures stay as they are.
# --mem: plotsr loads three SyRI outputs and four genomes.
# Submit from assemblies/thille.hic: sbatch compare/v4/relabel2/sbatch_relabel2.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
R=compare/v4/relabel2
S=compare/v4/syri
GG=/gpool/bin/singularityimages/graphgenomics.sif
AB=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
free -g
mkdir -p $R/.mpl
singularity exec --no-home "$AB" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$R/.mpl
python3 compare/plot_syri_overlay.py $S/westindian/westindiansyri.out thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai \
  $S/westindian/westindian.chr.samestrand.fasta.fai $R/thille_vs_westindian_syri_overlay --y-name 'Unnamed Lowland'
for mode in orient raw; do
  extra=''; suffix=''
  if [[ \$mode == raw ]]; then extra='--no-orient'; suffix='_raw'; fi
  python3 compare/plot_genome_dotplot.py compare/v4/westindian/westindian_vs_thille.asm5.c.paf \
    compare/v4/genome_dotplots/thille_chroms.order.tsv compare/v4/westindian/westindian_chroms.order.tsv \
    $R/thille_vs_westindian_genome\$suffix --x-side target --x-name Thille --y-name 'Unnamed Lowland' \
    --title 'Thille vs Unnamed Lowland (XTBG_Pame_1.0), minimap2 asm5 -c' \$extra
done
B=compare/v4/wg_wi/plots/thille_wg_westindian_wg_rc/thille_wg_westindian_wg_rc_COMPARE.bedpe
python3 compare/plot_wg_compare.py \$B compare/v4/wg/thille_wg.boundaries.tsv compare/v4/wg_wi/westindian_wg_rc.boundaries.tsv \
  $R/thille_vs_westindian_wg_rc_boundaries --y-name 'Unnamed Lowland' \
  --note 'Unnamed Lowland Pa01,04,05,06,10,11 reverse-complemented (strand convention)'
"
singularity exec --no-home "$GG" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$R/.mpl
for ext in pdf png; do
  plotsr --sr $S/wi_ref/wi_refsyri.out --sr $S/hass/hasssyri.out --genomes $R/genomes_threeway.txt \
    --chrord $S/threeway/chrord.txt -o $R/westindian_thille_hass_plotsr.\$ext -s 100000 -H 14 -W 8
  plotsr --sr $S/wi_ref/wi_refsyri.out --sr $S/puebla/pueblasyri.out --genomes $R/genomes_wi_thille_puebla.txt \
    --chrord $S/wi_thille_puebla/chrord.txt -o $R/westindian_thille_puebla_plotsr.\$ext -s 100000 -H 14 -W 8
  plotsr --sr $S/wi_ref/wi_refsyri.out --sr $S/hass/hasssyri.out --sr $S/hass_puebla/hass_pueblasyri.out \
    --genomes $R/genomes_fourway.txt --chrord $S/fourway/chrord.txt \
    -o $R/westindian_thille_hass_puebla_plotsr.\$ext -s 100000 -H 18 -W 8
done
for ext in pdf png; do
  plotsr --sr $S/westindian/westindiansyri.out --genomes $R/genomes_pair.txt -o $R/thille_vs_westindian_plotsr.\$ext -s 100000 -H 10 -W 7
done
for r in R3:chr3:7598588-11696806 R10:chr10:10061230-15097080; do
  id=\${r%%:*}; reg=\${r#*:}
  for ext in pdf png; do
    plotsr --sr $S/wi_ref/wi_refsyri.out --sr $S/hass/hasssyri.out --sr $S/hass_puebla/hass_pueblasyri.out \
      --genomes $R/genomes_ribbons.txt --reg Thille:\$reg --rtr -o $R/\${id}_fourway_ribbon.\$ext -s 10000 -H 6 -W 9
  done
done
"
ls -la $R
