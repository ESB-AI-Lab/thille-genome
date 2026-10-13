#!/bin/bash
#SBATCH -J notext_figs2
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o compare/v4/notext/notext2.%j.log
# Plan v18 F2, second part: the panels job 13303 did not reach (it stopped at the Hi-C map, whose strict-filter
# weights live in all.variants.mcool, not all.mcool) plus Figure S2 redrawn after the axis-label fix. Same commands as
# sbatch_notext.sh. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/compare/v4/notext/.mpl
O=compare/v4/notext; G=compare/v4/genome_dotplots
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
CP=/epool/edwin/pamericana/avoflow/annotation/annotation_plots
F=/gpool/edwin/avocados/thille_paper/flowering_beds/v4/igv
L=/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/consensi.fa.classified
H=/epool/edwin/pamericana/avoflow/z_GCA_029852735.1/ncbi_dataset/data/GCA_029852735.1
echo "node=$(hostname) job=$SLURM_JOB_ID"
for f in $O/S2_Gwen.png $O/all_genome_500kb_V2_strict-filters_notext.png compare/v4/inversions/figures_v2_notext; do [[ -e $f ]] && { echo "exists: $f"; exit 1; }; done
run() { singularity exec --no-home -B $CP:$CP:ro -B $F:$F:ro -B $L:$L:ro -B $H:$H:ro "$SIF" "$@"; }
run python3 compare/plot_genome_dotplot.py compare/v4/gwen_c/gwen_vs_thille.asm5.c.paf $G/thille_chroms.order.tsv \
  $G/gwen_contigs_by_thille_asm5c.order.tsv $O/S2_Gwen --x-side target --x-name "Thille chr1 to chr12" \
  --y-name "Gwen contigs" --no-text
run python3 hic/plot_contacts.py hic/v4/all/all.variants.mcool all $O hic/v4/all/all.pairstats.tsv --genome-only \
  --weight-name weight_strict --suffix _V2_strict-filters_notext --no-text
run python3 karyotype/render_karyotype.py "$CP" hic/v4/sizes/chr12.sizes \
  thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff \
  $O/karyotype --calls satellitefirst:centromere/v4/thille_centromeres.satellitefirst --styles split-idiogram \
  --colors lavender:#b39ddb --dots fst_check/karyotype_dots.tsv --no-text
run python3 compare/v4/inversions/plot_inversions.py --outdir figures_v2_notext --no-text
run python3 compare/v4/inversions/summarize_support.py --no-text --fig-out compare/v4/inversions/figures_v2_notext/overview_split_fraction
echo "[$(date)] done"
