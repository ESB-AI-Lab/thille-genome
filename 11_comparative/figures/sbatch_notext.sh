#!/bin/bash
#SBATCH -J notext_figs
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o compare/v4/notext/notext.%j.log
# Plan v18 F2 (Edwin, 2026-10-09): redraw this session's manuscript panels without titles, subtitles or notes
# (--no-text); that text moves to manuscript/figure_legends/Figure_<id>.md. Same inputs and options as the current
# panels; outputs only in compare/v4/notext/ and compare/v4/inversions/figures_v2_notext/ (nothing overwritten).
# --mem: the Hi-C genome map and the karyotype load the whole genome. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/compare/v4/notext/.mpl
mkdir -p $MPLCONFIGDIR
O=compare/v4/notext; S=compare/v4/syri; G=compare/v4/genome_dotplots
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
CP=/epool/edwin/pamericana/avoflow/annotation/annotation_plots
F=/gpool/edwin/avocados/thille_paper/flowering_beds/v4/igv
L=/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/consensi.fa.classified
H=/epool/edwin/pamericana/avoflow/z_GCA_029852735.1/ncbi_dataset/data/GCA_029852735.1
A=/epool/edwin/pamericana/avoflow/annotation/tiberius/thille.hic
echo "node=$(hostname) job=$SLURM_JOB_ID"
for f in ${CHECK:-$O/Fig5A_Hass.png $O/S15A_Hass.png}; do [[ -e $f ]] && { echo "exists: $f"; exit 1; }; done
run() { singularity exec --no-home -B $CP:$CP:ro -B $F:$F:ro -B $L:$L:ro -B $H:$H:ro "$SIF" "$@"; }
# Figure 5A: SyRI dot plots
run python3 compare/plot_syri_overlay.py $S/hass/hasssyri.out thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai \
  $S/hass/hass.chr.samestrand.fasta.fai $O/Fig5A_Hass --y-name 'Hass UoQ' --no-text
run python3 compare/plot_syri_overlay.py $S/westindian/westindiansyri.out thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai \
  $S/westindian/westindian.chr.samestrand.fasta.fai $O/Fig5A_Unnamed_Lowland --y-name 'Unnamed Lowland' --no-text
run python3 compare/plot_syri_overlay.py $S/puebla/pueblasyri.out thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai \
  compare/v4/puebla/puebla.chr.samestrand.fasta.fai $O/Fig5A_Puebla --y-name Puebla --no-text
# Figure S2: Thille vs Gwen contigs
run python3 compare/plot_genome_dotplot.py compare/v4/gwen_c/gwen_vs_thille.asm5.c.paf $G/thille_chroms.order.tsv \
  $G/gwen_contigs_by_thille_asm5c.order.tsv $O/S2_Gwen --x-side target --x-name "Thille chr1 to chr12" \
  --y-name "Gwen contigs" --no-text
# Figure S15: ModDotPlot
run python3 compare/plot_wg_compare.py compare/v4/wg_rc/plots/thille_wg_hass_wg_rc/thille_wg_hass_wg_rc_COMPARE.bedpe \
  compare/v4/wg/thille_wg.boundaries.tsv compare/v4/wg_rc/hass_wg_rc.boundaries.tsv $O/S15A_Hass --y-name 'Hass UoQ' --no-text
run python3 compare/plot_wg_compare.py compare/v4/wg_wi/plots/thille_wg_westindian_wg_rc/thille_wg_westindian_wg_rc_COMPARE.bedpe \
  compare/v4/wg/thille_wg.boundaries.tsv compare/v4/wg_wi/westindian_wg_rc.boundaries.tsv $O/S15B_Unnamed_Lowland --y-name 'Unnamed Lowland' --no-text
run python3 compare/plot_wg_compare.py compare/v4/wg_puebla/plots/thille_wg_puebla_wg_rc/thille_wg_puebla_wg_rc_COMPARE.bedpe \
  compare/v4/wg/thille_wg.boundaries.tsv compare/v4/wg_puebla/puebla_wg_rc.boundaries.rclabels.tsv $O/S15C_Puebla --y-name Puebla --no-text
# Figure S9: region tracks (host python, as for the current panels)
for r in R3 R10; do python3 compare/plot_region_tracks_v4.py $r $O/S9_${r}_tracks --no-text 2>&1 | grep -v Matplotlib; done
# Figure 3B: BUSCO
python3 $A/busco_plot.py -i $A/busco/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.proteins_embryophyta_odb12.2/short_summary.specific.embryophyta_odb12.2.thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.proteins_embryophyta_odb12.2.txt \
  $A/busco/thille.hifiv4.asm.hic.scaffold.masked.chr.only_embryophyta_odb12.2/short_summary.specific.embryophyta_odb12.2.thille.hifiv4.asm.hic.scaffold.masked.chr.only_embryophyta_odb12.2.txt \
  --labels "Thille Proteins" "Thille HiC Scaffolds" --no-title --annotation none -o $O/Fig3B_busco
# Figure 4: Hi-C genome map (V2 strict filters), no text
run python3 hic/plot_contacts.py hic/v4/all/all.variants.mcool all $O hic/v4/all/all.pairstats.tsv --genome-only \
  --weight-name weight_strict --suffix _V2_strict-filters_notext --no-text
# Figure 5C: karyotype
run python3 karyotype/render_karyotype.py "$CP" hic/v4/sizes/chr12.sizes \
  thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff \
  $O/karyotype --calls satellitefirst:centromere/v4/thille_centromeres.satellitefirst --styles split-idiogram \
  --colors lavender:#b39ddb --dots fst_check/karyotype_dots.tsv --no-text
# Figure S12: inversion breakpoints and overview
run python3 compare/v4/inversions/plot_inversions.py --outdir figures_v2_notext --no-text
run python3 compare/v4/inversions/summarize_support.py --no-text --fig-out compare/v4/inversions/figures_v2_notext/overview_split_fraction
echo "[$(date)] done"
