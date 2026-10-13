#!/usr/bin/env bash
# Coverage figures on the Thille v4 assembly (complete chr8; Edwin, 2026-10-04): the same figures as
# make_fst_coverage_figures.sh, make_pigment_coverage_figures.sh and the Figure 8 tandem-array figures, drawn with the
# v4 Thille tracks (<read set>.thille_v4), v4 Tiberius GTF and TE GFF, and the v4-relabelled BEDs and title maps
# (flowering_beds/v4/, summaries/v4/fst_alleles/). Hass and Lowland (WI) reference panels use their existing
# alignments; they are redrawn only because their titles carry Thille gene IDs. Thille-reference figures carry a fourth
# track, Puebla HiFi on Thille v4 (SLURM 13173; Edwin, 2026-10-05); Puebla is not aligned to Hass UoQ or Unnamed Lowland.
# No-text variant (Edwin, 2026-10-09): same figures with --no-text; titles and medians go to <stem>_notes.tsv for the legend md.
# Usage (from the repo root): bash scripts/make_coverage_figures_v4_notext.sh [outdir, default Results_v4_notext]
set -euo pipefail
O=${1:-Results_v4_notext}
B=flowering_beds/v4/igv; T=summaries/v4/fst_alleles
GTF=annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf
TE=annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff
TH="--tracks-b gwen.clr.thille_v4 thille.hifi.thille_v4 westindian.hifi.thille_v4 --labels-b Gwen_→_Thille Thille_→_Thille WI_→_Thille"
TH=$(echo "$TH" | sed 's/_→_/ → /g')
COMMON="--name-fields 2 --class-strip --png-dpi 300 --no-text"
G=$(grep -v '^#' $T/coverage_title_map_thille.tsv | cut -f1 | tr '\n' ' ')
GW=$(awk -F'\t' '!/^track/{split($4,a,"|"); print a[1]"|"a[2]}' $B/flowering_genes_wi_fst.bed | sort -u \
     | grep -Fx -f <(grep -v '^#' $T/coverage_title_map_wi.tsv | cut -f1) | tr '\n' ' ')
GP=$(grep -v '^#' $T/coverage_title_map_pigment_thille.tsv | cut -f1 | tr '\n' ' ')
thille() {  # $1 bed, $2 genes, $3 title map, $4 outdir
  python3 plot_gene_stacks.py -b "$1" --tracks-b gwen.clr.thille_v4 thille.hifi.thille_v4 westindian.hifi.thille_v4 puebla.hifi.thille_v4 \
    --labels-b "Gwen → Thille" "Thille → Thille" "Unnamed Lowland → Thille" "Puebla → Thille" --genes-gtf-b $GTF --te-gff-b $TE \
    --genes $2 $COMMON --plot-jobs 6 --title-map "$3" --outdir "$4"; }
thille $B/flowering_genes_thille.bed "$G" $T/coverage_title_map_thille.tsv $O/fst_coverage/thille_reference
thille $B/pigment_genes_thille.bed "$GP" $T/coverage_title_map_pigment_thille.tsv $O/pigment_coverage/thille_reference
python3 plot_gene_stacks.py -b $B/flowering_genes_hass.bed \
  --tracks-b gwen.clr.hass_uoq thille.hifi.hass_uoq westindian.hifi.hass_uoq --labels-b "Gwen → Hass" "Thille → Hass" "Unnamed Lowland → Hass" \
  --contig-names-b hass_cm_chromosomes.tsv --genes-gtf-b annotation/GCA_029852735.1_ASM2985273v1_genes.gtf \
  --genes $G $COMMON --plot-jobs 6 --title-map $T/coverage_title_map_hass.tsv --outdir $O/fst_coverage/hass_reference &
python3 plot_gene_stacks.py -b $B/flowering_genes_wi_fst.bed \
  --tracks-b gwen.clr.wi thille.hifi.wi westindian.hifi.wi --labels-b "Gwen → Unnamed Lowland" "Thille → Unnamed Lowland" "Unnamed Lowland → Unnamed Lowland" \
  --contig-names-b wi_cm_chromosomes.tsv --genes-gtf-b annotation/wi_GCA_051132755.1.genes.gtf \
  --genes $GW $COMMON --plot-jobs 6 --title-map $T/coverage_title_map_wi.tsv --outdir $O/fst_coverage/wi_reference &
python3 plot_gene_stacks.py -b $B/pigment_genes_hass.bed \
  --tracks-b gwen.clr.hass_uoq thille.hifi.hass_uoq westindian.hifi.hass_uoq --labels-b "Gwen → Hass" "Thille → Hass" "Unnamed Lowland → Hass" \
  --contig-names-b hass_cm_chromosomes.tsv --genes-gtf-b annotation/GCA_029852735.1_ASM2985273v1_genes.gtf \
  --genes $GP $COMMON --plot-jobs 5 --title-map $T/coverage_title_map_pigment_hass.tsv --outdir $O/pigment_coverage/hass_reference &
python3 plot_gene_stacks.py -b $B/pigment_genes_wi.bed \
  --tracks-b gwen.clr.wi thille.hifi.wi westindian.hifi.wi --labels-b "Gwen → Unnamed Lowland" "Thille → Unnamed Lowland" "Unnamed Lowland → Unnamed Lowland" \
  --contig-names-b wi_cm_chromosomes.tsv --genes-gtf-b annotation/wi_GCA_051132755.1.genes.gtf \
  --genes $GP $COMMON --plot-jobs 5 --title-map $T/coverage_title_map_pigment_wi.tsv --outdir $O/pigment_coverage/wi_reference &
wait
# Figure 8: selected tandem arrays (array copies track, other genes, TEs)
python3 plot_gene_stacks.py --tandem-figures --tandem-b $B/tandem_arrays_selected_thille.bed \
  --tracks-b gwen.clr.thille_v4 thille.hifi.thille_v4 westindian.hifi.thille_v4 puebla.hifi.thille_v4 \
  --labels-b "Gwen → Thille" "Thille → Thille" "Unnamed Lowland → Thille" "Puebla → Thille" \
  --genes-gtf-b $GTF --te-gff-b $TE --member-models $B/tandem_arrays_selected_members.tsv \
  --class-strip --png-dpi 300 --no-text --plot-jobs 4 --outdir $O/tandem_arrays
echo FIGURES_DONE
