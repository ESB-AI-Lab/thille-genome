#!/usr/bin/env bash
# Coverage plots (no figure text) for EVERY flowering gene in the v4 BEDs (Edwin, 2026-10-12: the publication plots
# folder holds all flowering genes, not only those in the paper), one folder per reference genome.
# Same tracks and options as make_coverage_figures_v4_notext.sh; no --genes filter.
# Usage (repo root, SLURM): bash scripts/make_coverage_figures_v4_all.sh [outdir, default Results_v4_notext_all]
set -euo pipefail
O=${1:-Results_v4_notext_all}
B=flowering_beds/v4/igv
GTF=annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf
TE=annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff
COMMON="--name-fields 2 --class-strip --png-dpi 300 --no-text --plot-jobs 6"
python3 plot_gene_stacks.py -b $B/flowering_genes_thille.bed \
  --tracks-b gwen.clr.thille_v4 thille.hifi.thille_v4 westindian.hifi.thille_v4 puebla.hifi.thille_v4 \
  --labels-b "Gwen → Thille" "Thille → Thille" "Unnamed Lowland → Thille" "Puebla → Thille" \
  --genes-gtf-b $GTF --te-gff-b $TE $COMMON --outdir $O/thille_reference
python3 plot_gene_stacks.py -b $B/flowering_genes_hass.bed \
  --tracks-b gwen.clr.hass_uoq thille.hifi.hass_uoq westindian.hifi.hass_uoq --labels-b "Gwen → Hass" "Thille → Hass" "Unnamed Lowland → Hass" \
  --contig-names-b hass_cm_chromosomes.tsv --genes-gtf-b annotation/GCA_029852735.1_ASM2985273v1_genes.gtf \
  $COMMON --outdir $O/hass_reference
python3 plot_gene_stacks.py -b $B/flowering_genes_wi_all.bed \
  --tracks-b gwen.clr.wi thille.hifi.wi westindian.hifi.wi --labels-b "Gwen → Unnamed Lowland" "Thille → Unnamed Lowland" "Unnamed Lowland → Unnamed Lowland" \
  --contig-names-b wi_cm_chromosomes.tsv --genes-gtf-b annotation/wi_GCA_051132755.1.genes.gtf \
  $COMMON --outdir $O/unnamed_lowland_reference
echo ALL_DONE
