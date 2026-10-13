#!/bin/bash
# Plan v19 Q2 prep: published biallelic PCA sites as a bcftools targets file (CHROM POS REF,ALT; bgzip + tabix, per the
# bcftools manual) and the contig list with lengths for chunking. Run from assemblies/thille.hic.
set -euo pipefail
O=puebla_pca/panel
BC="singularity exec --no-home /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif"
[[ -e $O/sites_pub.tsv.gz ]] && { echo "exists"; exit 1; }
awk 'BEGIN{OFS="\t"} !/^#/{print $1,$2,$4","$5}' $O/panel_pub.pvar | $BC bgzip -c > $O/sites_pub.tsv.gz
$BC tabix -s1 -b2 -e2 $O/sites_pub.tsv.gz
# the Gwen contigs are linked into puebla_pca/ref/ and indexed there (samtools faidx); nothing is written next to the Gwen data
mkdir -p puebla_pca/ref
[[ -e puebla_pca/ref/gwen_primary.fasta ]] || ln -s /epool/edwin/pamericana/gwen/scaffolding/pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta puebla_pca/ref/gwen_primary.fasta
[[ -e puebla_pca/ref/gwen_primary.fasta.fai ]] || /gpool/bin/samtools-1.24/samtools faidx puebla_pca/ref/gwen_primary.fasta
cut -f1,2 puebla_pca/ref/gwen_primary.fasta.fai > $O/contigs_pub.tsv
echo "sites: $(zcat $O/sites_pub.tsv.gz | wc -l); contigs: $(wc -l < $O/contigs_pub.tsv)"
