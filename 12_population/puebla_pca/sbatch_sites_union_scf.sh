#!/bin/bash
#SBATCH -J scf_sites
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -o puebla_pca/scf_sites.%j.log
# Plan v19: one bcftools targets file for genotyping on the Gwen scaffolds = union of the noOVL and noOVL.onlypsch
# published sites (CHROM POS REF,ALT; sites present in both carry the same alleles). Submit from assemblies/thille.hic.
set -euo pipefail
O=puebla_pca/scf
BC="singularity exec --no-home /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif"
[[ -e $O/sites_union.tsv.gz ]] && { echo "exists"; exit 1; }
zcat $O/panel_noOVL/sites_pub.tsv.gz $O/panel_noOVL.onlypsch/sites_pub.tsv.gz | sort -k1,1V -k2,2n -u -S 20G --parallel 4 \
  | $BC bgzip -c > $O/sites_union.tsv.gz
$BC tabix -s1 -b2 -e2 $O/sites_union.tsv.gz
echo "union sites: $(zcat $O/sites_union.tsv.gz | wc -l); positions with conflicting alleles: $(zcat $O/sites_union.tsv.gz | cut -f1,2 | uniq -d | wc -l)"
