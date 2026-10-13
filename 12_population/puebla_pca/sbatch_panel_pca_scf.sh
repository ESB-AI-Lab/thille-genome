#!/bin/bash
#SBATCH -J panel_pca_scf
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 16
#SBATCH --array=0-1
#SBATCH -o puebla_pca/panel_pca_scf.%A_%a.log
# Plan v19 Q3 on the published 2022 scaffold PCAs (task 0 noOVL, task 1 noOVL.onlypsch): PLINK 2 pgen from the panel
# VCF restricted to the published .map variants (biallelic), --freq counts, --pca allele-wts with --read-freq, and a
# bcftools targets file of the sites. --mem: ~26 M variants x 30 samples. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/epool/edwin/pamericana/gwen/scaffolding/scaffolds/plink/amcfilter/newanalysis
V=(noOVL noOVL.onlypsch); v=${V[$SLURM_ARRAY_TASK_ID]}
B=$P/gwen_avocado.scaffolds.haplo.amcfilter.$v.recodeall.maf.recode
O=puebla_pca/scf/panel_$v; mkdir -p $O
PL="singularity exec --no-home -B $P:$P:ro /gpool/bin/singularityimages/plink2/plink2_2.0.0a.7.8.sif plink2"
BC="singularity exec --no-home /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif"
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} panel=$v"
[[ -e $O/panel_pcs.eigenvec ]] && { echo "exists"; exit 1; }
# PLINK 2 reads scaffold names Chr1..Chr12 as chromosome codes 1..12 (ChrU stays), so @:# IDs are "10:1778" for
# "Chr10:1778"; convert the published .map IDs the same way (the first run matched only ChrU, 6.6 M variants)
cut -f2 $B.map | sed -E 's/^[Cc]hr([0-9]+):/\1:/' > $O/published_variant_ids.txt
echo "published ids: $(wc -l < $O/published_variant_ids.txt)"
$PL --vcf $B.vcf.gz --set-all-var-ids @:# --allow-extra-chr --threads 16 --memory 100000 \
  --extract $O/published_variant_ids.txt --max-alleles 2 --make-pgen --out $O/panel_pub
$PL --pfile $O/panel_pub --allow-extra-chr --threads 16 --memory 100000 --freq counts --out $O/panel_pub
$PL --pfile $O/panel_pub --allow-extra-chr --threads 16 --memory 100000 --read-freq $O/panel_pub.acount \
  --pca allele-wts vcols=chrom,ref,alt --out $O/panel_pcs
awk 'BEGIN{OFS="\t"} !/^#/{c=$1; if (c ~ /^[0-9]+$/) c="Chr"c; print c,$2,$4","$5}' $O/panel_pub.pvar | $BC bgzip -c > $O/sites_pub.tsv.gz
$BC tabix -s1 -b2 -e2 $O/sites_pub.tsv.gz
n=$(grep -vc "^#" $O/panel_pub.pvar); echo "[$(date)] done: $n variants kept of $(wc -l < $O/published_variant_ids.txt) published; $(zcat $O/sites_pub.tsv.gz | wc -l) sites"
[[ $n -eq $(wc -l < $O/published_variant_ids.txt) ]] || { echo "ERROR: not all published variants matched"; exit 1; }
