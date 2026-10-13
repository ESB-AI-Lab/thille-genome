#!/bin/bash
#SBATCH -J panel_pca
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 16
#SBATCH -o puebla_pca/panel_pca.%j.log
# Plan v19 Q2 prep and Q3. (1) Targets file of the Gwen-paper panel sites with their alleles, for
# bcftools call -C alleles (bcftools 1.24 manual: CHROM, POS, REF,ALT; bgzip + tabix -s1 -b2 -e2).
# (2) PLINK 2 on the same panel genotypes (noOG VCF, 30 samples): --freq counts and --pca allele-wts, as in the
# PLINK 2 PCA-projection docs; Q3 checks that the 30 sample PCs reproduce the published PLINK 1.9 eigenvec.
# Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/epool/edwin/pamericana/gwen/snpcalling/plink
VCF=$P/avocado.haplo.noOG.recodeall.maf.recode.vcf.gz
O=puebla_pca/panel
BC="singularity exec --no-home -B $P:$P:ro /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif"
PL="singularity exec --no-home -B $P:$P:ro /gpool/bin/singularityimages/plink2/plink2_2.0.0a.7.8.sif plink2"
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK job=$SLURM_JOB_ID"; free -g | head -2
mkdir -p $O; [[ -e $O/sites.tsv.gz ]] && { echo "exists: $O/sites.tsv.gz"; exit 1; }
$BC bash -c "bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $VCF | bgzip -@ 8 -c > $O/sites.tsv.gz && tabix -s1 -b2 -e2 $O/sites.tsv.gz"
echo "sites: $(zcat $O/sites.tsv.gz | wc -l)"
$PL --vcf $VCF --set-all-var-ids @:# --allow-extra-chr --threads $SLURM_CPUS_PER_TASK --memory 120000 \
  --make-pgen --out $O/panel
$PL --pfile $O/panel --allow-extra-chr --threads $SLURM_CPUS_PER_TASK --memory 120000 \
  --freq counts --pca allele-wts vcols=chrom,ref,alt --out $O/panel_pcs
echo "[$(date)] done"
