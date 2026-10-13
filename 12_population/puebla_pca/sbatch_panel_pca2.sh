#!/bin/bash
#SBATCH -J panel_pca2
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 16
#SBATCH -o puebla_pca/panel_pca2.%j.log
# Plan v19 Q3, second attempt. Job 13327 made panel.pgen and the sites file, but PLINK 2 refuses --pca with fewer
# than 50 samples unless allele frequencies come from --read-freq. Here: (1) restrict to the 34,472,447 variants of
# the published PLINK 1.9 PCA (IDs from its .map, biallelic), (2) --freq counts, (3) --pca allele-wts with
# --read-freq of those counts. --mem: the panel has 34 M variants x 30 samples. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
P=/epool/edwin/pamericana/gwen/snpcalling/plink
O=puebla_pca/panel
PL="singularity exec --no-home -B $P:$P:ro /gpool/bin/singularityimages/plink2/plink2_2.0.0a.7.8.sif plink2"
echo "node=$(hostname) job=$SLURM_JOB_ID"
[[ -e $O/panel_pcs.eigenvec ]] && { echo "exists"; exit 1; }
cut -f2 $P/avocado.haplo.noOG.recodeall.maf.recode.map > $O/published_variant_ids.txt
echo "published ids: $(wc -l < $O/published_variant_ids.txt), unique: $(sort -u $O/published_variant_ids.txt | wc -l)"
$PL --pfile $O/panel --allow-extra-chr --threads 16 --memory 120000 --extract $O/published_variant_ids.txt \
  --max-alleles 2 --make-pgen --out $O/panel_pub
$PL --pfile $O/panel_pub --allow-extra-chr --threads 16 --memory 120000 --freq counts --out $O/panel_pub
$PL --pfile $O/panel_pub --allow-extra-chr --threads 16 --memory 120000 --read-freq $O/panel_pub.acount \
  --pca allele-wts vcols=chrom,ref,alt --out $O/panel_pcs
echo "[$(date)] done"
