#!/bin/bash
#SBATCH -J pca_geno_scf
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 2
#SBATCH --mem=48G   # bcftools -T holds the target sites in memory
#SBATCH --array=0-38
#SBATCH -o puebla_pca/scf/geno.%A_%a.log
# Plan v19 Q2 on the Gwen scaffolds: one sample (task / 13: puebla.hifi, thille.hifi, gwen.clr) on one of the 13
# scaffold sequences (task % 13), at the union of the published scaffold PCA sites, same bcftools 1.24 commands as the
# contig genotyping (sbatch_genotype.sh). Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1
O=puebla_pca/scf; G=$O/geno; mkdir -p $G
REF=/epool/edwin/pamericana/gwen/scaffolding/scaffolds/gwen_scaffolds.fasta
SAMPLES=(puebla.hifi thille.hifi gwen.clr); s=${SAMPLES[$((SLURM_ARRAY_TASK_ID / 13))]}
seq=$(sed -n "$((SLURM_ARRAY_TASK_ID % 13 + 1))p" $REF.fai | cut -f1)
BC="singularity exec --no-home -B /epool/edwin/pamericana/gwen/scaffolding:/epool/edwin/pamericana/gwen/scaffolding:ro /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif bcftools"
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} sample=$s seq=$seq"
OUT=$G/$s.$seq.vcf.gz; [[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
$BC mpileup -f $REF -r $seq -T $O/sites_union.tsv.gz -q 20 -a FORMAT/AD,FORMAT/DP -d 500 -Ou $O/$s.gwen_scaffolds.bam \
  | $BC call -m -C alleles -T $O/sites_union.tsv.gz -i -Ou \
  | $BC filter -S . -e 'FMT/DP<5' -Oz -o $OUT.tmp
mv $OUT.tmp $OUT
echo "[$(date)] done: $($BC view -H $OUT | wc -l) sites"
