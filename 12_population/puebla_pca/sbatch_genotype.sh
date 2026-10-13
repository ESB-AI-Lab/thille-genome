#!/bin/bash
#SBATCH -J pca_geno
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 2
#SBATCH --mem=48G   # bcftools -T holds the 34 M target sites in memory (job 13331 was OOM-killed at the default)
#SBATCH --array=0-23
#SBATCH -o puebla_pca/geno.%A_%a.log
# Plan v19 Q2. Genotype one long-read sample (task / 8: puebla.hifi, thille.hifi, gwen.clr) on one of 8 contig chunks
# (task % 8) at the 34,472,447 published PCA sites, with the panel alleles: bcftools 1.24 mpileup -T sites (MAPQ >= 20,
# AD/DP annotated) | call -m -C alleles -T sites -i (sites with no coverage written as missing), then genotypes with
# DP < 5 set to missing. Submit from assemblies/thille.hic after the alignments.
set -euo pipefail
export PYTHONUNBUFFERED=1
O=puebla_pca; G=$O/geno; mkdir -p $G
# reference: link to the Gwen primary contigs, indexed here (no files written next to the Gwen data)
REF=$O/ref/gwen_primary.fasta
RD=/epool/edwin/pamericana/gwen/scaffolding
SAMPLES=(puebla.hifi thille.hifi gwen.clr); s=${SAMPLES[$((SLURM_ARRAY_TASK_ID / 8))]}; c=$((SLURM_ARRAY_TASK_ID % 8))
BC="singularity exec --no-home -B $RD:$RD:ro /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif bcftools"
echo "node=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} sample=$s chunk=$c"
# sites file restricted to the published biallelic variants, and the chunk's contigs (round robin over contig order)
if [[ ! -s $O/panel/sites_pub.tsv.gz.tbi ]]; then echo "missing $O/panel/sites_pub.tsv.gz.tbi"; exit 1; fi
awk -v c=$c 'NR % 8 == c {print $1"\t0\t"$2}' $O/panel/contigs_pub.tsv > $G/$s.chunk$c.bed
OUT=$G/$s.chunk$c.vcf.gz; [[ -e $OUT ]] && { echo "exists: $OUT"; exit 1; }
$BC mpileup -f $REF -R $G/$s.chunk$c.bed -T $O/panel/sites_pub.tsv.gz -q 20 -a FORMAT/AD,FORMAT/DP -d 500 -Ou $O/$s.gwen_primary.bam \
  | $BC call -m -C alleles -T $O/panel/sites_pub.tsv.gz -i -Ou \
  | $BC filter -S . -e 'FMT/DP<5' -Oz -o $OUT.tmp
mv $OUT.tmp $OUT
echo "[$(date)] done: $($BC view -H $OUT | wc -l) sites"
