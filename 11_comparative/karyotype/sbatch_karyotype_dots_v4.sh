#!/bin/bash
#SBATCH -J karyo.dots.v4
#SBATCH -o logs/karyo.dots.v4.o%j
#SBATCH -e logs/karyo.dots.v4.e%j
#SBATCH -c 4
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --export=ALL
# v4 copy (2026-10-04) of sbatch_karyotype_dots.sh and sbatch_karyotype_dots_hicfirst.sh: Figure 5C karyotype
# (split-idiogram, lavender centromere, flowering-gene dots) for both centromere call sets on the v4 genome.
# The v4 TE GFF already uses the assembly's chromosome names, so no renaming step is needed. The dots
# (fst_check/karyotype_dots.tsv) are coordinates; none is on chr8 and the other chromosomes are identical in v4.
# Writes karyotype/v4/figures/thille_karyotype_{satellitefirst,hicfirst}_split-idiogram_lavender_flowering.{pdf,png}.
# Submit from assemblies/thille.hic: sbatch karyotype/sbatch_karyotype_dots_v4.sh
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/hic/.mplconfig
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
CP=/epool/edwin/pamericana/avoflow/annotation/annotation_plots
mkdir -p karyotype/v4/figures
echo "node=$(hostname) job=${SLURM_JOB_ID}"
for m in satellitefirst hicfirst; do
  singularity exec --no-home -B $CP:$CP:ro "$SIF" python3 karyotype/render_karyotype.py "$CP" hic/v4/sizes/chr12.sizes \
    thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf \
    thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff \
    karyotype/v4/figures --calls $m:centromere/v4/thille_centromeres.$m --styles split-idiogram \
    --colors lavender:#b39ddb --dots fst_check/karyotype_dots.tsv
done
