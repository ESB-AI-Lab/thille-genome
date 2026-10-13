#!/bin/bash
#SBATCH -J notext_karyo
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -o compare/v4/notext/notext_karyo.%j.log
# Plan v18: Figure 5C redrawn after --no-text also strips ChromoPlot's in-panel labels. Submit from assemblies/thille.hic.
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/compare/v4/notext/.mpl
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
CP=/epool/edwin/pamericana/avoflow/annotation/annotation_plots
[[ -e compare/v4/notext/karyotype ]] && { echo "exists"; exit 1; }
singularity exec --no-home -B $CP:$CP:ro "$SIF" python3 karyotype/render_karyotype.py "$CP" hic/v4/sizes/chr12.sizes \
  thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff \
  compare/v4/notext/karyotype --calls satellitefirst:centromere/v4/thille_centromeres.satellitefirst --styles split-idiogram \
  --colors lavender:#b39ddb --dots fst_check/karyotype_dots.tsv --no-text
