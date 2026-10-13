#!/bin/bash
#SBATCH -J inv.plot
#SBATCH -o logs/inv.plot.o%j
#SBATCH -e logs/inv.plot.e%j
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --export=ALL
# Draw the 152 inversion breakpoint figures (compare/v4/inversions/plot_inversions.py). Read-only binds: the v4
# flowering BEDs, the TE library headers and the Hass UoQ gene BED directory. --mem: the Thille TE GFF is read whole.
# Submit from assemblies/thille.hic: sbatch compare/v4/inversions/sbatch_plot_inversions.sh
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/hic/.mplconfig
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
F=/gpool/edwin/avocados/thille_paper/flowering_beds/v4/igv
L=/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/consensi.fa.classified
H=/epool/edwin/pamericana/avoflow/z_GCA_029852735.1/ncbi_dataset/data/GCA_029852735.1
singularity exec --no-home -B $F:$F:ro -B $L:$L:ro -B $H:$H:ro \
  /gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif python3 compare/v4/inversions/plot_inversions.py --outdir figures_v2
echo "[$(date)] done"
