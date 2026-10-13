#!/bin/bash
#SBATCH --job-name=hapsolo.thille
#SBATCH --partition=gcluster,gcbigmem
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/hapsolo_thille.%j.log
# HapSolo2 (hapsolo2.sif, v2.0.0; minimap2 2.31, miniprot 0.18) on the Thille hifiasm Hi-C-mode primary contigs that
# went into YaHS (2,021 contigs, 989.0 Mb), to check for retained haplotigs (Dario Cantu comment #49; Edwin chose HapSolo
# over purge_dups, 2026-10-08). The assembly is not changed; the result reports how much sequence HapSolo would call
# secondary. Steps as in the image's scripts/sbatch_pipeline.sh: preprocess, self-alignment, ortholog search
# (embryophyta_odb12.2, miniprot splice model 1 = general, for plants), hill-climbing training.
# Submit from the repo root: sbatch scripts/sbatch_hapsolo_thille.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node $(hostname)  cpus ${SLURM_CPUS_PER_TASK}  start $(date)"
W=/gpool/edwin/avocados/thille_paper/hapsolo_thille
IMG=/gpool/bin/singularityimages/hapsolo2.sif
SRC=/epool/edwin/pamericana/avoflow/hifiasm/thille/hifi.hic.gcluster/thille.hifiv3.asm.hic.p_ctg.fasta
LIN=/epool/edwin/pamericana/avoflow/annotation/tiberius/hass.hic/busco_downloads/lineages/embryophyta_odb12.2
cd $W
[ -e thille_contigs.fasta ] || cp $SRC thille_contigs.fasta   # HapSolo rejects absolute input paths; a local copy
[ -e lineage ] || cp -r $LIN lineage                           # the OrthoDB lineage files only
H="singularity exec --bind $W:/work --pwd /work $IMG python3 /opt/HapSolo/hapsolo_cli.py"
T=$SLURM_CPUS_PER_TASK
[ -e thille_contigs_new.fasta ] || $H preprocess -i thille_contigs.fasta
[ -e thille_contigs_new_self_align.paf.gz ] || $H align -i thille_contigs_new.fasta -t $T
[ -d ortholog_output ] || $H search -i thille_contigs_new.fasta -l lineage -o ortholog_output -t $T -j 1
$H train -i thille_contigs_new.fasta --paf thille_contigs_new_self_align.paf.gz -b ortholog_output -t $T -n 2000
echo "done $(date)"; ls -la asms/ 2>/dev/null
