#!/bin/bash
#SBATCH -J ragtag.val
#SBATCH -o logs/ragtag.val.o%A.%a
#SBATCH -e logs/ragtag.val.e%A.%a
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster
#SBATCH --array=0-2
#SBATCH --export=ALL
# Validation of the Thille-guided scaffolds (PLAN.md v16 R5). Task 0 gwen_primary, 1 gwen_alt, 2 hass_uoq_broken:
#  1. chr1..chr12 of the scaffolded assembly -> <name>.chr.fasta (+ .fai)
#  2. minimap2 2.30 -x asm5 of those chromosomes against Thille v4 and a whole-genome dot plot
#     (assemblies/thille.hic/compare/plot_genome_dotplot.py, Thille on x)
#  3. BUSCO 6.1 genome mode, embryophyta_odb12.2, offline (same image, lineage and download folder as the Thille v4 run)
# Binds: read-only Thille FASTA and plotting scripts, the BUSCO download folder; outputs under the working directory.
# Submit from assemblies/ragtag_thille_v4: sbatch sbatch_validate.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
NAMES=(gwen_primary gwen_alt hass_uoq_broken); N=${NAMES[$SLURM_ARRAY_TASK_ID]}
T=/epool/edwin/pamericana/avoflow/assemblies/thille.hic
REF=$T/thille.hifiv4.asm.hic.scaffold.chr.only.fasta
SAM=/gpool/bin/samtools-1.24/samtools; MM=/gpool/bin/minimap2-2.30_x64-linux/minimap2
AB=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
# /dpool is not mounted on the compute nodes; this copy of the same lineage is (a BUSCO database, not Hass data).
BU=/gpool/bin/singularityimages/busco6.1.sif; DL=/epool/edwin/pamericana/avoflow/annotation/tiberius/hass.hic/busco_downloads
echo "node=$(hostname) cpus=$SLURM_CPUS_PER_TASK mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} name=$N"
free -g
cd "$N"
F=$N.thille_v4_ragtag.fasta; C=$N.chr.fasta
[[ -s $C ]] || { $SAM faidx "$F"; $SAM faidx "$F" $(seq -f 'chr%g' 1 12) > "$C.tmp"; mv "$C.tmp" "$C"; $SAM faidx "$C"; }
P=$N.chr_vs_thille_v4.asm5.paf
[[ -s $P.done ]] || { $MM -x asm5 -t "$SLURM_CPUS_PER_TASK" --secondary=no "$REF" "$C" > "$P.tmp"; mv "$P.tmp" "$P"; date > "$P.done"; }
awk 'BEGIN{OFS="\t"; print "name","length","group"} {print $1,$2,$1}' "$C.fai" > $N.chr.order.tsv
mkdir -p .mpl
[[ -s ${N}_vs_thille_v4_dotplot.png ]] || singularity exec --no-home -B "$T/compare:$T/compare:ro" "$AB" bash -c "
export MPLCONFIGDIR=\$PWD/.mpl
python3 $T/compare/plot_genome_dotplot.py $P $T/compare/v4/genome_dotplots/thille_chroms.order.tsv $N.chr.order.tsv \
  ${N}_vs_thille_v4_dotplot --x-side target --x-name Thille --y-name '$N (RagTag on Thille v4)' \
  --title 'Thille v4 vs $N scaffolded on Thille v4 (RagTag 2.1.1), minimap2 asm5' --no-orient
"
export SINGULARITYENV_TINI_SUBREAPER=1
singularity exec --no-home -B "$DL:$DL:ro" "$BU" busco -i "$C" -l embryophyta_odb12.2 -o busco_${N}_chr \
  --out_path . -m genome -c "$SLURM_CPUS_PER_TASK" --download_path "$DL" --offline
grep -h "C:" busco_${N}_chr/short_summary.specific.embryophyta_odb12.2.busco_${N}_chr.txt
echo "[$(date)] done"
