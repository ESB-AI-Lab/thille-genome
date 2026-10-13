#!/bin/bash
#SBATCH -J compare.wi.wg.v4
#SBATCH -o logs/compare.wi.wg.v4.o%j
#SBATCH -e logs/compare.wi.wg.v4.e%j
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Whole-genome ModDotPlot, Thille vs West Indian T2T (GCA_051132755.1), same settings as the Hass run
# (sbatch_cmp_hass_wg.sh): chr1..12 concatenated in order with 100 kb N spacers, --compare-only, defaults otherwise.
# West Indian Pa01,04,05,06,10,11 are reverse-complemented first: 98.8-99.95% of their minimap2 matches to Thille
# are on the minus strand (compare/genome_dotplots/thille_vs_westindian_genome.flipped.tsv), i.e. strand convention.
# --mem: the Hass whole-genome run peaked at 109 GB.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_cmp_westindian_wg_v4.sh
# Outputs: compare/v4/wg_wi/ (concatenated FASTA, boundaries, ModDotPlot bedpe/plots, redrawn figure, orientation TSV).
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
WI=westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna
OUT=compare/v4/wg_wi
QRY=compare/v4/wg/thille_wg.fasta
REF=$OUT/westindian_wg_rc.fasta
B=$OUT/plots/thille_wg_westindian_wg_rc/thille_wg_westindian_wg_rc_COMPARE.bedpe

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
free -g
mkdir -p "$OUT/plots"
[[ -s $REF ]] || python3 compare/concat_chroms.py "$WI" westindian_wg_rc 100000 "$REF" "$OUT/westindian_wg_rc.boundaries.tsv" \
  --revcomp chr1,chr4,chr5,chr6,chr10,chr11

/usr/bin/time -v singularity exec --no-home "$SIF" \
  moddotplot static -f "$QRY" "$REF" --compare-only -o "$OUT/plots"

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
python3 compare/plot_wg_compare.py $B compare/v4/wg/thille_wg.boundaries.tsv $OUT/westindian_wg_rc.boundaries.tsv \
  $OUT/thille_vs_westindian_wg_rc_boundaries --y-name 'West Indian T2T' \
  --note 'West Indian Pa01,04,05,06,10,11 reverse-complemented (strand convention)'
python3 compare/wg_orientation.py $B compare/v4/wg/thille_wg.boundaries.tsv $OUT/westindian_wg_rc.boundaries.tsv \
  $OUT/westindian_orientation_after_rc.tsv
"
cat "$OUT/westindian_orientation_after_rc.tsv"
echo "[$(date)] done"
