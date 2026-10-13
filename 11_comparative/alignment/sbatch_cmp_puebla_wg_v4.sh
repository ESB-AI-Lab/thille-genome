#!/bin/bash
#SBATCH -J compare.puebla.wg.v4
#SBATCH -o logs/compare.puebla.wg.v4.o%j
#SBATCH -e logs/compare.puebla.wg.v4.e%j
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -p gcbigmem,gcluster
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Whole-genome ModDotPlot, Thille v4 vs Puebla (GCA_060308345.1), same settings as the Hass and Unnamed Lowland runs
# (sbatch_cmp_westindian_wg_v4.sh): chr1..12 concatenated in order with 100 kb N spacers, --compare-only, defaults
# otherwise (audit A4, 2026-10-07). Input is the same-strand Puebla FASTA from sbatch_puebla_prep_v4.sh, in which
# chr1, 4, 5, 6, 10, 11 are already reverse-complemented (98.2-100% of their matches to Thille on the minus strand,
# compare/v4/puebla/puebla_strand.tsv), so no --revcomp here.
# --mem: the Hass whole-genome run peaked at 109 GB.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_cmp_puebla_wg_v4.sh
# Outputs: compare/v4/wg_puebla/
set -euo pipefail
export PYTHONUNBUFFERED=1

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
PUEBLA=compare/v4/puebla/puebla.chr.samestrand.fasta
OUT=compare/v4/wg_puebla
QRY=compare/v4/wg/thille_wg.fasta
REF=$OUT/puebla_wg_rc.fasta
B=$OUT/plots/thille_wg_puebla_wg_rc/thille_wg_puebla_wg_rc_COMPARE.bedpe

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
free -g
[[ -e $OUT/plots ]] && { echo "exists: $OUT/plots"; exit 1; }
mkdir -p "$OUT/plots"
[[ -s $REF ]] || python3 compare/concat_chroms.py "$PUEBLA" puebla_wg_rc 100000 "$REF" "$OUT/puebla_wg_rc.boundaries.tsv"

/usr/bin/time -v singularity exec --no-home "$SIF" \
  moddotplot static -f "$QRY" "$REF" --compare-only -o "$OUT/plots"

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
python3 compare/plot_wg_compare.py $B compare/v4/wg/thille_wg.boundaries.tsv $OUT/puebla_wg_rc.boundaries.tsv \
  $OUT/thille_vs_puebla_wg_rc_boundaries --y-name 'Puebla' \
  --note 'Puebla chr1, 4, 5, 6, 10, 11 reverse-complemented (strand convention)'
python3 compare/wg_orientation.py $B compare/v4/wg/thille_wg.boundaries.tsv $OUT/puebla_wg_rc.boundaries.tsv \
  $OUT/puebla_orientation_after_rc.tsv
"
cat "$OUT/puebla_orientation_after_rc.tsv"
echo "[$(date)] done"
