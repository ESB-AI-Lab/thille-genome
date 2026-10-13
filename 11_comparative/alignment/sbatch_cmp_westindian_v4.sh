#!/bin/bash
#SBATCH -J compare.westindian.v4
#SBATCH -o logs/compare.westindian.v4.o%j
#SBATCH -e logs/compare.westindian.v4.e%j
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# West Indian T2T avocado (Yang et al. 2024; GenBank GCA_051132755.1, XTBG_Pame_1.0, Pa01-Pa12) -> Thille chr1-12,
# base-level alignment (-c) like the Gwen run, then a MUMmer-style whole-genome dot plot with Thille on x and
# West Indian chromosomes on y, each flipped to its majority strand. --mem because the unconstrained Gwen/Hass
# alignments ran out of memory at the 20 GB default.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_cmp_westindian_v4.sh
# Outputs: compare/v4/westindian/westindian_vs_thille.asm5.c.paf, compare/v4/genome_dotplots/thille_vs_westindian_genome*.
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
WI=westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna
SEQREP=westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/sequence_report.jsonl
T=${SLURM_CPUS_PER_TASK}
OUT=compare/v4/westindian
PAF=$OUT/westindian_vs_thille.asm5.c.paf
PLOTS=compare/v4/genome_dotplots

echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
echo "minimap2 $($MINIMAP2 --version)"
free -g
mkdir -p "$OUT" "$PLOTS"

# Order file: West Indian chromosomes in Pa01..Pa12 order, grouped by their own names.
python3 - "$SEQREP" > "$OUT/westindian_chroms.order.tsv" <<'EOF'
import json, sys
rows = [json.loads(l) for l in open(sys.argv[1])]
rows.sort(key=lambda r: r["chrName"])
print("name\tlength\tgroup")
for r in rows:
    print(f'{r["genbankAccession"]}\t{r["length"]}\t{r["chrName"]}')
EOF
cat "$OUT/westindian_chroms.order.tsv"

if [[ -s $PAF.done ]]; then
  echo "[$(date)] $PAF exists with .done, skipping minimap2"
else
  echo "[$(date)] minimap2 -c -x asm5 start"
  $MINIMAP2 -c -x asm5 --secondary=no -t "$T" "$THILLE" "$WI" > "$PAF.tmp"
  mv "$PAF.tmp" "$PAF"
  date > "$PAF.done"
fi
echo "[$(date)] PAF lines: $(wc -l < "$PAF")"

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
for mode in orient raw; do
  extra=''; suffix=''
  if [[ \$mode == raw ]]; then extra='--no-orient'; suffix='_raw'; fi
  python3 compare/plot_genome_dotplot.py $PAF $PLOTS/thille_chroms.order.tsv $OUT/westindian_chroms.order.tsv \
    $PLOTS/thille_vs_westindian_genome\$suffix --x-side target --x-name Thille --y-name 'West Indian (T2T)' \
    --title 'Thille vs West Indian T2T (XTBG_Pame_1.0), minimap2 asm5 -c' \$extra
done
"
echo "[$(date)] done"
