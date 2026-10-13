#!/bin/bash
#SBATCH -J compare.puebla.prep.v4
#SBATCH -o logs/compare.puebla.prep.v4.o%j
#SBATCH -e logs/compare.puebla.prep.v4.e%j
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Puebla (GCA_060308345.1, ASM6030834v1; chromosomes CM191109.1-CM191120.1 = chr1-12) against Thille, step 1:
# approximate minimap2 asm5 alignment of Puebla chr1-12 to Thille chr1-12 (same settings as the Hass genome dot plot),
# per-chromosome majority strand and best Thille chromosome from the data, a same-strand Puebla FASTA renamed chr1..chr12
# for SyRI (chromosomes with mostly minus-strand matches reverse-complemented), and the MUMmer-style genome dot plots.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_puebla_prep_v4.sh
# Outputs: compare/v4/puebla/ (PAF, strand table, puebla.chr.samestrand.fasta), compare/v4/genome_dotplots/thille_vs_puebla_genome*
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1
MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SAMTOOLS=/gpool/bin/samtools-1.24/samtools
SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
PDIR=/gpool/edwin/avocados/thille_paper/genome/puebla/GCA_060308345.1/ncbi_dataset/data/GCA_060308345.1
PUEBLA=$PDIR/GCA_060308345.1_ASM6030834v1_genomic.fna
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
OUT=compare/v4/puebla
PLOTS=compare/v4/genome_dotplots
T=${SLURM_CPUS_PER_TASK}
echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
echo "minimap2 $($MINIMAP2 --version)  $($SAMTOOLS --version | head -1)"
free -g
mkdir -p "$OUT"

# chr1..12 only, original ids
[[ -s $OUT/puebla.chr.orig.fasta ]] || $SAMTOOLS faidx "$PUEBLA" CM191109.1 CM191110.1 CM191111.1 CM191112.1 CM191113.1 \
  CM191114.1 CM191115.1 CM191116.1 CM191117.1 CM191118.1 CM191119.1 CM191120.1 > "$OUT/puebla.chr.orig.fasta"
$SAMTOOLS faidx "$OUT/puebla.chr.orig.fasta"
python3 - "$OUT/puebla.chr.orig.fasta.fai" > "$OUT/puebla_chroms.order.tsv" <<'EOF'
import sys
print("name\tlength\tgroup")
for i, l in enumerate(open(sys.argv[1])):
    n, L = l.split()[:2]
    print(f"{n}\t{L}\tchr{i + 1}")
EOF

PAF=$OUT/puebla_vs_thille.asm5.paf
if [[ ! -s $PAF.done ]]; then
  echo "[$(date)] minimap2 -x asm5 start"
  $MINIMAP2 -x asm5 --secondary=no -t "$T" "$THILLE" "$OUT/puebla.chr.orig.fasta" > "$PAF.tmp"
  mv "$PAF.tmp" "$PAF"; date > "$PAF.done"
fi

# strand and chromosome correspondence from blocks >= 10 kb
python3 - "$PAF" "$OUT/puebla_chroms.order.tsv" > "$OUT/puebla_strand.tsv" <<'EOF'
import collections, sys
order = [l.rstrip("\n").split("\t") for l in open(sys.argv[2]) if not l.startswith("name")]
grp = {r[0]: r[2] for r in order}
m = collections.defaultdict(lambda: collections.Counter())
for l in open(sys.argv[1]):
    f = l.split("\t")
    if int(f[10]) >= 10000:
        m[f[0]][(f[5], f[4])] += int(f[9])
print("puebla_id\tpuebla_chr\tbest_thille_chr\tfrac_matches_on_best\tminus_frac_on_best\trevcomp")
for pid, _, g in [(r[0], r[1], r[2]) for r in order]:
    c = m[pid]; tot = sum(c.values())
    per = collections.Counter()
    for (t, s), v in c.items():
        per[t] += v
    best = per.most_common(1)[0][0]
    minus = c[(best, "-")] / per[best]
    print(f"{pid}\t{g}\t{best}\t{per[best] / tot:.3f}\t{minus:.4f}\t{'yes' if minus > 0.5 else 'no'}")
EOF
cat "$OUT/puebla_strand.tsv"
awk -F'\t' 'NR>1 && $2!=$3 {bad=1; print "ERROR: " $2 " best matches Thille " $3 > "/dev/stderr"} END{exit bad}' "$OUT/puebla_strand.tsv"

# same-strand FASTA renamed chr1..chr12 (input for SyRI)
QRY=$OUT/puebla.chr.samestrand.fasta
: > "$QRY.tmp"
while IFS=$'\t' read -r pid chr best frac minus rc; do
  [[ $pid == puebla_id ]] && continue
  if [[ $rc == yes ]]; then flag=-i; else flag=; fi
  $SAMTOOLS faidx $flag --mark-strand no "$OUT/puebla.chr.orig.fasta" "$pid" | sed "1s/^>.*/>$chr/" >> "$QRY.tmp"
done < "$OUT/puebla_strand.tsv"
mv "$QRY.tmp" "$QRY"; $SAMTOOLS faidx "$QRY"

singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$OUT/.mpl
python3 compare/plot_genome_dotplot.py $PAF $PLOTS/thille_chroms.order.tsv $OUT/puebla_chroms.order.tsv $PLOTS/thille_vs_puebla_genome \
  --x-side target --x-name Thille --y-name 'Puebla' --title 'Thille vs Puebla (ASM6030834v1), minimap2 asm5'
python3 compare/plot_genome_dotplot.py $PAF $PLOTS/thille_chroms.order.tsv $OUT/puebla_chroms.order.tsv $PLOTS/thille_vs_puebla_genome_raw \
  --x-side target --x-name Thille --y-name 'Puebla (as stored)' --title 'Thille vs Puebla (ASM6030834v1), minimap2 asm5' --no-orient
"
echo "[$(date)] done"
