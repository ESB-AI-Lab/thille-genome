#!/bin/bash
#SBATCH -J compare.prep.v4
#SBATCH -o logs/compare.prep.v4.o%j
#SBATCH -e logs/compare.prep.v4.e%j
#SBATCH -c 32
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Build inputs for the comparison dot plots (host tools only, no container):
#   1. Thille per-chr FASTAs
#   2. minimap2 Gwen contigs -> Thille (asm5 and asm10), Gwen pseudo-chromosomes per preset
#   3. whole-genome Thille and Hass concatenations (chr1..12, 100 kb N spacers)
# Submit from assemblies/thille.hic: sbatch compare/sbatch_prep_v4.sh
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
HASS=hass_UoQ_chr_only.fasta
GWEN=pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta
SPACER=100000
T=${SLURM_CPUS_PER_TASK}

echo "node=$(hostname) cpus=$T mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_JOB_ID}"
echo "minimap2 $($MINIMAP2 --version)  python $(python3 --version)"
free -g

# v4: compare/v4/thille_chrs/ already written from the same input on 2026-10-04 (11 byte-identical to v3); not rewritten.
# python3 compare/split_fasta.py "$THILLE" compare/v4/thille_chrs

mkdir -p compare/v4/gwen
for PRESET in asm5 asm10; do
  PAF=compare/v4/gwen/gwen_vs_thille.$PRESET.paf
  if [[ -s $PAF.done ]]; then
    echo "[$(date)] $PAF exists with .done, skipping minimap2"
  else
    echo "[$(date)] minimap2 -x $PRESET start"
    $MINIMAP2 -x "$PRESET" -t "$T" "$THILLE" "$GWEN" > "$PAF.tmp"
    mv "$PAF.tmp" "$PAF"
    date > "$PAF.done"
  fi
  echo "[$(date)] build pseudo-chr $PRESET"
  python3 compare/build_gwen_pseudochr.py "$PAF" "$GWEN" "compare/v4/gwen/$PRESET/gwen_$PRESET" "gwen_$PRESET" --min-mapq 20
  cat "compare/v4/gwen/$PRESET/gwen_$PRESET.summary.tsv"
done

python3 compare/concat_chroms.py "$THILLE" thille_wg $SPACER compare/v4/wg/thille_wg.fasta compare/v4/wg/thille_wg.boundaries.tsv
python3 compare/concat_chroms.py "$HASS"   hass_wg   $SPACER compare/v4/wg/hass_wg.fasta   compare/v4/wg/hass_wg.boundaries.tsv
echo "[$(date)] prep done"
