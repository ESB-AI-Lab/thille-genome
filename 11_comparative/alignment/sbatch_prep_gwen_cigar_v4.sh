#!/bin/bash
#SBATCH -J compare.prep.gwen_c.v4
#SBATCH -o logs/compare.prep.gwen_c.v4.o%A.%a
#SBATCH -e logs/compare.prep.gwen_c.v4.e%A.%a
#SBATCH -c 32
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --array=0-1
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Gwen contigs -> Thille with base-level alignment (-c), one array task per preset (0 asm5, 1 asm10).
# Without -c minimap2 only chains seeds, and asm5/asm10 share the same seeding/chaining, so their
# approximate PAFs are identical; the presets only differ in alignment scoring, which -c applies.
# Submit from assemblies/thille.hic: sbatch compare/sbatch_prep_gwen_cigar_v4.sh
# Outputs: compare/v4/gwen_c/gwen_vs_thille.<preset>.c.paf, compare/v4/gwen_c/<preset>/gwen_c_<preset>.*
# v4 copy (2026-10-04): Thille v4 chr-only input, outputs under compare/v4/.
set -euo pipefail
export PYTHONUNBUFFERED=1

MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
THILLE=thille.hifiv4.asm.hic.scaffold.chr.only.fasta
GWEN=pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta
PRESETS=(asm5 asm10)
PRESET=${PRESETS[$SLURM_ARRAY_TASK_ID]}
T=${SLURM_CPUS_PER_TASK}
OUT=compare/v4/gwen_c
PAF=$OUT/gwen_vs_thille.$PRESET.c.paf

echo "node=$(hostname) cpus=$T mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} preset=$PRESET"
echo "minimap2 $($MINIMAP2 --version)  python $(python3 --version)"
free -g
mkdir -p "$OUT"

if [[ -s $PAF.done ]]; then
  echo "[$(date)] $PAF exists with .done, skipping minimap2"
else
  echo "[$(date)] minimap2 -c -x $PRESET start"
  $MINIMAP2 -c -x "$PRESET" -t "$T" "$THILLE" "$GWEN" > "$PAF.tmp"
  mv "$PAF.tmp" "$PAF"
  date > "$PAF.done"
fi
echo "[$(date)] build pseudo-chr $PRESET"
python3 compare/build_gwen_pseudochr.py "$PAF" "$GWEN" "$OUT/$PRESET/gwen_c_$PRESET" "gwen_c_$PRESET" --min-mapq 20
cat "$OUT/$PRESET/gwen_c_$PRESET.summary.tsv"
echo "[$(date)] done"
