#!/bin/bash
#SBATCH -J hic.balvar.v4
#SBATCH -o logs/hic.balvar.v4.o%j
#SBATCH -e logs/hic.balvar.v4.e%j
#SBATCH -c 8
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Default ICE balancing of the all-scaffolds map diverged at every resolution, so compare four ways of
# showing it (details in hic/hic-readme.md):
#   V1 raw only | V2 stricter bin filters | V3 blacklist scaffolds < 1 Mb | V4 cis-only balancing
# Works on a copy (hic/all/all.variants.mcool); the original all.mcool and its 'weight' column are untouched.
# v4 copy (2026-10-04): reads and writes under hic/v4/.
# Submit from assemblies/thille.hic: sbatch hic/sbatch_balance_variants_v4.sh [scope]
set -euo pipefail
export PYTHONUNBUFFERED=1 MPLCONFIGDIR=$PWD/hic/.mplconfig

SIF=/gpool/bin/singularityimages/asmbenchmarkers/asmbenchmarkers_20260925.sif
SCOPE=${1:-all}
SRC=hic/v4/$SCOPE/$SCOPE.mcool
MC=hic/v4/$SCOPE/$SCOPE.variants.mcool
RES=(50000 100000 250000 500000 1000000)
MINLEN=1000000
BL=hic/v4/$SCOPE/blacklist_lt1Mb.bed
TSV=hic/v4/$SCOPE/balance_variants.tsv
T=${SLURM_CPUS_PER_TASK}
mkdir -p "$MPLCONFIGDIR"
run() { singularity exec --no-home "$SIF" "$@"; }

echo "node=$(hostname) cpus=$T mem_per_cpu=${SLURM_MEM_PER_CPU:-NA}MB job=${SLURM_JOB_ID} scope=$SCOPE"
free -g
[[ -e $MC ]] && { echo "ERROR: $MC exists; move it aside first" >&2; exit 1; }
cp "$SRC" "$MC.tmp" && mv "$MC.tmp" "$MC"

# Whole-sequence blacklist for every sequence shorter than MINLEN.
awk -v m=$MINLEN 'BEGIN{OFS="\t"} $2 < m {print $1, 0, $2}' "hic/v4/sizes/$SCOPE.sizes" > "$BL"
echo "blacklist: $(wc -l < "$BL") sequences, $(awk '{s+=$3} END{print s}' "$BL") bp"

declare -A ARGS=(
  [weight_strict]="--min-nnz 50 --mad-max 3"
  [weight_bl1mb]="--blacklist $BL"
  [weight_cis]="--cis-only"
)
for name in weight_strict weight_bl1mb weight_cis; do
  for r in "${RES[@]}"; do
    echo "[$(date)] balance $name ${ARGS[$name]} @ $r"
    run cooler balance -p "$T" --name "$name" ${ARGS[$name]} "$MC::resolutions/$r"
  done
done

run python3 - "$MC" "$TSV" "${RES[@]}" <<'EOF'
import sys, cooler, numpy as np
mc, out, res = sys.argv[1], sys.argv[2], [int(r) for r in sys.argv[3:]]
with open(out, "w") as fh:
    fh.write("weight\tresolution\tconverged\tvar\tmasked_bins\ttotal_bins\n")
    for name in ("weight", "weight_strict", "weight_bl1mb", "weight_cis"):
        for r in res:
            c = cooler.Cooler(f"{mc}::resolutions/{r}")
            with c.open("r") as h:
                a, w = dict(h[f"bins/{name}"].attrs), h[f"bins/{name}"][:]
            # --cis-only stores converged/var per chromosome (arrays): converged only if all did; report the max var.
            conv = bool(np.all(a.get("converged")))
            var = float(np.nanmax(np.atleast_1d(a.get("var", np.nan))))
            fh.write(f"{name}\t{r}\t{conv}\t{var:.3e}\t{int(np.isnan(w).sum())}\t{len(w)}\n")
EOF
cat "$TSV"

PS=hic/v4/$SCOPE/$SCOPE.pairstats.tsv
P="python3 hic/plot_contacts.py $MC $SCOPE hic/v4/figures $PS --genome-only"
run $P --raw-only --suffix _V1_raw-only \
  --variant-label "V1: raw counts only (default ICE balancing did not converge for this scope)"
run $P --weight-name weight_strict --suffix _V2_strict-filters \
  --variant-label "V2: ICE with stricter bin filters (--min-nnz 50 --mad-max 3)"
run $P --weight-name weight_bl1mb --suffix _V3_blacklist-lt1Mb \
  --variant-label "V3: ICE with all scaffolds < 1 Mb blacklisted (shown white in balanced panel)"
run $P --weight-name weight_cis --suffix _V4_cis-only \
  --variant-label "V4: ICE weights from intra-scaffold (cis) contacts only (--cis-only)"
echo "[$(date)] done"
