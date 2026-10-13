#!/usr/bin/env python3
"""Summarize read support at the inversion breakpoints (support/task*.tsv from inv_pileup.py).

Per inversion, read set and genome, the spanning and split-read counts of the two breakpoints are summed and scored:
  collinear     >= MIN spanning reads and < MIN split reads: the reads follow that assembly's orientation
  inverted      >= MIN split reads and < MIN spanning reads: the reads carry the other orientation
  both          >= MIN of each (expected where the read donor is heterozygous for the inversion)
  low coverage  neither reaches MIN
The split fraction split / (split + spanning) is the overview color, in five bands with 0.5 highlighted. For inversions shorter than the 5 kb spanning
flank a read can cross the whole inversion, so "collinear" there means no split alignment was made, which minimap2
may also do for a short inverted segment (it can align through it with mismatches); treat short inversions with care.
Writes support_summary_v2.tsv and figures_v2/overview_split_fraction.{pdf,png} (with the Puebla read sets).

Usage (inside the asmbenchmarkers container, from assemblies/thille.hic): summarize_support.py [--min 3]
"""
import argparse
import os
import collections
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402

D = "compare/v4/inversions"
INK, INK2, MUTED, BASE, NA = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#f0efec"
NAME = {"hass": "Hass UoQ", "westindian": "Unnamed Lowland", "puebla": "Puebla"}
RS = [("thille_hifi", "Thille HiFi"), ("lowland_hifi", "Unnamed Lowland HiFi"), ("gwen_clr", "Gwen CLR"),
      ("puebla_hifi", "Puebla HiFi")]
# Five bands so that 0.5 (both orientations, one per haplotype) has its own color. Anchors from the dataviz reference
# palette (blue #2a78d6, yellow #eda100, red #e34948; validated all-pairs for color-vision deficiency) with lighter
# steps between; cells called "both" also carry a dot, because yellow is below 3:1 contrast on white.
BOUNDS = [0, 0.15, 0.35, 0.65, 0.85, 1.0]
CMAP = ListedColormap(["#2a78d6", "#9ec5f4", "#eda100", "#f2a3a2", "#e34948"])
NORM = BoundaryNorm(BOUNDS, CMAP.N)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5, "pdf.fonttype": 42, "text.color": INK})


def verdict(span, split, m):
    if span >= m and split >= m:
        return "both"
    if split >= m:
        return "inverted"
    if span >= m:
        return "collinear"
    return "low coverage"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min", type=int, default=3)
    ap.add_argument("--no-text", action="store_true", help="overview without notes and long titles (text goes in the legend); default off")
    ap.add_argument("--fig-out", help="write only the overview figure, to this prefix (support_summary_v2.tsv and figures_v2/ are not touched)")
    a = ap.parse_args()
    inv = [l.rstrip("\n").split("\t") for l in list(open(f"{D}/inversions.tsv"))[1:]]
    tot = collections.defaultdict(lambda: [0, 0])
    for p in glob.glob(f"{D}/support/task*.tsv"):
        for l in list(open(p))[1:]:
            f = l.rstrip("\n").split("\t")
            k = (f[0], "thille" if f[2] == "thille" else "other", f[7])
            tot[k][0] += int(f[9])
            tot[k][1] += int(f[11])
    cols = [(g, rs, lab) for g in ("thille", "other") for rs, lab in RS]
    with open(os.devnull if a.fig_out else f"{D}/support_summary_v2.tsv", "w") as fh:
        fh.write("inv_id\tpair\tthille_len\tother_len\t" + "\t".join(
            f"{rs}_on_{g}_{x}" for g, rs, _ in cols for x in ("span", "split", "call")) + "\n")
        for r in inv:
            vals = []
            for g, rs, _ in cols:
                if (r[0], g, rs) in tot:
                    sp, sl = tot[(r[0], g, rs)]
                    vals += [sp, sl, verdict(sp, sl, a.min)]
                else:
                    vals += ["", "", "no reads"]
            fh.write("\t".join(map(str, r[:2] + [r[5], r[9]] + vals)) + "\n")
    # overview: one panel per compared genome, rows = inversions by Thille size, cells = split fraction
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 10.5), gridspec_kw=dict(wspace=0.75))
    for ax, pair in zip(axes, ("hass", "westindian", "puebla")):
        rows = sorted((r for r in inv if r[1] == pair), key=lambda r: -int(r[5]))
        for i, r in enumerate(rows):
            for j, (g, rs, _) in enumerate(cols):
                k = (r[0], g, rs)
                if k not in tot or sum(tot[k]) < a.min:
                    color = NA
                else:
                    sp, sl = tot[k]
                    color = CMAP(NORM(sl / (sp + sl)))
                ax.add_patch(plt.Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, color=color, lw=0))
                if k in tot and verdict(*tot[k], a.min) == "both":
                    ax.plot(j + 0.5, i + 0.5, marker="o", ms=1.8, color=INK, mew=0)
        ax.set_xlim(0, len(cols))
        ax.set_ylim(len(rows), 0)
        ax.set_yticks([i + 0.5 for i in range(len(rows))])
        ax.set_yticklabels([f"{r[0].split('_', 1)[1]}  {int(r[5]) / 1e3:,.1f} kb" for r in rows], fontsize=5.2,
                           color=INK2)
        ax.set_xticks([j + 0.5 for j in range(len(cols))])
        ax.set_xticklabels([f"{lab} on {'Thille' if g == 'thille' else NAME[pair]}" for g, _, lab in cols],
                           rotation=60, ha="right", fontsize=5.6, color=INK2)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.axvline(len(RS), color=BASE, lw=0.8)
        ax.set_title(NAME[pair] if a.no_text else f"Thille vs {NAME[pair]} ({len(rows)} inversions)", fontsize=7.5, loc="left", color=INK,
                     fontweight="bold")
    cax = fig.add_axes([0.08, 0.045, 0.3, 0.01])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=BOUNDS)
    cb.set_label("split fraction" if a.no_text else "split reads / (split + spanning reads) at the two breakpoints", fontsize=6.5, color=INK2)
    for x, t in ([] if a.no_text else ((0.075, "follows\nassembly"), (0.5, "both orientations\n(~0.5: one haplotype each)"),
                 (0.925, "other\norientation"))):
        cax.text(x, 2.4, t, transform=cax.transAxes, ha="center", va="bottom", fontsize=5.8, color=INK2)
    cb.outline.set_visible(False)
    if not a.no_text:
      fig.text(0.45, 0.05, f"Dot: >= {a.min} spanning and >= {a.min} split reads (both orientations present).\nGray: "
             f"fewer than {a.min} informative reads, or no alignments of that read set to that genome.", fontsize=6.2,
             color=INK2, va="top")
    fig.subplots_adjust(left=0.08, right=0.99, top=0.965, bottom=0.215)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.fig_out}.{ext}" if a.fig_out else f"{D}/figures_v2/overview_split_fraction.{ext}", dpi=200, facecolor="white")
    calls = collections.Counter()
    for r in inv:
        for g, rs, _ in cols:
            k = (r[0], g, rs)
            if k in tot:
                calls[(r[1], g, rs, verdict(*tot[k], a.min))] += 1
    for k in sorted(calls):
        print("\t".join(k), calls[k], sep="\t")


if __name__ == "__main__":
    main()
