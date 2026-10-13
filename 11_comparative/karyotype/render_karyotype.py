#!/usr/bin/env python3
"""Thille karyotype figures with ChromoPlot: round centromere at the called centromere point, semicircular telomere
caps, genes (forest green) and repeats (baby blue), the called centromere region shaded as pericentromere.

Loads the annotations once per centromere call set, then renders every style x centromere color as PDF and PNG.

Usage (from assemblies/thille.hic, inside the asmbenchmarkers container):
  render_karyotype.py <chromoplot_dir> <chrom.sizes> <genes.gtf> <repeats.gff> <outdir> \\
      --calls satellitefirst:centromere/thille_centromeres.satellitefirst \\
      --calls hicfirst:centromere/thille_centromeres.hicfirst \\
      [--styles split-idiogram,dual-track,composite,banded] [--colors darkgrey:#4a4a4a,lavender:#b39ddb]
Each --calls prefix needs <prefix>.peak.bed (centromere point) and <prefix>.bed (region).

--dots TSV (columns chrom, start, end, gene, set; header line) adds small dots on the gene-side edge of each
chromosome at every gene midpoint, one row of dots per set a gene belongs to, in --dot-sets order with
--dot-colors (split-idiogram and dual-track styles; outputs get the suffix --dot-suffix). Without --dots the figures are unchanged.
"""
import argparse
import copy
import os
import re
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("chromoplot_dir")
    ap.add_argument("sizes")
    ap.add_argument("genes")
    ap.add_argument("repeats")
    ap.add_argument("outdir")
    ap.add_argument("--calls", action="append", required=True, metavar="NAME:PREFIX")
    ap.add_argument("--styles", default="split-idiogram,dual-track,composite,banded")
    ap.add_argument("--colors", default="darkgrey:#4a4a4a,lavender:#b39ddb")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--no-text", action="store_true", help="omit the figure title and the legend box title (legend text goes in the figure legend); default off")
    ap.add_argument("--dots", help="TSV of genes to mark: chrom, start, end, gene, set")
    ap.add_argument("--dot-sets", default="fst,coverage,florid", help="set names in --dots, in row order")
    ap.add_argument("--dot-labels", default="Flowering genes in top 1% Fst windows,Flowering genes selected by coverage,"
                    "FLOR-ID flowering-time genes")
    ap.add_argument("--dot-colors", default="#D55E00,#0072B2,#E69F00",
                    help="Okabe-Ito vermillion, blue, orange; validated colorblind-safe as a set")
    ap.add_argument("--dot-suffix", default="_flowering")
    a = ap.parse_args()
    if a.dots and not set(a.styles.split(",")) <= set(DOT_GEOMETRY):
        sys.exit(f"ERROR: --dots is implemented for --styles {','.join(DOT_GEOMETRY)} only")

    sys.path.insert(0, a.chromoplot_dir)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from chromoplot import (ChromosomePlotter, load_centromeres, load_features, load_pericentromeres,
                            load_scaffolds)
    from chromoplot.parser import natural_sort_key

    t0 = time.time()
    base = load_scaffolds(a.sizes)
    load_features(a.genes, feature_type="gene", chromosomes=base)
    load_features(a.repeats, feature_type="repeat", chromosomes=base)
    n_g = sum(len(c.genes) for c in base.values())
    n_r = sum(len(c.repeats) for c in base.values())
    print(f"loaded {len(base)} chromosomes, {n_g} genes, {n_r} repeats in {time.time() - t0:.0f} s", flush=True)
    os.makedirs(a.outdir, exist_ok=True)

    for spec in a.calls:
        name, prefix = spec.split(":", 1)
        chroms = copy.deepcopy(base)
        load_centromeres(prefix + ".peak.bed", chromosomes=chroms)
        load_pericentromeres(prefix + ".bed", chromosomes=chroms)
        ordered = [chroms[k] for k in sorted(chroms, key=natural_sort_key)]
        for cspec in a.colors.split(","):
            cname, color = cspec.split(":", 1)
            plotter = ChromosomePlotter(chromosomes=ordered, centromere_color=color,
                                        centromere_shape="circle", round_telomeres=True)
            for style in a.styles.split(","):
                t = time.time()
                out = os.path.join(a.outdir, f"thille_karyotype_{name}_{style}_{cname}")
                fig = plotter.plot(output_path=None if a.dots else out + ".png", style=style, layout="horizontal", dpi=a.dpi,
                                   title="" if a.no_text else f"Thille chromosomes: genes, repeats and centromeres ({name} calls)")
                if a.no_text:
                    for leg in fig.legends:
                        leg.set_title(None)
                if a.dots:
                    out += a.dot_suffix
                    add_dots(fig, ordered, a, style)
                if a.no_text:
                    strip_text(fig)
                if a.dots:
                    fig.savefig(out + ".png", dpi=a.dpi, bbox_inches="tight")
                fig.savefig(out + ".pdf", bbox_inches="tight")
                plt.close(fig)
                print(f"{out}.png/.pdf in {time.time() - t:.0f} s", flush=True)


def strip_text(fig):
    """--no-text: drop ChromoPlot's in-panel words (Pericentromere, p-arm, q-arm), reduce the chromosome labels to
    "chrN" (size and morphology go to the legend) and drop parenthetical notes and counts from the keys."""
    for ax in fig.axes:
        for t in list(ax.texts):
            if t.get_text().strip() in ("Pericentromere", "p-arm", "q-arm", "CEN"):
                t.remove()
        labels = [re.sub(r"\s*\n.*", "", lab.get_text(), flags=re.S) for lab in ax.get_yticklabels()]
        if any(labels):
            ax.set_yticks(ax.get_yticks(), labels)
        legs = ([ax.get_legend()] if ax.get_legend() else [])
        for leg in legs:
            for txt in leg.get_texts():
                txt.set_text(re.sub(r"\s*\([^)]*\)", "", txt.get_text()))
    for leg in fig.legends:
        for txt in leg.get_texts():
            txt.set_text(re.sub(r"\s*\([^)]*\)", "", txt.get_text()))


# Per style: (row spacing, first dot row above y_mid, step between set rows). ChromoPlot (plotter.py,
# _plot_horizontal) puts chromosome i at y = (n - 1 - i) * spacing. split-idiogram: outline 0.72 high, genes in its
# upper half, top edge y_mid + 0.36, pericentromere bracket y_mid + 0.65. dual-track: gene track ends at
# y_mid + 0.38 / 2 + 0.08 + 0.70 = y_mid + 0.97, bracket y_mid + 1.25.
DOT_GEOMETRY = {"split-idiogram": (2.1, 0.45, 0.09), "dual-track": (3.5, 1.06, 0.10)}


def add_dots(fig, ordered, a, style):
    """Dots on the gene-side edge of each chromosome, one row per set."""
    spacing, base, step = DOT_GEOMETRY[style]
    from matplotlib.lines import Line2D
    sets, labels, colors = a.dot_sets.split(","), a.dot_labels.split(","), a.dot_colors.split(",")
    ax = fig.axes[0]
    n = len(ordered)
    ymid = {c.name: (n - 1 - i) * spacing for i, c in enumerate(ordered)}
    seen = {}
    for line in open(a.dots):
        f = line.rstrip("\n").split("\t")
        if f[0] == "chrom":
            continue
        if f[0] not in ymid or f[4] not in sets:
            sys.exit(f"ERROR: unknown chromosome or set in --dots: {f}")
        seen.setdefault(f[4], set()).add((f[0], f[1]))
        k = sets.index(f[4])
        ax.plot((int(f[1]) + int(f[2])) / 2, ymid[f[0]] + base + step * k, marker="o", ms=4.2, ls="",
                mfc=colors[k], mec="#1f1f1f", mew=0.4, zorder=12)
    handles = [Line2D([], [], marker="o", ls="", ms=6, mfc=colors[i], mec="#1f1f1f", mew=0.4,
                      label=f"{labels[i]} ({len(seen.get(s, ()))})") for i, s in enumerate(sets)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.52, 0.915 if style == "dual-track" else 0.905), ncol=3, fontsize=9, frameon=False)
    print("dots per set: " + ", ".join(f"{s}={len(seen.get(s, ()))}" for s in sets), flush=True)


if __name__ == "__main__":
    main()
