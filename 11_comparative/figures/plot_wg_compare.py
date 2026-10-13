#!/usr/bin/env python3
"""Redraw a ModDotPlot --compare-only whole-genome plot from its bedpe, with chromosome boundaries.

The bedpe holds one row per window pair (query = x axis, reference = y axis) and its identity (%).
Windows are drawn as one image at native window resolution, colored like ModDotPlot's default
(Spectral_11, 11 equal bins from the identity cutoff to 100%, high identity = dark red).
Black lines mark every chromosome start and end on both axes; ticks at chromosome midpoints.
Chromosomes whose boundaries TSV has strand '-' are labelled 'chrN (rc)'.

Usage:
  plot_wg_compare.py <compare.bedpe> <x.boundaries.tsv> <y.boundaries.tsv> <out_prefix>
                     [--x-name Thille] [--y-name "Hass UoQ"] [--note TEXT] [--min-identity 86]
Writes <out_prefix>.pdf (image embedded 1 pixel per window) and <out_prefix>.png (300 dpi).
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import BoundaryNorm  # noqa: E402

INK, MUTED = "#1f1f1f", "#5f5f5f"
plt.rcParams.update({"font.size": 8, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "pdf.fonttype": 42})


def load_bounds(path):
    rows = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        rows.append((int(r["chrom_number"]), int(r["start_0based"]), int(r["end"]), r.get("strand", "+")))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bedpe")
    ap.add_argument("x_bounds")
    ap.add_argument("y_bounds")
    ap.add_argument("out_prefix")
    ap.add_argument("--x-name", default="Thille")
    ap.add_argument("--y-name", default="Hass UoQ")
    ap.add_argument("--note", default="", help="extra title line")
    ap.add_argument("--min-identity", type=float, default=86.0,
                    help="identity cutoff used in the ModDotPlot run (default 86, ModDotPlot's default)")
    ap.add_argument("--no-text", action="store_true", help="omit the title and explanatory notes (legend text goes in the figure legend); default off")
    a = ap.parse_args()

    xs, xe, ys, ye, ident = [], [], [], [], []
    with open(a.bedpe) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            xs.append(int(f[1])); xe.append(int(f[2])); ys.append(int(f[4])); ye.append(int(f[5]))
            ident.append(float(f[6]))
    xs, xe, ys, ye, ident = map(np.asarray, (xs, xe, ys, ye, ident))
    if ident.size == 0:
        raise SystemExit(f"ERROR: no rows in {a.bedpe}")
    w = int(np.median(xe - xs + 1))
    xb, yb = load_bounds(a.x_bounds), load_bounds(a.y_bounds)
    xlen, ylen = max(e for _, _, e, _ in xb), max(e for _, _, e, _ in yb)
    nx, ny = int(np.ceil(xlen / w)), int(np.ceil(ylen / w))
    img = np.full((ny, nx), np.nan)
    ci, ri = np.minimum((xs - 1) // w, nx - 1), np.minimum((ys - 1) // w, ny - 1)
    img[ri, ci] = ident

    cmap = plt.get_cmap("Spectral_r", 11).copy()
    cmap.set_bad("#ffffff")
    norm = BoundaryNorm(np.linspace(a.min_identity, 100, 12), cmap.N)

    gb = 1e9
    fig, ax = plt.subplots(figsize=(8.6, 8), constrained_layout=True)
    im = ax.imshow(img, cmap=cmap, norm=norm, interpolation="none", origin="lower",
                   extent=(0, nx * w / gb, 0, ny * w / gb), aspect="auto")
    for _, s, e, _ in xb:
        for b in (s, e):
            ax.axvline(b / gb, color="black", lw=0.6)
    for _, s, e, _ in yb:
        for b in (s, e):
            ax.axhline(b / gb, color="black", lw=0.6)
    ax.set_xticks([(s + e) / 2 / gb for _, s, e, _ in xb], [f"chr{n}" for n, *_ in xb], rotation=90)
    ax.set_yticks([(s + e) / 2 / gb for _, s, e, _ in yb],
                  [f"chr{n}" + (" (rc)" if st == "-" else "") for n, _, _, st in yb])
    ax.tick_params(length=0)
    ax.set_xlim(0, xlen / gb)
    ax.set_ylim(0, ylen / gb)
    sfx = "" if a.no_text else " (chr1 to chr12 joined, 100 kb N spacers)"
    ax.set_xlabel(f"{a.x_name}{sfx}")
    ax.set_ylabel(f"{a.y_name}{sfx}")
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.01)
    cb.set_label("window identity (%)", color=MUTED)
    title = f"{a.x_name} vs {a.y_name}: ModDotPlot identity, {w / 1e3:.0f} kb windows (≥{a.min_identity:g}%)"
    if not a.no_text:
        fig.suptitle(title + (f"\n{a.note}" if a.note else ""), color=INK)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out_prefix}.{ext}", dpi=300)
    plt.close(fig)
    print(f"wrote {a.out_prefix}.png/.pdf: {ident.size} windows, image {nx}x{ny}, window {w} bp")


if __name__ == "__main__":
    main()
