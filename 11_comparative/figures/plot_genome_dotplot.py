#!/usr/bin/env python3
"""MUMmer-style whole-genome alignment dot plot from a minimap2 PAF.

Both genomes are laid end to end in the order given by an order file per axis. A thin grey grid line marks
every fragment boundary (chromosome/scaffold/contig) and a black line marks every group boundary (for example
all Gwen contigs placed on Thille chr3 form group chr3); ticks label the groups. Alignment blocks are drawn
as segments, blue for '+' and red for '-'.

Orientation (default, like mummerplot --layout): each y-axis fragment whose aligned matches (PAF col 10) are mostly
on the '-' strand is drawn reverse-complemented, because a whole fragment stored on the other strand is a strand
convention of that assembly, not an inversion. Only blocks against the fragment's majority strand stay red, so red
means a real internal inversion. Flipped fragments are listed in <out_prefix>.flipped.tsv and single-fragment groups
are labelled "(rc)". --no-orient draws the raw PAF.

Order file (TSV with header): name, length, group — rows in axis order. Every PAF name on that axis that is not
in the order file is ignored (counted in the summary).

Usage:
  plot_genome_dotplot.py <aln.paf> <x_order.tsv> <y_order.tsv> <out_prefix> --x-side {query,target}
                         [--x-name X] [--y-name Y] [--min-block 10000] [--min-mapq 0] [--title TEXT] [--no-orient]
Writes <out_prefix>.pdf and .png (300 dpi).
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

PLUS, MINUS, INK, GRID = "#1f5fa8", "#c0392b", "#1f1f1f", "#c8c8c8"
plt.rcParams.update({"font.size": 8, "text.color": INK, "pdf.fonttype": 42})


def load_order(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    off, pos, groups = {}, 0, []
    for r in rows:
        off[r["name"]] = pos
        n = int(r["length"])
        if groups and groups[-1][0] == r["group"]:
            groups[-1][2] = pos + n
        else:
            groups.append([r["group"], pos, pos + n])
        pos += n
    return off, [off[r["name"]] for r in rows] + [pos], groups, pos, len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf")
    ap.add_argument("x_order")
    ap.add_argument("y_order")
    ap.add_argument("out_prefix")
    ap.add_argument("--x-side", choices=("query", "target"), required=True, help="which PAF side goes on the x axis")
    ap.add_argument("--x-name", default="x")
    ap.add_argument("--y-name", default="y")
    ap.add_argument("--min-block", type=int, default=10000, help="minimum alignment block length (PAF col 11), default 10 kb")
    ap.add_argument("--min-mapq", type=int, default=0, help="minimum MAPQ (default 0 = all)")
    ap.add_argument("--title", default="")
    ap.add_argument("--no-orient", action="store_true", help="draw the raw PAF strands (no majority-strand flipping)")
    ap.add_argument("--no-text", action="store_true", help="omit the title and explanatory notes (legend text goes in the figure legend); default off")
    a = ap.parse_args()
    xoff, xb, xg, xtot, xn = load_order(a.x_order)
    yoff, yb, yg, ytot, yn = load_order(a.y_order)
    ylen = {r["name"]: int(r["length"]) for r in csv.DictReader(open(a.y_order), delimiter="\t")}
    ygroup = {r["name"]: r["group"] for r in csv.DictReader(open(a.y_order), delimiter="\t")}

    blocks, skipped = [], 0
    with open(a.paf) as fh:
        for line in fh:
            f = line.split("\t")
            if int(f[10]) < a.min_block or int(f[11]) < a.min_mapq:
                continue
            q = (f[0], int(f[2]), int(f[3]))
            t = (f[5], int(f[7]), int(f[8]))
            x, y = (q, t) if a.x_side == "query" else (t, q)
            if x[0] not in xoff or y[0] not in yoff:
                skipped += 1
                continue
            blocks.append((x, y, f[4], int(f[9])))

    # Majority strand per y fragment, by residue matches over all its drawn blocks.
    tally = {}
    for _, y, st, m in blocks:
        tally.setdefault(y[0], {"+": 0, "-": 0})[st] += m
    flipped = set() if a.no_orient else {n for n, v in tally.items() if v["-"] > v["+"]}

    segs = {"+": [], "-": []}
    for x, y, st, m in blocks:
        x0, x1 = xoff[x[0]] + x[1], xoff[x[0]] + x[2]
        ys, ye = y[1], y[2]
        if y[0] in flipped:
            ys, ye = ylen[y[0]] - y[2], ylen[y[0]] - y[1]
            st = "+" if st == "-" else "-"
        y0, y1 = yoff[y[0]] + ys, yoff[y[0]] + ye
        # For '-' the y sequence runs backwards against x.
        if st == "-":
            y0, y1 = y1, y0
        segs[st].append([(x0 / 1e6, y0 / 1e6), (x1 / 1e6, y1 / 1e6)])
    kept = len(blocks)
    with open(f"{a.out_prefix}.flipped.tsv", "w") as out:
        out.write("fragment\tgroup\tlength\tmatch_plus_bp\tmatch_minus_bp\tflipped\n")
        for n in sorted(tally, key=lambda n: yoff[n]):
            out.write(f"{n}\t{ygroup[n]}\t{ylen[n]}\t{tally[n]['+']}\t{tally[n]['-']}\t{'yes' if n in flipped else 'no'}\n")
    single = {}
    for n, g in ygroup.items():
        single.setdefault(g, []).append(n)

    fig, ax = plt.subplots(figsize=(9, 8.6), constrained_layout=True)
    for b in xb:
        ax.axvline(b / 1e6, color=GRID, lw=0.25, zorder=0)
    for b in yb:
        ax.axhline(b / 1e6, color=GRID, lw=0.25, zorder=0)
    for _, s, e in xg:
        for b in (s, e):
            ax.axvline(b / 1e6, color="black", lw=0.6, zorder=1)
    for _, s, e in yg:
        for b in (s, e):
            ax.axhline(b / 1e6, color="black", lw=0.6, zorder=1)
    for st, col in (("+", PLUS), ("-", MINUS)):
        ax.add_collection(LineCollection(segs[st], colors=col, linewidths=1.0, capstyle="butt", zorder=2))
    ax.set_xticks([(s + e) / 2e6 for _, s, e in xg], [g for g, _, _ in xg], rotation=90)
    ax.set_yticks([(s + e) / 2e6 for _, s, e in yg],
                  [g + (" (rc)" if len(single[g]) == 1 and single[g][0] in flipped else "") for g, _, _ in yg])
    ax.tick_params(length=0)
    ax.set_xlim(0, xtot / 1e6)
    ax.set_ylim(0, ytot / 1e6)
    ax.set_xlabel(a.x_name if a.no_text else f"{a.x_name}: {xn} fragments, {xtot / 1e6:.0f} Mb")
    ax.set_ylabel(a.y_name if a.no_text else f"{a.y_name}: {yn} fragments, {ytot / 1e6:.0f} Mb")
    if a.no_orient:
        ax.plot([], [], color=PLUS, lw=2, label="same strand (+)")
        ax.plot([], [], color=MINUS, lw=2, label="opposite strand (−)")
    else:
        ax.plot([], [], color=PLUS, lw=2, label="collinear")
        ax.plot([], [], color=MINUS, lw=2, label="inverted relative to rest of fragment")
    ax.legend(loc="upper left", frameon=True, fontsize=7)
    t = a.title or f"{a.x_name} vs {a.y_name}"
    orient_txt = ("raw PAF strands" if a.no_orient else
                  f"{len(flipped)} of {len(tally)} y-axis fragments reverse-complemented to their majority strand "
                  "(whole-fragment reversal = strand convention)")
    if not a.no_text:
        fig.suptitle(f"{t}\nminimap2 blocks ≥ {a.min_block // 1000} kb"
                     + (f", MAPQ ≥ {a.min_mapq}" if a.min_mapq else "")
                     + "; grey = fragment boundaries, black = chromosome boundaries"
                     + f"\n{orient_txt}", color=INK, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out_prefix}.{ext}", dpi=300)
    plt.close(fig)
    print(f"wrote {a.out_prefix}.pdf/.png: {kept} blocks drawn, {skipped} skipped (sequence not in an order file)")


if __name__ == "__main__":
    main()
