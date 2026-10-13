#!/usr/bin/env python3
"""Whole-genome dot plot of a SyRI run with its structural-rearrangement calls overlaid.

Draws every SyRI alignment row (SYNAL, INVAL, TRANSAL, INVTRAL, DUPAL, INVDPAL) from <name>syri.out as a segment,
Thille (reference) on x and the query on y, both chr1..chr12 in numeric order. The query is the same-strand FASTA
SyRI was run on, so whole chromosomes stored on the other strand are already flipped. Then:
  events        : SyRI inversions, translocations, duplications, deletions, insertions and copy gains/losses
                  >= --min-size. --marker box (boxes spanning the event on both genomes, enlarged to at least
                  --min-box of each axis so small events stay visible) or --marker circle (default; open circles at the
                  event midpoint, size grows with event length up to 2 Mb)
  red strips    : NOTAL regions >= --min-size, under the x axis for Thille-only sequence and beside the y axis for
                  query-only sequence

Usage: plot_syri_overlay.py <syri.out> <ref.fai> <qry.fai> <out_prefix> --y-name NAME [--min-size 100000]
                             [--marker circle|box] [--min-box 0.004]
Writes <out_prefix>.pdf and .png (300 dpi).
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

INK, ALN, ABSENT = "#1f1f1f", "#4a78b5", "#c0392b"
EVENT = {"INV": ("#e08a00", "inversion"), "TRANS": ("#2e9e44", "translocation"),
         "INVTR": ("#2e9e44", "translocation"), "DUP": ("#8e44ad", "duplication"),
         "INVDP": ("#8e44ad", "duplication"), "DEL": ("#1f1f1f", "deletion"), "INS": ("#17a2b8", "insertion"),
         "CPG": ("#8c564b", "copy gain or loss"), "CPL": ("#8c564b", "copy gain or loss")}
ALIGN_ROWS = {"SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"}
plt.rcParams.update({"font.size": 8, "text.color": INK, "pdf.fonttype": 42})


def offsets(fai):
    off, lens, pos = {}, {}, 0
    for line in sorted(open(fai), key=lambda l: int(l.split()[0][3:])):
        c, n = line.split()[:2]
        off[c], lens[c] = pos, int(n)
        pos += int(n)
    return off, lens, pos


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("syri"); ap.add_argument("ref_fai"); ap.add_argument("qry_fai"); ap.add_argument("out")
    ap.add_argument("--y-name", required=True)
    ap.add_argument("--min-size", type=int, default=100000)
    ap.add_argument("--marker", choices=("circle", "box"), default="circle")
    ap.add_argument("--min-box", type=float, default=0.004, help="minimum box side as a fraction of the axis")
    ap.add_argument("--no-text", action="store_true", help="omit the title and explanatory notes (legend text goes in the figure legend); default off")
    a = ap.parse_args()
    ox, lx, tx = offsets(a.ref_fai)
    oy, ly, ty = offsets(a.qry_fai)

    segs, events, only_x, only_y = [], [], [], []
    for line in open(a.syri):
        f = line.rstrip("\n").split("\t")
        t = f[10]
        if t in ALIGN_ROWS:
            segs.append([(ox[f[0]] + int(f[1]), oy[f[5]] + int(f[6])), (ox[f[0]] + int(f[2]), oy[f[5]] + int(f[7]))])
        elif t in EVENT:
            if f[1] == "-" or f[6] == "-":
                continue
            r0, r1 = sorted((int(f[1]), int(f[2])))
            q0, q1 = sorted((int(f[6]), int(f[7])))
            size = max(r1 - r0, q1 - q0)
            if size >= a.min_size:
                events.append((ox[f[0]] + r0, ox[f[0]] + r1, oy[f[5]] + q0, oy[f[5]] + q1, t, size))
        elif t == "NOTAL":
            if f[1] != "-" and int(f[2]) - int(f[1]) >= a.min_size:
                only_x.append((ox[f[0]] + int(f[1]), int(f[2]) - int(f[1])))
            elif f[6] != "-" and int(f[7]) - int(f[6]) >= a.min_size:
                only_y.append((oy[f[5]] + int(f[6]), int(f[7]) - int(f[6])))

    fig, ax = plt.subplots(figsize=(9, 9.4))
    ax.add_collection(LineCollection(segs, colors=ALN, linewidths=1.8, zorder=1))
    ax.broken_barh(only_x, (-ty * 0.019, ty * 0.014), color=ABSENT, zorder=3)
    for s, w in only_y:
        ax.broken_barh([(-tx * 0.019, tx * 0.014)], (s, w), color=ABSENT, zorder=3)
    for x0, x1, y0, y1, t, size in events:
        if a.marker == "circle":
            ax.scatter((x0 + x1) / 2, (y0 + y1) / 2, s=18 + 90 * min(size / 2e6, 1), facecolors="none",
                       edgecolors=EVENT[t][0], linewidths=1.3, zorder=4)
        else:
            w, h = max(x1 - x0, tx * a.min_box), max(y1 - y0, ty * a.min_box)
            ax.add_patch(Rectangle(((x0 + x1 - w) / 2, (y0 + y1 - h) / 2), w, h, fill=False, ec=EVENT[t][0],
                                   lw=1.2, zorder=4))
    for c in ox:
        ax.axvline(ox[c], color=INK, lw=0.6)
    for c in oy:
        ax.axhline(oy[c], color=INK, lw=0.6)
    ax.set_xticks([ox[c] + lx[c] / 2 for c in ox]); ax.set_xticklabels(list(ox), rotation=90)
    ax.set_yticks([oy[c] + ly[c] / 2 for c in oy]); ax.set_yticklabels(list(oy))
    ax.set_xlim(-tx * 0.021, tx); ax.set_ylim(-ty * 0.021, ty)
    ax.set_xlabel("Thille chr1 to chr12" + ("" if a.no_text else f" ({tx / 1e6:.0f} Mb)"))
    ax.set_ylabel(f"{a.y_name} chr1 to chr12" + ("" if a.no_text else f", same strand as Thille ({ty / 1e6:.0f} Mb)"))
    kb = a.min_size // 1000
    how = "boxes span each event in both genomes" if a.marker == "box" else "circle size grows with event length"
    if not a.no_text:
        ax.set_title(f"Thille vs {a.y_name}: SyRI alignments with rearrangements and unaligned sequence >= {kb} kb\n"
                     f"{how}; red strips = sequence with no alignment in the other genome", fontsize=9)
    handles = [Line2D([], [], color=ALN, lw=1.5, label="alignment")]
    for col, lab in dict.fromkeys(EVENT.values()):
        handles.append(Line2D([], [], marker="s" if a.marker == "box" else "o", ls="", mfc="none", mec=col, mew=1.3,
                              label=lab))
    handles.append(Patch(color=ABSENT, label=f"not aligned (x: Thille only, y: {a.y_name} only)"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7, frameon=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    print(f"wrote {a.out}.pdf/.png: {len(segs)} alignments, {len(events)} rearrangements, "
          f"{len(only_x)} Thille-only and {len(only_y)} {a.y_name}-only regions >= {kb} kb")


if __name__ == "__main__":
    main()
