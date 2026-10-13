#!/usr/bin/env python3
"""Same-chromosome pairing of high-identity ModDotPlot windows and the identity of off-chromosome windows.
Usage: moddotplot_stats.py <COMPARE.bedpe> <x.boundaries.tsv> <y.boundaries.tsv> [--hi 96] [--lo 93]"""
import argparse, csv, statistics
ap = argparse.ArgumentParser(); ap.add_argument("bedpe"); ap.add_argument("xb"); ap.add_argument("yb")
ap.add_argument("--hi", type=float, default=96); ap.add_argument("--lo", type=float, default=93); a = ap.parse_args()
def bounds(p):
    return [(int(r["chrom_number"]), int(r["start_0based"]), int(r["end"])) for r in csv.DictReader(open(p), delimiter="\t")]
def which(b, p):
    for n, s, e in b:
        if s <= p < e:
            return n
xb, yb = bounds(a.xb), bounds(a.yb)
hi = same = 0; off = []
for l in open(a.bedpe):
    if l.startswith("#"):
        continue
    f = l.split("\t"); x = (int(f[1]) + int(f[2])) // 2; y = (int(f[4]) + int(f[5])) // 2; idn = float(f[6])
    cx, cy = which(xb, x), which(yb, y)
    if cx is None or cy is None:
        continue
    if idn >= a.hi:
        hi += 1; same += cx == cy
    if cx != cy:
        off.append(idn)
lt = sum(1 for v in off if v < a.lo)
print(f">= {a.hi}%: {same} of {hi} ({100 * same / hi:.1f}%) same chromosome; different-chromosome windows: "
      f"{100 * lt / len(off):.1f}% below {a.lo}% (median {statistics.median(off):.1f}%; n={len(off)})")
