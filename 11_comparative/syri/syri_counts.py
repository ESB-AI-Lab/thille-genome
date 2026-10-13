#!/usr/bin/env python3
"""Counts quoted in the Results SyRI paragraph, from <name>syri.out: syntenic Thille Mb (and % of the Thille chromosome
length), and events >= --min-size (size = larger of the two spans) by SyRI type: inversions (INV), translocations
(TRANS + INVTR, inverted = INVTR), duplications (DUP + INVDP, inverted = INVDP), deletions (DEL), insertions (INS),
copy gains (CPG) and losses (CPL), and NOTAL regions on each side (count, Mb).
Usage: syri_counts.py <syri.out> <thille_total_bp> [--min-size 100000]
"""
import argparse
import collections

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("syri"); ap.add_argument("total", type=int); ap.add_argument("--min-size", type=int, default=100000)
a = ap.parse_args()
syn, n, notal = 0, collections.Counter(), collections.defaultdict(lambda: [0, 0])
span = lambda s, e: (int(e) - int(s) + 1) if s not in ("-", "") and e not in ("-", "") else 0  # noqa: E731
for line in open(a.syri):
    f = line.rstrip("\n").split("\t")
    t = f[10]
    rs, qs = span(f[1], f[2]), span(f[6], f[7])
    if t == "SYN":
        syn += rs
        continue
    size = max(rs, qs)
    if t == "NOTAL":
        side = "thille" if rs else "other"
        if size >= a.min_size:
            notal[side][0] += 1; notal[side][1] += size
    elif t in ("INV", "TRANS", "INVTR", "DUP", "INVDP", "DEL", "INS", "CPG", "CPL") and size >= a.min_size:
        n[t] += 1
print(f"syntenic Thille {syn / 1e6:.1f} Mb ({100 * syn / a.total:.1f}%)")
print(f"inversions {n['INV']}; translocations {n['TRANS'] + n['INVTR']} ({n['INVTR']} inverted); duplications "
      f"{n['DUP'] + n['INVDP']} ({n['INVDP']} inverted); deletions {n['DEL']}; insertions {n['INS']}; copy gains "
      f"{n['CPG']}; copy losses {n['CPL']}")
print(f"NOTAL >= {a.min_size // 1000} kb: Thille-only {notal['thille'][0]} ({notal['thille'][1] / 1e6:.1f} Mb), "
      f"other-only {notal['other'][0]} ({notal['other'][1] / 1e6:.1f} Mb)")
