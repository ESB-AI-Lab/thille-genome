#!/usr/bin/env python3
"""Thille sequence >= 100 kb with no counterpart in a compared genome (SyRI NOTAL with a Thille interval, DEL, CPL),
by chromosome and relative to the Thille centromere calls. Usage: thille_only_stats.py <syri.out> <centromere.bed>"""
import collections, sys
syri, cenf = sys.argv[1:]
cen = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in open(cenf)}
ev = []
for l in open(syri):
    f = l.rstrip("\n").split("\t"); t = f[10]
    if f[1] == "-":
        continue
    s, e = int(f[1]) - 1, int(f[2])
    if t == "NOTAL" and f[5] != "-":
        continue
    if t in ("NOTAL", "DEL", "CPL") and e - s >= 100000:
        ev.append((f[0], s, e, t))
cnt = collections.Counter(t for *_, t in ev)
tot = sum(e - s for _, s, e, _ in ev)
by = collections.Counter(); near = 0
for c, s, e, t in ev:
    by[c] += e - s
    cs, ce = cen[c]
    lo, hi = cs - 5_000_000, ce + 5_000_000
    near += max(0, min(e, hi) - max(s, lo))
print(f"{dict(cnt)}; total {tot / 1e6:.1f} Mb; within 5 Mb of a centromere call {near / 1e6:.1f} Mb ({100 * near / tot:.0f}%)")
print("by chromosome:", ", ".join(f"{c} {v / 1e6:.1f}" for c, v in by.most_common()))
for c in sorted(by, key=lambda x: int(x[3:])):
    big = max((x for x in ev if x[0] == c), key=lambda x: x[2] - x[1])
    cs, ce = cen[c]
    ov = big[1] < ce and big[2] > cs
    print(f"  largest {c}: {(big[2] - big[1]) / 1e6:.2f} Mb {big[1] / 1e6:.2f}-{big[2] / 1e6:.2f} {'overlaps centromere' if ov else ''}")
c7 = [x for x in ev if x[0] == "chr7" and x[1] >= 31e6 and x[2] <= 51e6]
print(f"chr7 31-51 Mb: {len(c7)} regions, {sum(e - s for _, s, e, _ in c7) / 1e6:.1f} Mb")
