#!/usr/bin/env python3
"""Gene and TE density at chromosome ends and around centromere calls (manuscript Gene Annotation paragraph).
Genes counted by midpoint. Ends = terminal 10% of each chromosome at each end; interior = the rest. Near-centromere =
within W bp of the centromere region (region extended by W on each side); TE coverage = union of TE intervals.
Usage: density_stats.py <genes.gtf> <TE.gff> <centromere.bed> <fai> [--w 5000000]
"""
import argparse
import bisect
import collections

ap = argparse.ArgumentParser(); ap.add_argument("gtf"); ap.add_argument("te"); ap.add_argument("cen"); ap.add_argument("fai")
ap.add_argument("--w", type=int, default=5_000_000); a = ap.parse_args()
L = {l.split()[0]: int(l.split()[1]) for l in open(a.fai)}
mids = collections.defaultdict(list)
for l in open(a.gtf):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "gene":
        mids[f[0]].append((int(f[3]) - 1 + int(f[4])) / 2)
cen = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in open(a.cen)}
te = collections.defaultdict(list)
for l in open(a.te):
    if l.startswith("#"):
        continue
    f = l.split("\t")
    if f[0] in L:
        te[f[0]].append((int(f[3]) - 1, int(f[4])))


def union(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered(u, s, e):
    return sum(max(0, min(e, b) - max(s, a)) for a, b in u)


endn = endbp = inn = inbp = 0; denser = 0
nearn = nearbp = nearte = 0; farte = 0; tot_bp = 0
for c in sorted(L, key=lambda x: int(x[3:])):
    Lc = L[c]; m = sorted(mids[c]); t = Lc * 0.1
    e_n = sum(1 for x in m if x < t or x >= Lc - t); i_n = len(m) - e_n
    endn += e_n; endbp += 2 * t; inn += i_n; inbp += Lc - 2 * t
    denser += (e_n / (2 * t)) > (i_n / (Lc - 2 * t))
    s, e = max(0, cen[c][0] - a.w), min(Lc, cen[c][1] + a.w)
    nearn += sum(1 for x in m if s <= x < e); nearbp += e - s
    u = union(te[c]); nt = covered(u, s, e); nearte += nt; farte += covered(u, 0, Lc) - nt
    tot_bp += Lc
ngen = sum(len(v) for v in mids.values())
print(f"ends {endn / endbp * 1e6:.1f} genes/Mb, interiors {inn / inbp * 1e6:.1f}; ends denser on {denser} of {len(L)}")
print(f"within {a.w / 1e6:.0f} Mb of centromeres {nearn / nearbp * 1e6:.1f} genes/Mb (genome-wide {ngen / tot_bp * 1e6:.1f}); "
      f"TE cover {100 * nearte / nearbp:.1f}% vs elsewhere {100 * farte / (tot_bp - nearbp):.1f}%")
