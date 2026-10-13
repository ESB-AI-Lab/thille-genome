#!/usr/bin/env python3
"""Summarize inter-chromosomal (trans) Hi-C contacts in the chr1-12 matrix (read-only copy in hic/chr12/).

Per chromosome pair: observed trans pairs / expected from the product of each chromosome's trans totals.
Per bin pair: ICE-balanced frequency / median balanced trans bin (cooler weights stored in the mcool).
Reads the mcool with h5py (cooler is not installed); raw counts and weights are used as stored.

Usage (from the repo root):
    python3 scripts/hic_trans_contacts.py [--mcool hic/chr12/chr12.mcool] [--resolution 500000] [--top 25]
Writes summaries/hic_trans_chrom_pairs.tsv and summaries/hic_trans_top_bins.tsv.
"""
import argparse
import h5py
import numpy as np

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--mcool", default="hic/chr12/chr12.mcool")
p.add_argument("--resolution", default="500000")
p.add_argument("--top", type=int, default=25)
a = p.parse_args()

g = h5py.File(a.mcool, "r")[f"resolutions/{a.resolution}"]
ch, st, w = g["bins/chrom"][:], g["bins/start"][:], g["bins/weight"][:]
names = [n.decode() for n in g["chroms/name"][:]]
b1, b2, c = g["pixels/bin1_id"][:], g["pixels/bin2_id"][:], g["pixels/count"][:].astype(float)
trans = ch[b1] != ch[b2]
n = len(names)
M = np.zeros((n, n))
np.add.at(M, (ch[b1[trans]], ch[b2[trans]]), c[trans])
M = M + M.T
row = M.sum(1)
E = np.outer(row, row) / row.sum()
with open("summaries/hic_trans_chrom_pairs.tsv", "w") as out:
    out.write(f"# {a.mcool} {a.resolution} bp; total pairs {c.sum():.0f}; trans fraction {c[trans].sum() / c.sum():.4f}\n")
    out.write("chrom1\tchrom2\ttrans_pairs\tobs_over_exp\n")
    for i in range(n):
        for j in range(i + 1, n):
            out.write(f"{names[i]}\t{names[j]}\t{M[i, j]:.0f}\t{M[i, j] / E[i, j]:.3f}\n")
bal = c * w[b1] * w[b2]
t = trans & np.isfinite(bal)
med = np.median(bal[t])
with open("summaries/hic_trans_top_bins.tsv", "w") as out:
    out.write(f"# balanced median trans bin {med:.6g}\n")
    out.write("chrom1\tstart1\tchrom2\tstart2\traw_count\tbalanced\tbalanced_over_median\n")
    for k in np.flatnonzero(t)[np.argsort(-bal[t])][: a.top]:
        out.write(f"{names[ch[b1[k]]]}\t{st[b1[k]]}\t{names[ch[b2[k]]]}\t{st[b2[k]]}\t{c[k]:.0f}\t{bal[k]:.6g}\t{bal[k] / med:.1f}\n")
