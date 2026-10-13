#!/usr/bin/env python3
"""Per-read QV of a FASTQ (Edwin, 2026-10-02): QV = -10 log10(mean base error probability), the error probability of
each base being 10^(-Q/10) from its Phred score; the arithmetic mean of Phred scores is also kept for comparison.
Writes <out>.tsv.gz (read, length, QV, mean Phred), <out>.summary.tsv and a histogram <out>.{png,pdf}.
Usage: pigz -dc reads.fastq.gz | python3 scripts/read_qv_hist.py <out prefix> "<title>"
"""
import gzip, sys, math
import numpy as np
out, title = sys.argv[1], sys.argv[2]
LUT = 10 ** (-np.arange(128) / 10.0)
qv, mp, ln = [], [], []
with gzip.open(out + ".tsv.gz", "wt", compresslevel=3) as tsv:
    tsv.write("read\tlength\tqv\tmean_phred\n")
    it = iter(sys.stdin.buffer)
    for h in it:
        s = next(it); next(it); q = next(it).rstrip(b"\n")
        a = np.frombuffer(q, dtype=np.uint8).astype(np.int16) - 33
        e = LUT[a].mean(); v = -10 * math.log10(e) if e > 0 else 99.0; m = float(a.mean())
        qv.append(v); mp.append(m); ln.append(len(a))
        tsv.write(f"{h[1:].split()[0].decode()}\t{len(a)}\t{v:.2f}\t{m:.2f}\n")
qv, mp, ln = np.array(qv), np.array(mp), np.array(ln)
w = ln / ln.sum()
with open(out + ".summary.tsv", "w") as f:
    f.write("reads\tbases\tmedian_qv\tmean_qv\tbase_weighted_qv\tpct_reads_qv_ge_20\tpct_reads_qv_ge_30\tmedian_mean_phred\n")
    bw = -10 * math.log10((LUT[0] * 0 + (10 ** (-qv / 10)) * w).sum())
    f.write(f"{len(qv)}\t{ln.sum()}\t{np.median(qv):.2f}\t{qv.mean():.2f}\t{bw:.2f}\t{100*(qv>=20).mean():.2f}\t{100*(qv>=30).mean():.2f}\t{np.median(mp):.2f}\n")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.hist(qv, bins=np.arange(0, max(61, qv.max() + 1), 0.5), color="#3b6ea5")
ax.axvline(np.median(qv), color="black", lw=1, ls="--", label=f"median QV {np.median(qv):.1f}")
ax.set_xlabel("Per-read QV (from base error probabilities)"); ax.set_ylabel("Reads"); ax.set_title(title); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".pdf")
print(open(out + ".summary.tsv").read())
