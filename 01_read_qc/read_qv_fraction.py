#!/usr/bin/env python3
"""Fraction of reads at or above Q20 and Q30 (Edwin, 2026-10-03), two per-read definitions:
QV from the mean base error probability, -10 log10(mean 10^(-Q/10)), as for the HiFi histogram, and the arithmetic
mean Phred score. Reads are processed in batches of equal length for speed.
Usage: pigz -dc reads.fastq.gz | python3 scripts/read_qv_fraction.py <label>"""
import sys, numpy as np
lab = sys.argv[1]; LUT = 10 ** (-np.arange(128) / 10.0)
n = q20 = q30 = m20 = m30 = 0; batch = []
def flush(b):
    global n, q20, q30, m20, m30
    by = {}
    for q in b: by.setdefault(len(q), []).append(q)
    for L, qs in by.items():
        a = np.frombuffer(b"".join(qs), dtype=np.uint8).reshape(-1, L).astype(np.int16) - 33
        qv = -10 * np.log10(LUT[a].mean(1)); mp = a.mean(1)
        n += len(qs); q20 += int((qv >= 20).sum()); q30 += int((qv >= 30).sum()); m20 += int((mp >= 20).sum()); m30 += int((mp >= 30).sum())
it = iter(sys.stdin.buffer)
for h in it:
    next(it); next(it); batch.append(next(it).rstrip(b"\n"))
    if len(batch) == 1_000_000: flush(batch); batch = []
flush(batch)
print("label\treads\tpct_qv_ge20\tpct_qv_ge30\tpct_meanphred_ge20\tpct_meanphred_ge30")
print(f"{lab}\t{n}\t{100*q20/n:.2f}\t{100*q30/n:.2f}\t{100*m20/n:.2f}\t{100*m30/n:.2f}")
