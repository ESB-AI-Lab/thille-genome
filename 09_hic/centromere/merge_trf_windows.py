#!/usr/bin/env python3
"""Merge windowed TRF -ngs outputs back to chromosome coordinates.

Each window file is named <chr>_w<i>_<start>.fa.dat (written by sbatch_trf_windows.sh); its array coordinates are
shifted by <start>. Every listed file must have a .done marker. Arrays found in two overlapping windows are reported
once (identical shifted start, end and period). Arrays crossing a window edge can appear as partial pieces; for the
position tracks the caller merges overlapping pieces, so positions are unaffected.

Usage: merge_trf_windows.py <chr> <out.dat> <window.fa.dat> [<window.fa.dat> ...]
"""
import os
import re
import sys

if len(sys.argv) < 4:
    sys.exit(__doc__)
chrom, out, files = sys.argv[1], sys.argv[2], sys.argv[3:]
missing = [f for f in files if not os.path.exists(f + ".done")]
if missing:
    sys.exit(f"ERROR: unfinished windows: {missing}")
seen, rows, spans = set(), [], []
for f in files:
    m = re.search(r"_w\d+_(\d+)\.fa\.dat$", f)
    if not m:
        sys.exit(f"ERROR: cannot read window start from {f}")
    off = int(m.group(1))
    spans.append(off)
    for line in open(f):
        if line.startswith("@"):
            continue
        p = line.split()
        s, e = int(p[0]) + off, int(p[1]) + off
        key = (s, e, p[2])
        if key in seen:
            continue
        seen.add(key)
        rows.append((s, e, p[2:]))
rows.sort()
with open(out + ".tmp", "w") as fh:
    fh.write(f"@{chrom}\n")
    for s, e, rest in rows:
        fh.write(f"{s} {e} {' '.join(rest)}\n")
os.replace(out + ".tmp", out)
with open(out + ".done", "w") as fh:
    fh.write("merged from windows starting at: " + ",".join(map(str, sorted(spans))) + "\n")
print(f"{chrom}: {len(rows)} arrays from {len(files)} windows -> {out}")
