#!/usr/bin/env python3
"""Find clipped-read breakpoints in a region of one BAM and count the reads that break vs span each one.

A breakpoint is a reference position where >= --min-reads alignments (primary or supplementary) end with a right
soft/hard clip >= --min-clip bp, or start with a left clip >= --min-clip bp, within +-5 bp. For each breakpoint:
reads breaking there, primary reads spanning it (aligned >= 200 bp on both sides, no clip at it), and where the
clipped parts align (SA tags: chromosome:position, count). Samtools 1.21 in containers/aligners_v3.sif.

Usage: python3 scripts/breakpoint_scan.py <bam> <chrom:start-end> [--min-clip 100] [--min-reads 3]
"""
import argparse
import collections
import re
import subprocess

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("bam"); p.add_argument("region")
p.add_argument("--min-clip", type=int, default=100); p.add_argument("--min-reads", type=int, default=3)
a = p.parse_args()
chrom, s, e = re.fullmatch(r"(\w+):(\d+)-(\d+)", a.region).groups()
out = subprocess.run(["singularity", "exec", "--bind", "/epool", "containers/aligners_v3.sif", "samtools", "view",
                      "-F", "0x104", a.bam, a.region], check=True, text=True, capture_output=True).stdout
alns, events = [], collections.defaultdict(list)
for l in out.splitlines():
    x = l.split("\t"); flag = int(x[1]); pos = int(x[3])
    ops = re.findall(r"(\d+)([MIDNSHX=])", x[5])
    end = pos + sum(int(n) for n, o in ops if o in "MDN=X") - 1
    lc = int(ops[0][0]) if ops[0][1] in "SH" else 0
    rc = int(ops[-1][0]) if ops[-1][1] in "SH" else 0
    sa = re.search(r"SA:Z:(\S+)", l)
    sat = [f"{t.split(',')[0]}:{int(t.split(',')[1]) // 1000}kb" for t in sa.group(1).split(";") if t] if sa else []
    alns.append((pos, end, lc, rc, flag & 2048))
    if rc >= a.min_clip: events[end].append(("R", x[0], sat))
    if lc >= a.min_clip: events[pos].append(("L", x[0], sat))
# merge event positions within 5 bp
keys = sorted(events); groups = []
for k in keys:
    if groups and k - groups[-1][-1] <= 5: groups[-1].append(k)
    else: groups.append([k])
print("position\tbreaking_reads\tside\tspanning_primary\tclipped_parts_align_to")
for g in groups:
    ev = [v for k in g for v in events[k]]
    if len({r for _, r, _ in ev}) < a.min_reads: continue
    bp = g[len(g) // 2]
    span = sum(1 for p0, e0, lc, rc, sup in alns if not sup and p0 <= bp - 200 and e0 >= bp + 200)
    sides = collections.Counter(side for side, _, _ in ev)
    targets = collections.Counter(t for _, _, ts in ev for t in ts if not t.startswith(f"{chrom}:{bp // 1000}"))
    print(f"{chrom}:{bp}\t{len({r for _, r, _ in ev})}\t{dict(sides)}\t{span}\t" + ", ".join(f"{t} x{c}" for t, c in targets.most_common(4)))
