#!/usr/bin/env python3
"""Per-read alignment identity over primary alignments of a minimap2 BAM (Edwin, 2026-10-03).
BLAST identity = matching bases / alignment columns, columns = M + I + D bases, mismatches = NM - I - D
(standard definition: every inserted or deleted base is one column). Also gap-compressed identity = 1 - de.
Reports mean, median and pooled (sum matches / sum columns) BLAST identity, and mean/median gap-compressed identity.
Usage: samtools view -F 0x904 <bam> | python3 scripts/identity_stats.py <label>"""
import sys, re, statistics
OPS = re.compile(r"(\d+)([MIDNSHP=X])")
lab = sys.argv[1]; ids = []; gcs = []; mt = cl = 0
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if int(f[1]) & 0x904 or f[5] == "*": continue
    m = i = d = 0
    for n, o in OPS.findall(f[5]):
        n = int(n)
        if o in "M=X": m += n
        elif o == "I": i += n
        elif o == "D": d += n
    tags = "\t".join(f[11:])
    nm = int(re.search(r"NM:i:(\d+)", tags).group(1)); x = nm - i - d
    match = m - x; cols = m + i + d
    ids.append(match / cols); mt += match; cl += cols
    g = re.search(r"de:f:([0-9.eE-]+)", tags)
    if g: gcs.append(1 - float(g.group(1)))
print("label\tprimary_alignments\tmean_blast_identity_pct\tmedian_blast_identity_pct\tpooled_blast_identity_pct\tmean_gapcomp_identity_pct\tmedian_gapcomp_identity_pct")
print(f"{lab}\t{len(ids)}\t{100*statistics.fmean(ids):.3f}\t{100*statistics.median(ids):.3f}\t{100*mt/cl:.3f}\t{100*statistics.fmean(gcs):.3f}\t{100*statistics.median(gcs):.3f}")
