#!/usr/bin/env python3
"""Cluster indels across reads in a region and report the well-supported ones, optionally per haplotype.

Every CIGAR insertion/deletion >= --min-len bp in primary alignments is collected. Indels from different reads are
merged into one cluster when they are the same type, start within --pos-tol bp of the cluster and differ in length by
<= max(--len-tol bp, 10% of the length). Support = reads carrying the cluster / reads spanning it (aligned over
cluster start - 50 to end + 50). Random noisy-read indels scatter in position and length and stay below the support
cut; real indels stack. With --top N only the N best-supported clusters per locus (within --pos-tol) are kept; by default all are kept,
since random noisy-read indels are mostly 1 bp and do not stack.
With --groups (read<TAB>group from sdmyb_haplotypes.py) or --hap-sites (pos:alt,... read is group 1 when it carries
>= 60% of the listed alt alleles it covers, group 0 when <= 20%, needs >= 2 covered sites) support is also given per
group, so each indel can be placed on a haplotype.

Usage (from the repo root):
    python3 scripts/indel_clusters.py <bam> <chrom:start-end> <label> [--groups file | --hap-sites 10311829:T,...]
Writes <outdir>/<label>.indels.tsv (default outdir summaries/sdmyb_haplotypes)
"""
import argparse
import collections
import re
import subprocess
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("bam"); p.add_argument("region"); p.add_argument("label")
p.add_argument("--min-len", type=int, default=2); p.add_argument("--pos-tol", type=int, default=10)
p.add_argument("--len-tol", type=int, default=2); p.add_argument("--min-support", type=float, default=0.25)
p.add_argument("--min-reads", type=int, default=5)
p.add_argument("--top", type=int, default=0, help="keep only the top N clusters per locus (default 0 = keep all)")
p.add_argument("--groups"); p.add_argument("--hap-sites")
p.add_argument("--outdir", default="summaries/sdmyb_haplotypes")
p.add_argument("--group-min-support", type=float,
               help="also keep a cluster carried by >= this fraction of one group's spanning reads (needs --groups or "
                    "--hap-sites); a minority-haplotype indel can fall below --min-support over all reads")
a = p.parse_args()

out = subprocess.run(["singularity", "exec", "--no-home", "--bind", "/epool/edwin/pamericana/avoflow/mappingback:/epool/edwin/pamericana/avoflow/mappingback:ro", "containers/aligners_v3.sif", "samtools", "view",
                      "-F", "0x904", a.bam, a.region], check=True, text=True, capture_output=True).stdout
reads, events = {}, []
for line in out.splitlines():
    x = line.split("\t"); name, r = x[0], int(x[3]); bases = {}
    start = r; q = 0
    for n, op in re.findall(r"(\d+)([MIDNSHX=])", x[5]):
        n = int(n)
        if op in "M=X":
            for k in range(n):
                bases[r + k] = x[9][q + k]
            r += n; q += n
        elif op == "I":
            if n >= a.min_len: events.append((name, "INS", r - 1, n))
            q += n
        elif op in "DN":
            if n >= a.min_len: events.append((name, "DEL", r - 1, n))
            r += n
        elif op == "S":
            q += n
    reads[name] = (start, r - 1, bases)

group = {}
if a.groups:
    group = {l.split()[0]: int(l.split()[1]) for l in open(a.groups)}
elif a.hap_sites:
    alts = {int(s.split(":")[0]): s.split(":")[1] for s in a.hap_sites.split(",")}
    for n, (s, e, b) in reads.items():
        cov = [pos for pos in alts if pos in b and b[pos] != "-"]
        if len(cov) >= 2:
            f = sum(b[pos] == alts[pos] for pos in cov) / len(cov)
            if f >= 0.6: group[n] = 1
            elif f <= 0.2: group[n] = 0

clusters = []  # [type, pos, len, set(reads)]
for name, typ, pos, ln in sorted(events, key=lambda z: (z[1], z[2])):
    for c in clusters:
        if c[0] == typ and abs(c[1] - pos) <= a.pos_tol and abs(c[2] - ln) <= max(a.len_tol, 0.1 * c[2]):
            c[3].add(name); break
    else:
        clusters.append([typ, pos, ln, {name}])

rows = []
for typ, pos, ln, rs in clusters:
    end = pos + (ln if typ == "DEL" else 1)
    span = [n for n, (s, e, _) in reads.items() if s <= pos - 50 and e >= end + 50]
    if len(rs) < a.min_reads or not span:
        continue
    sup = len(rs & set(span)) / len(span)
    g, gf = {}, []
    for k in (0, 1):
        sp = [n for n in span if group.get(n) == k]
        g[k] = f"{len(rs & set(sp))}/{len(sp)}"
        if len(sp) >= 3:
            gf.append(len(rs & set(sp)) / len(sp))
    if sup < a.min_support and not (a.group_min_support and gf and max(gf) >= a.group_min_support):
        continue
    rows.append((pos, typ, ln, len(rs & set(span)), len(span), sup, g[0], g[1]))

# top two per locus
rows.sort()
kept, i = [], 0
while i < len(rows):
    j = i
    while j + 1 < len(rows) and rows[j + 1][0] - rows[i][0] <= a.pos_tol:
        j += 1
    loc = sorted(rows[i:j + 1], key=lambda r: -r[5])
    kept += loc[:a.top] if a.top else loc
    i = j + 1
Path(a.outdir).mkdir(parents=True, exist_ok=True)
with open(f"{a.outdir}/{a.label}.indels.tsv", "w") as fh:
    fh.write(f"# {a.bam} {a.region}; min len {a.min_len}, pos tol {a.pos_tol}, min support {a.min_support}; "
             f"reads grouped: {sum(1 for v in group.values() if v == 0)} group0, {sum(1 for v in group.values() if v == 1)} group1\n")
    fh.write("pos\ttype\tlen\treads_with\treads_spanning\tsupport\tgroup0_with/spanning\tgroup1_with/spanning\n")
    for r in kept:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]:.2f}\t{r[6]}\t{r[7]}\n")
print(f"{a.label}: {len(kept)} indel clusters (support >= {a.min_support}); groups {collections.Counter(group.values())}")
