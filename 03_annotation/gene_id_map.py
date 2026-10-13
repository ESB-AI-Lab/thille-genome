#!/usr/bin/env python3
"""Map Thille Tiberius gene IDs from v3 to v4 (Edwin, 2026-10-04).
Genes are matched by chromosome, strand and the full CDS structure of the longest transcript. v3 chr8 coordinates are
lifted to v4 by +25,531,594 (v4 chr8 = scaffold_13 1-25,531,494 rc + 100 bp gap + v3 chr8, from the curation AGP);
chr8 genes are matched after the lift, and any v3 chr8 gene without an identical v4 model is listed as 'changed' with the
v4 gene(s) overlapping its lifted span.
Usage: python3 scripts/gene_id_map.py > summaries/gene_id_map_v3_v4.tsv"""
import re, sys, collections
V3 = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
V4 = "annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf"
SHIFT = {"chr8": 25531594}

def load(path, shift=None):
    genes, cds, tx_gene = {}, collections.defaultdict(list), {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        s = (shift or {}).get(f[0], 0)
        if f[2] == "gene":
            gid = re.search(r'gene_id "([^"]+)"', f[8]) or re.match(r"(\S+)", f[8])
            genes[gid.group(1)] = (f[0], int(f[3]) + s, int(f[4]) + s, f[6])
        elif f[2] == "CDS":
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx_gene[tid] = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            cds[tid].append((int(f[3]) + s, int(f[4]) + s))
    key = {}
    for tid, blocks in cds.items():
        g = tx_gene[tid]; c, a, b, st = genes[g]
        key.setdefault(g, (c, st, tuple(sorted(blocks))))
    return genes, key

g3, k3 = load(V3, SHIFT); g4, k4 = load(V4)
by4 = collections.defaultdict(list)
for g, k in k4.items(): by4[k].append(g)
by_chr4 = collections.defaultdict(list)
for g, (c, a, b, st) in g4.items(): by_chr4[c].append((a, b, g))
print("v3_gene\tv4_gene\tchrom\tv3_start_v4coords\tv3_end_v4coords\tstrand\tstatus")
n = collections.Counter()
for g, (c, a, b, st) in sorted(g3.items(), key=lambda x: int(x[0][1:])):
    hits = by4.get(k3.get(g), [])
    if len(hits) == 1:
        st_ = "identical"; v4 = hits[0]
    elif len(hits) > 1:
        st_ = "ambiguous"; v4 = ",".join(hits)
    else:
        ov = [x[2] for x in by_chr4[c] if x[0] <= b and x[1] >= a]
        st_ = "changed"; v4 = ",".join(ov) if ov else "none"
    n[(c == "chr8", st_)] += 1
    print(f"{g}\t{v4}\t{c}\t{a}\t{b}\t{st}\t{st_}")
for k, v in sorted(n.items()): print(f"# chr8={k[0]} {k[1]}: {v}", file=sys.stderr)
