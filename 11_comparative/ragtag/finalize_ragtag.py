#!/usr/bin/env python3
"""Final FASTA/AGP for a RagTag scaffold run on Thille v4: chrN_RagTag -> chrN; unplaced contigs keep their own name
(the _RagTag suffix removed). Checks that every input contig is recovered exactly from the output through the AGP (reverse-complemented where
placed on the minus strand), each used once in full. Writes <prefix>.fasta (chr1..chr12 first, then unplaced by
length), <prefix>.agp and <prefix>.stats.tsv (per chromosome: placed contigs, bp, joins).
Usage: finalize_ragtag.py <ragtag_scaffold_dir> <input_contigs.fasta> <prefix>
"""
import collections
import re
import sys


def read(p):
    n, b = None, []
    for l in open(p):
        if l[0] == ">":
            if n is not None:
                yield n, "".join(b)
            n, b = l[1:].split()[0], []
        else:
            b.append(l.strip())
    if n is not None:
        yield n, "".join(b)


d, inp, pre = sys.argv[1:]
rename = lambda n: re.sub(r"_RagTag$", "", n)  # noqa: E731
seqs = {rename(n): s for n, s in read(f"{d}/ragtag.scaffold.fasta")}
# Exact check: every AGP component, read from the output (reverse-complemented when placed on the minus strand),
# must equal its input contig, and every input contig must be used exactly once in full.
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
contigs = dict(read(inp))
used = collections.Counter()
for l in open(f"{d}/ragtag.scaffold.agp"):
    if l.startswith("#"):
        continue
    f = l.rstrip("\n").split("\t")
    if f[4] != "W":
        continue
    seg = seqs[rename(f[0])][int(f[1]) - 1:int(f[2])]
    if f[8] == "-":
        seg = seg.translate(RC)[::-1]
    src = rename(f[5])
    if seg != contigs[src][int(f[6]) - 1:int(f[7])]:
        sys.exit(f"component {src} {f[6]}-{f[7]} does not match its input contig")
    used[src] += int(f[7]) - int(f[6]) + 1
bad = [c for c in contigs if used[c] != len(contigs[c])]
if bad:
    sys.exit(f"{len(bad)} input contigs not used exactly once in full, e.g. {bad[:3]}")
cin = sum(len(s) for s in contigs.values())
cin = collections.Counter({"bases": cin})
chroms = sorted((n for n in seqs if re.fullmatch(r"chr\d+", n)), key=lambda n: int(n[3:]))
rest = sorted((n for n in seqs if n not in chroms), key=lambda n: -len(seqs[n]))
with open(f"{pre}.fasta", "w") as fh:
    for n in chroms + rest:
        fh.write(f">{n}\n")
        s = seqs[n]
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")
st = collections.defaultdict(lambda: [0, 0, 0])
with open(f"{pre}.agp", "w") as fh:
    for l in open(f"{d}/ragtag.scaffold.agp"):
        if l.startswith("#"):
            fh.write(l)
            continue
        f = l.rstrip("\n").split("\t")
        f[0] = rename(f[0])
        if f[4] == "W":
            f[5] = rename(f[5])
            if f[0] in chroms:
                st[f[0]][0] += 1
                st[f[0]][1] += int(f[2]) - int(f[1]) + 1
        elif f[0] in chroms:
            st[f[0]][2] += 1
        fh.write("\t".join(f) + "\n")
with open(f"{pre}.stats.tsv", "w") as fh:
    fh.write("chrom\tlength\tcontigs\tcontig_bp\tjoins\n")
    for n in chroms:
        fh.write(f"{n}\t{len(seqs[n])}\t{st[n][0]}\t{st[n][1]}\t{st[n][2]}\n")
print(f"{pre}: {len(chroms)} chromosomes ({sum(len(seqs[n]) for n in chroms):,} bp), {len(rest)} unplaced "
      f"({sum(len(seqs[n]) for n in rest):,} bp); all {len(contigs):,} input contigs ({cin['bases']:,} bp) recovered exactly")
