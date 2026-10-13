#!/usr/bin/env python3
"""Split a scaffolded FASTA into contigs at every run of N (gaps), keeping each contig's origin.

Contig names: <scaffold>_<start>_<end> (0-based start, end exclusive, in scaffold coordinates). Writes the contigs
FASTA and a TSV of origins (contig, scaffold, start, end, gap_before_bp), then checks that rejoining the contigs with
their original gaps reproduces every scaffold exactly.
Usage: split_at_gaps.py <in.fasta> <out_contigs.fasta> <out_origins.tsv> [--min-gap 1]
"""
import argparse
import re


def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("fasta"); ap.add_argument("out_fa"); ap.add_argument("out_tsv")
ap.add_argument("--min-gap", type=int, default=1)
a = ap.parse_args()
gap = re.compile(r"[Nn]{%d,}" % a.min_gap)
nseq = ncontig = 0
with open(a.out_fa, "w") as fa, open(a.out_tsv, "w") as tsv:
    tsv.write("contig\tscaffold\tstart\tend\tgap_before_bp\n")
    for name, seq in read_fasta(a.fasta):
        nseq += 1
        pieces, pos = [], 0
        for m in gap.finditer(seq):
            if m.start() > pos:
                pieces.append((pos, m.start()))
            pos = m.end()
        if pos < len(seq):
            pieces.append((pos, len(seq)))
        rebuilt, last = [], 0
        for s, e in pieces:
            cid = f"{name}_{s}_{e}"
            fa.write(f">{cid}\n")
            for i in range(s, e, 60):
                fa.write(seq[i:min(i + 60, e)] + "\n")
            tsv.write(f"{cid}\t{name}\t{s}\t{e}\t{s - last if rebuilt else s}\n")
            rebuilt.append(seq[last:s] if rebuilt else seq[:s])
            rebuilt.append(seq[s:e])
            last = e
            ncontig += 1
        rebuilt.append(seq[last:])
        assert "".join(rebuilt) == seq, f"rebuild failed for {name}"
print(f"{nseq} scaffolds -> {ncontig} contigs; every scaffold rebuilds exactly from its contigs and gaps")
