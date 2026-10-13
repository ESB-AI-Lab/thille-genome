#!/usr/bin/env python3
"""Large (>= --min-size) presence/absence and copy-number candidates between Gwen contigs and Thille chromosomes.

Gwen is a contig assembly, so SyRI does not apply. From the Gwen -> Thille PAF (Thille = target), primary alignment
blocks of at least --min-block with MAPQ >= --min-mapq are projected onto Thille:
  uncovered : Thille intervals >= --min-size that no Gwen contig covers (Gwen deletion OR missing from the Gwen assembly)
  multi     : Thille intervals >= --min-size covered by >= 2 different Gwen contigs (Gwen duplication OR retained haplotig)
Each interval is then tested against Gwen read depth on Thille (samtools bedcov over the Gwen CLR BAM), relative to
the median depth of 1 Mb windows on that chromosome, using the coverage classes of the Methods: <= 0.25 absent,
0.25-0.75 half (one haplotype), 0.75-1.25 two copies, > 1.25 more than two copies.

Usage: gwen_sv_candidates.py <gwen_vs_thille.paf> <thille.fai> <gwen_reads_on_thille.bam> <out_prefix>
                             [--min-size 100000] [--min-block 10000] [--min-mapq 20] [--samtools PATH]
Writes <out_prefix>.tsv (one row per interval) and prints a summary.
"""
import argparse
import collections
import statistics
import subprocess
import sys


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def bedcov(samtools, bam, rows):
    """Mean depth per (chrom, start, end) via samtools bedcov (sum of per-base depth / length)."""
    bed = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in rows)
    r = subprocess.run([samtools, "bedcov", "/dev/stdin", bam], input=bed, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"samtools bedcov failed ({r.returncode}): {r.stderr.strip()}")
    return [int(l.split("\t")[3]) / (int(l.split("\t")[2]) - int(l.split("\t")[1])) for l in r.stdout.splitlines()]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf"); ap.add_argument("fai"); ap.add_argument("bam"); ap.add_argument("out")
    ap.add_argument("--min-size", type=int, default=100000)
    ap.add_argument("--min-block", type=int, default=10000)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--samtools", default="/gpool/bin/samtools-1.24/samtools")
    a = ap.parse_args()

    lens = {l.split()[0]: int(l.split()[1]) for l in open(a.fai)}
    cov = collections.defaultdict(list)                       # chrom -> [(s, e, contig)]
    for line in open(a.paf):
        f = line.split("\t")
        if int(f[10]) < a.min_block or int(f[11]) < a.min_mapq or "tp:A:P" not in line:
            continue
        cov[f[5]].append((int(f[7]), int(f[8]), f[0]))

    rows = []
    for c, L in lens.items():
        blocks = merge([(s, e) for s, e, _ in cov[c]])
        prev = 0
        for s, e in blocks + [[L, L]]:
            if s - prev >= a.min_size:
                rows.append([c, prev, s, "uncovered", 0, ""])
            prev = max(prev, e)
        # depth of distinct contigs along the chromosome (sweep over merged per-contig intervals)
        ev = []
        per = collections.defaultdict(list)
        for s, e, q in cov[c]:
            per[q].append((s, e))
        for q, iv in per.items():
            for s, e in merge(iv):
                ev += [(s, 1, q), (e, -1, q)]
        ev.sort(key=lambda x: (x[0], x[1]))
        active, start = collections.Counter(), None
        for pos, d, q in ev:
            was = sum(1 for v in active.values() if v > 0)
            active[q] += d
            now = sum(1 for v in active.values() if v > 0)
            if was < 2 <= now:
                start, names = pos, set()
            if now >= 2:
                names = names | {k for k, v in active.items() if v > 0}
            if was >= 2 > now and pos - start >= a.min_size:
                rows.append([c, start, pos, "multi", len(names), ",".join(sorted(names))])

    # chromosome median depth from 1 Mb windows
    win = [(c, s, min(s + 1000000, L)) for c, L in lens.items() for s in range(0, L, 1000000)]
    # one pass over the BAM for windows and candidate intervals together
    allcov = bedcov(a.samtools, a.bam, win + [(r[0], r[1], r[2]) for r in rows])
    wd, depths = allcov[:len(win)], allcov[len(win):]
    med = {c: statistics.median([d for (cc, _, _), d in zip(win, wd) if cc == c]) for c in lens}

    def cls(r):
        return "absent" if r <= 0.25 else "half" if r <= 0.75 else "two_copies" if r <= 1.25 else "gt_two_copies"

    order = {c: i for i, c in enumerate(lens)}
    rows.sort(key=lambda r: (order[r[0]], r[1]))
    with open(a.out + ".tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tsize_bp\tkind\tn_contigs\tcontigs\tgwen_read_depth\tchrom_median_depth\t"
                 "depth_ratio\tread_class\tat_chrom_end\n")
        for r, d in zip(rows, depths):
            ratio = d / med[r[0]] if med[r[0]] else float("nan")
            end = "yes" if r[1] == 0 or r[2] == lens[r[0]] else "no"
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[2]-r[1]}\t{r[3]}\t{r[4]}\t{r[5]}\t{d:.1f}\t{med[r[0]]:.1f}\t"
                     f"{ratio:.2f}\t{cls(ratio)}\t{end}\n")
    n = collections.Counter((r[3], cls(d / med[r[0]])) for r, d in zip(rows, depths))
    print(f"{len(rows)} intervals >= {a.min_size} bp -> {a.out}.tsv", file=sys.stderr)
    for k, v in sorted(n.items()):
        print(f"  {k[0]:9s} {k[1]:14s} {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
