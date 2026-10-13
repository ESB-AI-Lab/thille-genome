#!/usr/bin/env python3
"""Assign query contigs to reference chromosomes from a minimap2 PAF and build one
pseudo-chromosome per reference chromosome for per-chromosome dot plots.

Rule (user-approved):
  - keep PAF lines with MAPQ >= --min-mapq
  - contig -> reference chr with the most matching bases (PAF col 10)
  - contigs ordered along that chr by the matching-base-weighted midpoint of their alignments
  - contig reverse-complemented if minus-strand matching bases on that chr exceed plus-strand
  - concatenated without spacers; contigs with no qualifying alignment are listed as unassigned

Usage:
  build_gwen_pseudochr.py <aln.paf> <query.fasta> <out_prefix> <name_prefix> [--min-mapq 20]
  build_gwen_pseudochr.py compare/gwen/gwen_vs_thille.asm5.paf gwen.fasta compare/gwen/asm5/gwen_asm5 gwen_asm5
Outputs:
  <out_prefix>.<chr>.fasta (one record named <name_prefix>_<chr>), <out_prefix>.assign.tsv,
  <out_prefix>.unassigned.tsv, <out_prefix>.boundaries.tsv, <out_prefix>.summary.tsv
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fastaio import read_fasta, revcomp, write_record, die  # noqa: E402


def parse_paf(path, min_mapq):
    """Return per-contig, per-chr stats: matches, plus/minus matches, weighted midpoint sum."""
    stats = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0.0]))  # [match, plus, minus, mid*match]
    n_lines = n_kept = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                die(f"{path}: malformed PAF line {n_lines + 1}")
            n_lines += 1
            if int(f[11]) < min_mapq:
                continue
            n_kept += 1
            q, strand, t = f[0], f[4], f[5]
            tstart, tend, match = int(f[7]), int(f[8]), int(f[9])
            s = stats[q][t]
            s[0] += match
            s[1 if strand == "+" else 2] += match
            s[3] += (tstart + tend) / 2 * match
    print(f"PAF: {n_lines} lines, {n_kept} with MAPQ >= {min_mapq}", file=sys.stderr)
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf")
    ap.add_argument("query_fasta")
    ap.add_argument("out_prefix")
    ap.add_argument("name_prefix")
    ap.add_argument("--min-mapq", type=int, default=20, help="minimum PAF MAPQ (default 20)")
    a = ap.parse_args()

    stats = parse_paf(a.paf, a.min_mapq)
    seqs = {h.split()[0]: s for h, s in read_fasta(a.query_fasta)}
    unknown = set(stats) - set(seqs)
    if unknown:
        die(f"{len(unknown)} PAF query names not in {a.query_fasta}, e.g. {sorted(unknown)[:3]}")

    assign = {}  # contig -> (chr, match, frac, strand, midpoint)
    for q, per_chr in stats.items():
        total = sum(v[0] for v in per_chr.values())
        chrom, v = max(per_chr.items(), key=lambda kv: (kv[1][0], kv[0]))
        if v[0] == 0:
            continue
        strand = "-" if v[2] > v[1] else "+"
        assign[q] = (chrom, v[0], v[0] / total, strand, v[3] / v[0])

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    by_chr = defaultdict(list)
    for q, (chrom, match, frac, strand, mid) in assign.items():
        by_chr[chrom].append((mid, q))

    with open(a.out_prefix + ".assign.tsv", "w") as fa, open(a.out_prefix + ".boundaries.tsv", "w") as fb:
        fa.write("contig\tlength\tchrom\tmatch_bp\tfrac_of_contig_matches\tstrand\tweighted_midpoint\n")
        fb.write("pseudo_chrom\tcontig\tstrand\tstart_0based\tend\n")
        for chrom in sorted(by_chr, key=lambda c: (len(c), c)):
            name = f"{a.name_prefix}_{chrom}"
            parts, pos = [], 0
            for mid, q in sorted(by_chr[chrom]):
                _, match, frac, strand, _ = assign[q]
                seq = seqs[q] if strand == "+" else revcomp(seqs[q])
                fa.write(f"{q}\t{len(seq)}\t{chrom}\t{match}\t{frac:.4f}\t{strand}\t{mid:.0f}\n")
                fb.write(f"{name}\t{q}\t{strand}\t{pos}\t{pos + len(seq)}\n")
                parts.append(seq)
                pos += len(seq)
            out = f"{a.out_prefix}.{chrom}.fasta"
            with open(out + ".tmp", "w") as fh:
                write_record(fh, name, "".join(parts))
            os.replace(out + ".tmp", out)

    unassigned = [q for q in seqs if q not in assign]
    with open(a.out_prefix + ".unassigned.tsv", "w") as fh:
        fh.write("contig\tlength\n")
        for q in sorted(unassigned, key=lambda q: -len(seqs[q])):
            fh.write(f"{q}\t{len(seqs[q])}\n")

    total_bp = sum(len(s) for s in seqs.values())
    asg_bp = sum(len(seqs[q]) for q in assign)
    un_bp = sum(len(seqs[q]) for q in unassigned)
    if asg_bp + un_bp != total_bp:
        die("assigned + unassigned bp != total bp")
    with open(a.out_prefix + ".summary.tsv", "w") as fh:
        fh.write("chrom\tn_contigs\tbp\n")
        for chrom in sorted(by_chr, key=lambda c: (len(c), c)):
            fh.write(f"{chrom}\t{len(by_chr[chrom])}\t{sum(len(seqs[q]) for _, q in by_chr[chrom])}\n")
        fh.write(f"ASSIGNED\t{len(assign)}\t{asg_bp}\n")
        fh.write(f"UNASSIGNED\t{len(unassigned)}\t{un_bp}\n")
        fh.write(f"TOTAL\t{len(seqs)}\t{total_bp}\n")
    print(f"assigned {len(assign)} contigs ({asg_bp} bp), unassigned {len(unassigned)} ({un_bp} bp), total {total_bp} bp",
          file=sys.stderr)


if __name__ == "__main__":
    main()
