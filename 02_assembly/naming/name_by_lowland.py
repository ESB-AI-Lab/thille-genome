#!/usr/bin/env python3
"""Name the curated Thille scaffolds by the Lowland accession T2T chromosomes (GCA_051132755.1, Pa01-Pa12).

A scaffold of at least --min-chrom-len bp is named chrN when the Pa chromosome it aligns to most (union of aligned
scaffold bases, MAPQ >= --min-mapq) covers at least --min-frac of its length, and no two scaffolds claim the same Pa
chromosome; PaNN -> chrN. Every other scaffold becomes unplaced_scaffold_N, keeping N from the juicer post AGP so it
traces back. Orientation is not changed; the fraction of aligned bases on the plus strand is reported.

Usage (from curation/naming): name_by_lowland.py [--min-chrom-len 10000000] [--min-frac 0.5] [--min-mapq 1]
Writes, next to the curated FASTA: thille.hifiv4.asm.hic.p_ctg.fasta and .agp; here: curated_v2_names.tsv.
"""
import argparse
import collections
import os
import sys

CUR_DIR = "/epool/edwin/pamericana/avoflow/hifiasm/thille/hifi.hic.gcluster/yahs_thille.hifiv3.asm.hic.p_ctg.curation"
FA = f"{CUR_DIR}/out_JBAT.v2.FINAL.fa"
AGP = f"{CUR_DIR}/out_JBAT.v2.FINAL.agp"
OUT_FA = f"{CUR_DIR}/thille.hifiv4.asm.hic.p_ctg.fasta"
OUT_AGP = f"{CUR_DIR}/thille.hifiv4.asm.hic.p_ctg.agp"
ORDER = "../../compare/westindian/westindian_chroms.order.tsv"
PAF = "curated_v2_vs_lowland.asm5.paf"


def union_len(iv):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if ce is not None else 0)


def read_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-chrom-len", type=int, default=10_000_000)
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--min-mapq", type=int, default=1)
    a = ap.parse_args()

    for out in (OUT_FA, OUT_AGP):
        if os.path.exists(out):
            sys.exit(f"{out} exists; not overwriting")
    pa = {}
    for line in list(open(ORDER))[1:]:
        acc, _, grp = line.split()
        pa[acc] = grp
    seqs = dict(read_fasta(FA))
    order = list(seqs)
    print(f"{len(seqs)} curated scaffolds, {sum(map(len, seqs.values())):,} bp", file=sys.stderr)

    iv = collections.defaultdict(lambda: collections.defaultdict(list))
    plus = collections.Counter()
    alnb = collections.Counter()
    for line in open(PAF):
        f = line.split("\t")
        if int(f[11]) < a.min_mapq:
            continue
        q, qs, qe, strand, t = f[0], int(f[2]), int(f[3]), f[4], f[5]
        iv[q][t].append((qs, qe))
        alnb[q] += qe - qs
        if strand == "+":
            plus[q] += qe - qs

    rows, claimed = [], collections.defaultdict(list)
    for q in order:
        L = len(seqs[q])
        cov = sorted(((union_len(v), t) for t, v in iv[q].items()), reverse=True)
        top = cov[0] if cov else (0, "")
        sec = cov[1] if len(cov) > 1 else (0, "")
        pfrac = plus[q] / alnb[q] if alnb[q] else float("nan")
        is_chr = L >= a.min_chrom_len and top[0] / L >= a.min_frac
        if is_chr:
            claimed[top[1]].append(q)
        rows.append([q, L, pa.get(top[1], ""), top[0] / L, pa.get(sec[1], ""), sec[0] / L, pfrac, is_chr])

    dup = {t: qs for t, qs in claimed.items() if len(qs) > 1}
    if dup or len(claimed) != 12:
        sys.exit(f"naming is not one to one: {len(claimed)} Pa chromosomes claimed; shared: {dup}")

    new = {}
    for r in rows:
        num = r[0].split("_")[-1]
        new[r[0]] = f"chr{int(r[2][2:])}" if r[7] else f"unplaced_scaffold_{num}"

    with open("curated_v2_names.tsv", "w") as fh:
        fh.write("curated_name\tnew_name\tlength\ttop_lowland\ttop_frac\tsecond_lowland\tsecond_frac\tplus_strand_frac\n")
        for r in rows:
            fh.write(f"{r[0]}\t{new[r[0]]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]}\t{r[5]:.4f}\t{r[6]:.3f}\n")
    for r in rows:
        if r[7]:
            print(f"{r[0]:>14} -> {new[r[0]]:<5} {r[1]:>12,} bp  {r[2]} {r[3]:.1%}  second {r[4] or '-'} {r[5]:.1%}  "
                  f"plus strand {r[6]:.1%}", file=sys.stderr)

    chroms = sorted((q for q in order if new[q].startswith("chr")), key=lambda q: int(new[q][3:]))
    with open(OUT_FA + ".tmp", "w") as fh:
        for q in chroms + [q for q in order if not new[q].startswith("chr")]:
            fh.write(f">{new[q]}\n")
            s = seqs[q]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    with open(OUT_AGP + ".tmp", "w") as fh:
        for line in open(AGP):
            f = line.split("\t")
            f[0] = new[f[0]]
            fh.write("\t".join(f))
    os.rename(OUT_FA + ".tmp", OUT_FA)
    os.rename(OUT_AGP + ".tmp", OUT_AGP)
    print(f"wrote {OUT_FA}, {OUT_AGP}, curated_v2_names.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
