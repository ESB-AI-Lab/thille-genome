#!/usr/bin/env python3
"""Alignment-free test of the SDMYB alleles in raw reads (21-mers), reproducing summaries/sdmyb_kmers/ (commit
2606426, 2026-09-30) and adding Puebla HiFi reads (Edwin, 2026-10-11).

Design (subcommand `design`): for each site, every 21-mer that overlaps the site, on both strands, from
  allele 1 = the Thille sequence (chr10; v3 and v4 are identical here), and
  allele 2 = the Thille sequence with the Gwen haplotype-2 substitutions and deletions applied (calls with
             hap2_with / hap2_covering >= 0.5 in summaries/sdmyb_alleles/_pre_gwen_correction/variants.tsv, the
             allele-2 calls of 2026-09-30).
The ten coding substitutions rebuild the saved design (summaries/sdmyb_kmers/kmers.tsv) exactly; the eleventh site,
an exon-3 2 bp deletion of the earlier call that the corrected allele places elsewhere, does not, and the manuscript
uses the ten substitutions only. Counting uses the saved design (kmers.txt), so new counts are directly comparable.

Counting (sbatch_sdmyb_kmer_test.sh): zcat <reads> | LC_ALL=C grep -o -F -f kmers.txt | sort | uniq -c, i.e. the
non-overlapping matches grep -o reports per read line, as in the original run.

Summary (subcommand `summarize`): per read set and site, the summed counts of the allele-1 and allele-2 21-mers, and
totals over the ten substitutions (manuscript: Unnamed Lowland 701/0, Thille 974/1, Gwen 149/78).

Usage (from the repo root):
  python3 scripts/sdmyb_kmer_test.py design            # rebuilds the design and checks it against kmers.tsv
  python3 scripts/sdmyb_kmer_test.py summarize --counts-dir summaries/sdmyb_kmers_rerun --out summaries/sdmyb_kmers_rerun/allele_kmer_counts.tsv
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

SAVED = Path("summaries/sdmyb_kmers")
FA = "genome/thille.hifiv4.asm.hic.scaffold.masked.chr.only.fasta"
VARIANTS = "summaries/sdmyb_alleles/_pre_gwen_correction/variants.tsv"
ST = "/gpool/bin/samtools-1.24/samtools"
K = 21
DELETION_SITE = 10311778   # exon-3 2 bp deletion of the earlier allele-2 call; not used in the manuscript
NAMES = {10311778: "2bp del (exon 3)"}
READ_SETS = ["Thille_hifi", "WI_hifi", "Gwen_clr", "Puebla_hifi"]


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def window(lo, hi):
    out = subprocess.run([ST, "faidx", FA, f"chr10:{lo}-{hi}"], capture_output=True, text=True, check=True).stdout
    return "".join(out.split("\n")[1:]).upper()


def kmers_over(seq, pos, site):
    ks = {seq[i:i + K] for i in range(len(seq) - K + 1) if pos[i] <= site <= pos[i + K - 1]}
    return ks | {rc(k) for k in ks}


def design(_):
    saved = list(csv.DictReader(open(SAVED / "kmers.tsv"), delimiter="\t"))
    calls = [r for r in csv.DictReader(open(VARIANTS), delimiter="\t") if r["read_set"] == "Gwen"]
    hap2 = lambda r: int(r["hap2_covering"] or 0) and int(r["hap2_with"] or 0) / int(r["hap2_covering"]) >= 0.5
    snp = {int(r["pos"]): r["alt"] for r in calls if r["type"] == "SNP" and hap2(r)}
    dels = [(int(r["pos"]), int(r["length"])) for r in calls if r["type"] == "DEL" and hap2(r)]
    ok = True
    for site in sorted({int(r["site"]) for r in saved}):
        lo, hi = site - 2 * K, site + 2 * K
        ref = window(lo, hi)
        pos = list(range(lo, hi + 1))
        a1 = kmers_over(ref, pos, site)
        alt = [snp.get(p, b) for p, b in zip(pos, ref)]
        gone = {q for p, l in dels for q in range(p + 1, p + 1 + l)}
        s2 = "".join(b for p, b in zip(pos, alt) if p not in gone)
        p2 = [p for p in pos if p not in gone]
        a2 = kmers_over(s2, p2, site)
        w1 = {r["kmer"] for r in saved if int(r["site"]) == site and r["allele"].startswith("allele1")}
        w2 = {r["kmer"] for r in saved if int(r["site"]) == site and r["allele"].startswith("allele2")}
        same = (a1 == w1, a2 == w2)
        print(f"{site}\tallele1 {len(a1 & w1)}/{len(w1)}\tallele2 {len(a2 & w2)}/{len(w2)}"
              + ("" if all(same) else "\tdiffers" + (" (expected: earlier deletion call)" if site == DELETION_SITE else "")))
        if site != DELETION_SITE and not all(same):
            ok = False
    if not ok:
        sys.exit("design does not reproduce the saved k-mers at a substitution site")
    print("design reproduces the saved k-mers at all ten substitution sites")


def summarize(a):
    saved = list(csv.DictReader(open(SAVED / "kmers.tsv"), delimiter="\t"))
    labels = {int(r["pos"]): r["variant"] for r in csv.DictReader(open(SAVED / "allele_kmer_counts.tsv"), delimiter="\t")}
    rows, totals = [], []
    for rs in READ_SETS:
        f = Path(a.counts_dir) / f"counts_{rs}.txt"
        if not f.is_file():
            print(f"missing {f}", file=sys.stderr)
            continue
        c = {}
        for line in open(f):
            n, k = line.split()
            c[k] = int(n)
        t1 = t2 = 0
        for site in sorted({int(r["site"]) for r in saved}, key=lambda s: (s != DELETION_SITE, s)):
            h1 = sum(c.get(r["kmer"], 0) for r in saved if int(r["site"]) == site and r["allele"].startswith("allele1"))
            h2 = sum(c.get(r["kmer"], 0) for r in saved if int(r["site"]) == site and r["allele"].startswith("allele2"))
            rows.append([rs, labels.get(site, NAMES.get(site, "")), site, h1, h2, f"{h2 / (h1 + h2):.3f}" if h1 + h2 else "n/a"])
            if site != DELETION_SITE:
                t1 += h1; t2 += h2
        totals.append((rs, t1, t2))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("read_set\tvariant\tpos\tallele1_Thille_kmer_hits\tallele2_Gwen_kmer_hits\tallele2_fraction\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    for rs, t1, t2 in totals:
        print(f"{rs}: ten substitutions, SDMYB^a {t1}, SDMYB^b {t2}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("design").set_defaults(fn=design)
    s = sub.add_parser("summarize")
    s.add_argument("--counts-dir", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=summarize)
    a = p.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
