#!/usr/bin/env python3
"""Check loci with skewed exon allele frequencies (one phased haplotype at 60 to 80% of reads) for read-level signs of a
duplication on the majority haplotype or a deletion on the minority haplotype.

For each phased read set (Thille HiFi, Lowland/WI HiFi, Gwen CLR) whose exon AF (summaries/fst_alleles/allele_summary.tsv)
is 60 to 80% for one haplotype:
  - per-haplotype depth: phased reads of each haplotype spanning points every 250 bp over the gene and over flanks
    2 to 6 kb outside it; ratio gene / flank per haplotype and for all reads
  - clipped reads: alignments (primary and supplementary) over gene +- 10 kb with a soft or hard clip >= --min-clip bp;
    clip positions are clustered within 100 bp and a breakpoint is reported when >= 3 reads share it
Expected signatures: an extra copy on one of the two haplotypes raises the all-read depth over the gene to about 1.5x
the flanks, and loss of one haplotype lowers it to about 0.5x, with clip clusters at the junctions. Per-haplotype
ratios are reported but are biased upward (reads are phased only near the gene) and are not used for the call. Without
a depth signal, a 60 to 80% AF fits unequal recovery of the two haplotypes (for example fewer reads of a divergent
haplotype aligning to Thille), a somatic or tissue-mosaic allele, or phasing imbalance. Observation per locus, not a test.

Usage (from the repo root): python3 scripts/af_skew_check.py [--min-clip 200]
Writes summaries/fst_alleles/af_skew_check.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sdmyb_haplotypes as H  # noqa: E402
from fst_common import OUT, load_models  # noqa: E402
from gene_alleles_reads import SETS  # noqa: E402


def mean(x):
    return sum(x) / len(x) if x else float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--min-clip", type=int, default=200)
    p.add_argument("--lo", type=float, default=60); p.add_argument("--hi", type=float, default=80)
    a = p.parse_args()
    models = load_models()
    targets = defaultdict(dict)
    for r in csv.DictReader(open(OUT / "allele_summary.tsv"), delimiter="\t"):
        cv, h = r["haplotype"].split(" h")
        if cv in SETS and r["allele_freq_pct"] not in ("", "n.d."):
            targets[(r["label"], cv)][int(h)] = float(r["allele_freq_pct"])
    out = []
    for (lab, rs), af in sorted(targets.items()):
        if len(af) < 2 or not any(a.lo <= v <= a.hi for v in af.values()):
            continue
        d = OUT / lab
        g = models[lab]
        vt = list(csv.DictReader(open(d / f"{rs}.variants.tsv"), delimiter="\t"))
        if vt and str(vt[0].get("unphased", "")) in ("1", "1.0"):
            continue  # unphased: no haplotype groups to compare
        grp = {l.split()[0]: int(l.split()[1]) + 1 for l in open(d / f"{rs}.groups.tsv")}
        bam = SETS[rs][0].format(chrom=g.chrom)
        lo, hi = max(1, g.gene_start - 10000), g.gene_end + 10000
        reads = H.read_events(bam, f"{g.chrom}:{lo}-{hi}")
        gene_pts = list(range(g.gene_start, g.gene_end + 1, 250)) or [g.gene_start]
        flank_pts = [x for x in range(g.gene_start - 6000, g.gene_start - 1999, 250)] + \
                    [x for x in range(g.gene_end + 2000, g.gene_end + 6001, 250)]
        depth = lambda names, pts: mean([sum(1 for n in names if reads[n][0] <= x <= reads[n][1]) for x in pts])
        res = {"label": lab, "gene": g.gene, "read_set": rs, "af_h1": af[1], "af_h2": af[2]}
        for h in (1, 2):
            names = [n for n in reads if grp.get(n) == h]
            dg, df = depth(names, gene_pts), depth(names, flank_pts)
            res[f"h{h}_gene_depth"] = round(dg, 1); res[f"h{h}_flank_depth"] = round(df, 1)
            res[f"h{h}_gene_over_flank"] = round(dg / df, 2) if df else ""
        allr = list(reads)
        dg, df = depth(allr, gene_pts), depth(allr, flank_pts)
        res["all_gene_depth"] = round(dg, 1); res["all_flank_depth"] = round(df, 1)
        res["all_gene_over_flank"] = round(dg / df, 2) if df else ""
        # clipped alignments, primary and supplementary
        clips = []
        for line in H.sh("samtools", "view", "-F", "0x104", bam, f"{g.chrom}:{lo}-{hi}").splitlines():
            x = line.split("\t")
            pos, cig = int(x[3]), x[5]
            ops = re.findall(r"(\d+)([MIDNSHX=])", cig)
            if not ops:
                continue
            ref_len = sum(int(n) for n, o in ops if o in "MDN=X")
            if ops[0][1] in "SH" and int(ops[0][0]) >= a.min_clip:
                clips.append((pos, x[0]))
            if ops[-1][1] in "SH" and int(ops[-1][0]) >= a.min_clip:
                clips.append((pos + ref_len - 1, x[0]))
        clips.sort()
        clusters, cur = [], []
        for c in clips:
            if cur and c[0] - cur[-1][0] > 100:
                clusters.append(cur); cur = []
            cur.append(c)
        if cur:
            clusters.append(cur)
        bps = []
        for c in clusters:
            names = {n for _, n in c}
            if len(names) >= 3:
                hp = defaultdict(int)
                for n in names:
                    hp[grp.get(n, 0)] += 1
                bps.append(f"{c[len(c) // 2][0]}:{len(names)}reads(h1={hp[1]},h2={hp[2]},unphased={hp[0]})")
        res["clip_breakpoints"] = ";".join(bps) or "none"
        # per-haplotype gene/flank ratios are biased upward: reads are assigned a haplotype only when they cover
        # heterozygous sites near the gene, so phased reads are depleted in the flanks. The call uses all reads:
        # an extra copy on one of two haplotypes gives about 1.5, loss of one haplotype about 0.5.
        ar = res["all_gene_over_flank"]
        if ar != "" and ar >= 1.3:
            call = "extra copy (all-read depth >= 1.3x flanks)"
        elif ar != "" and ar <= 0.75:
            call = "one haplotype missing over part of the gene (all-read depth <= 0.75x flanks)"
        else:
            call = "no depth signal"
        res["call"] = call + ("; clip clusters present" if bps else "")
        out.append(res)
        print({k: res[k] for k in ("label", "read_set", "af_h1", "af_h2", "h1_gene_over_flank", "h2_gene_over_flank",
                                   "all_gene_over_flank", "call")}, "|", res["clip_breakpoints"][:160], flush=True)
    cols = list(out[0])
    with open(OUT / "af_skew_check.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(out)


if __name__ == "__main__":
    main()
