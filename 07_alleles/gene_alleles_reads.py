#!/usr/bin/env python3
"""Phased read-based variants and CDS for one Fst flowering gene and one read set (Thille HiFi, WI HiFi, Gwen CLR).

Same method as SDMYB (scripts/sdmyb_haplotypes.py, scripts/indel_clusters.py, scripts/sdmyb_allele_tables.py), with
the gene's own CDS exons and strand from summaries/fst_alleles/gene_models.tsv:
  - phasing: heterozygous sites over the gene +- --phase-flank bp, two-means on read allele vectors (Gwen CLR:
    --clr noise masking, substitutions only); groups 0/1 -> haplotypes h1/h2
  - CDS per haplotype: majority-base consensus of that haplotype's reads (sdmyb_haplotypes.py)
  - variants, carried allele per haplotype (Edwin, 2026-10-01): among one haplotype's reads covering a position
    (>= 3, gaps excluded), the allele in >= --carry (50%) of them, or, when no allele reaches 50%, the single most
    frequent allele. A read set that cannot be phased is one diploid pool: the most frequent allele is h1, the second
    most frequent is h2 when it is in >= --second-af (25%) of the reads, otherwise h2 = h1. Indel clusters (>= 2 bp,
    kept when >= 25% of all spanning reads or >= 50% of one haplotype's reads carry them) follow the same rules with
    two states (indel / no indel). hap1_carried / hap2_carried record the call.
    Reads carrying / reads covering are given per haplotype. (Heterozygous-site detection alone misses variants
    carried by both haplotypes and variants of a haplotype with < 30% of the reads.)
When phasing gives fewer than two groups (no heterozygous site, or no two reads sharing enough of them), phasing is
rerun with --min-sites 0 so every read is in h1, h2 repeats h1, and the read set is reported "not phased".

Usage (from the repo root):
    python3 scripts/gene_alleles_reads.py <label> <Thille|WI|Gwen>
Writes summaries/fst_alleles/<label>/<set>.{tsv,fa,groups.tsv,indels.tsv,variants.tsv}
"""
import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sdmyb_haplotypes as H  # noqa: E402
from fst_common import OUT, load_models  # noqa: E402

M = "/epool/edwin/pamericana/avoflow/mappingback"
# THILLE_ASSEMBLY=v4 reads the alignments to the v4 assembly (complete chr8, 2026-10-04); default v3 as published so far
SETS = {"v3": {"Thille": (M + "/minimap2/thille.hifi.thille.{chrom}.bam", "HiFi", False),
               "WI": (M + "/minimap2/westindian.hifi.thille.{chrom}.bam", "HiFi", False),
               "Gwen": (M + "/gwen.clr.thille.{chrom}.bam", "CLR", True)},
        "v4": {"Thille": (M + "/v4/thille.hifi.thille_v4.{chrom}.bam", "HiFi", False),
               "WI": (M + "/v4/westindian.hifi.thille_v4.{chrom}.bam", "HiFi", False),
               "Gwen": (M + "/v4/gwen.clr.thille_v4.{chrom}.bam", "CLR", True)}}[os.environ.get("THILLE_ASSEMBLY", "v3")]
PHASE_REGION = {"SDMYB": "chr10:10309000-10322000"}  # the region used for the published SDMYB analysis


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("label"); p.add_argument("read_set", choices=SETS)
    p.add_argument("--phase-flank", type=int, default=2000)
    p.add_argument("--carry", type=float, default=0.5)
    p.add_argument("--second-af", type=float, default=0.25,
                   help="unphased read sets: the second most frequent allele goes on h2 when in >= this fraction of reads")
    p.add_argument("--reuse-phasing", action="store_true", help="keep existing <set>.tsv/.fa/.groups.tsv")
    a = p.parse_args()
    g = load_models()[a.label]
    bam_t, tech, clr = SETS[a.read_set]
    bam = bam_t.format(chrom=g.chrom)
    d = OUT / a.label; d.mkdir(parents=True, exist_ok=True)
    region = PHASE_REGION.get(a.label, f"{g.chrom}:{max(1, g.gene_start - a.phase_flank)}-{g.gene_end + a.phase_flank}")
    window = f"{g.chrom}:{g.win_start}-{g.win_end}"
    cds_arg = ",".join(f"{s}-{e}" for s, e in g.exons) if g.exons else None
    cmd = [sys.executable, "scripts/sdmyb_haplotypes.py", bam, a.read_set, "--region", region, "--outdir", str(d),
           "--strand", g.strand]
    if cds_arg:
        cmd += ["--cds", cds_arg]
    else:  # no CDS model: phase only; the CDS output is meaningless and ignored
        cmd += ["--cds", f"{g.gene_start}-{g.gene_start + 2}"]
    if clr:  # CLR: seed reads must share >= 20 het sites; on SDMYB this reproduces the CDS-SNP seeded grouping (88/88 reads)
        cmd += ["--clr", "--seed-min-shared", "20"]
    gfile = d / f"{a.read_set}.groups.tsv"
    if not (a.reuse_phasing and gfile.exists()):
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    if len({l.split()[1] for l in open(gfile)}) < 2:
        subprocess.run(cmd + ["--min-sites", "0"], check=True, capture_output=True, text=True)
    subprocess.run([sys.executable, "scripts/indel_clusters.py", bam, window, a.read_set, "--groups",
                    str(d / f"{a.read_set}.groups.tsv"), "--outdir", str(d), "--group-min-support", str(a.carry)],
                   check=True, capture_output=True, text=True)

    ref = "".join(H.sh("samtools", "faidx", H.REF, window).splitlines()[1:]).upper()
    refb = lambda pos: ref[pos - g.win_start]
    reads = H.read_events(bam, window)
    grp = {l.split()[0]: int(l.split()[1]) + 1 for l in open(d / f"{a.read_set}.groups.tsv")}
    unphased = len(set(grp.values())) < 2
    sites = [(int(l.split("\t")[2]), l.split("\t")[3], l.split("\t")[4].strip())
             for l in open(d / f"{a.read_set}.tsv") if l.startswith("site\tB")]
    rows = []
    hap_reads = {h: [r for n, r in reads.items() if grp.get(n) == h] for h in (1, 2)}
    # positions each haplotype covers with >= 3 reads, so the figure can show "no data" instead of "reference"
    with open(d / f"{a.read_set}.covered.tsv", "w") as fh:
        fh.write("hap\tstart\tend\n")
        for h in (1, 2):
            src = hap_reads[1] if unphased else hap_reads[h]
            diff = [0] * (g.win_end - g.win_start + 2)
            for r in src:
                s, e = max(r[0], g.win_start), min(r[1], g.win_end)
                if s <= e:
                    diff[s - g.win_start] += 1; diff[e - g.win_start + 1] -= 1
            depth, start = 0, None
            for i, x in enumerate(diff):
                depth += x
                if depth >= 3 and start is None:
                    start = i
                elif depth < 3 and start is not None:
                    fh.write(f"{h}\t{g.win_start + start}\t{g.win_start + i - 1}\n"); start = None
    for pos in range(g.win_start, g.win_end + 1):
        rb = refb(pos)
        cnt = {}
        for h in (1, 2):
            c = {}
            for r in hap_reads[h]:
                if r[0] <= pos <= r[1]:
                    b = r[2].get(pos)
                    if b and b != "-":
                        c[b] = c.get(b, 0) + 1
            cnt[h] = c
        allele = {}
        if unphased:  # one diploid read pool: most frequent allele on h1, second on h2 when >= --second-af
            c = cnt[1]; cov = sum(c.values())
            if cov >= 3:
                ranked = sorted(c.items(), key=lambda kv: -kv[1])
                allele[1] = ranked[0][0]
                allele[2] = ranked[1][0] if len(ranked) > 1 and ranked[1][1] / cov >= a.second_af else ranked[0][0]
        else:  # each phased haplotype: allele in >= --carry of its reads, else its single most frequent allele
            for h in (1, 2):
                c = cnt[h]; cov = sum(c.values())
                if cov < 3:
                    continue
                ranked = sorted(c.items(), key=lambda kv: -kv[1])
                if ranked[0][1] / cov >= a.carry or len(ranked) == 1 or ranked[0][1] > ranked[1][1]:
                    allele[h] = ranked[0][0]
        for alt in sorted({b for b in allele.values() if b != rb}):
            rec = {"read_set": a.read_set, "technology": tech, "pos": pos, "region": g.region_of(pos), "type": "SNP",
                   "ref": rb, "alt": alt, "length": 1}
            for h in (1, 2):
                c = cnt[1] if unphased else cnt[h]
                rec[f"hap{h}_with"] = c.get(alt, 0)
                rec[f"hap{h}_covering"] = sum(c.values())
                rec[f"hap{h}_carried"] = int(allele.get(h) == alt)
            rows.append(rec)
    for r in csv.DictReader((l for l in open(d / f"{a.read_set}.indels.tsv") if not l.startswith("#")), delimiter="\t"):
        w1, c1 = map(int, r["group0_with/spanning"].split("/"))
        w2, c2 = map(int, r["group1_with/spanning"].split("/"))
        rec = {"read_set": a.read_set, "technology": tech, "pos": int(r["pos"]), "region": g.region_of(int(r["pos"])),
               "type": r["type"], "ref": "", "alt": "", "length": int(r["len"]),
               "hap1_with": w1, "hap1_covering": c1, "hap2_with": w2, "hap2_covering": c2}
        if unphased:  # one diploid read pool: majority state on h1, minority state on h2 when >= --second-af
            f = w1 / c1 if c1 else 0.0
            rec["hap1_carried"] = int(c1 > 0 and f >= 0.5)
            rec["hap2_carried"] = int(c1 > 0 and ((f >= 0.5 and 1 - f < a.second_af) or (f < 0.5 and f >= a.second_af)))
            rec["hap2_with"], rec["hap2_covering"] = w1, c1
        else:  # two states (indel / no indel): the most frequent is the one in >= 50%
            for h, (w, c) in ((1, (w1, c1)), (2, (w2, c2))):
                rec[f"hap{h}_carried"] = int(c > 0 and w / c >= a.carry)
        if rec["hap1_carried"] or rec["hap2_carried"]:
            rows.append(rec)
    cols = ["read_set", "technology", "pos", "region", "type", "ref", "alt", "length", "hap1_with", "hap1_covering",
            "hap2_with", "hap2_covering", "hap1_carried", "hap2_carried", "unphased"]
    with open(d / f"{a.read_set}.variants.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["pos"]):
            r["unphased"] = int(unphased)
            w.writerow(r)
    n_snp = sum(r["type"] == "SNP" for r in rows)
    print(f"{a.label} {a.read_set}: {len(sites)} het SNP sites; carried SNPs {n_snp}, indel clusters {len(rows) - n_snp}, "
          f"groups {sorted(set(grp.values()))}{' (not phased)' if unphased else ''}")


if __name__ == "__main__":
    main()
