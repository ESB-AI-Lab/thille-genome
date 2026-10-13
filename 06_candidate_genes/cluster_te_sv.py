#!/usr/bin/env python3
"""TE and SV density around the two clusters of flowering-function candidate genes on Thille v4, against the genome.

Edwin (2026-10-08): do the flowering-gene clusters sit in sequence with more TEs or more structural variants than the
genomic background? Clusters (exact locus spans, manuscript Results): chr3 73,340,678-91,615,352 (CIGR2 to LIP1) and
chr5 19,577,582-30,337,161 (MTERF2 to PDIL2-1). Regions: the cluster interval, and flanks of 100 kb, 1 Mb and 5 Mb on
both sides together (truncated at chromosome ends; the truncation is reported).
  TE density = TE bp / region bp, each nucleotide counted once (union of the v4 RepeatMasker GFF intervals, without
               simple-repeat targets "(N)n" and low-complexity "X-rich" targets).
  SV density = SyRI events >= 50 bp per Mb with Thille coordinates in the region (compare/v4/sv_bed/thille.sv.all.bed;
               comparisons with Hass (UoQ), Unnamed Lowland and Puebla, each and pooled; Gwen-only categories excluded).
Background: every non-overlapping window of the region's length tiled along chr1-12 (partial windows at chromosome ends
dropped); the percentile is the share of background windows with a density at or below the region's.
Usage (from the repo root): python3 scripts/cluster_te_sv.py   # writes summaries/v4/cluster_te_sv.tsv
"""
import bisect
import csv
import re
from collections import defaultdict
from pathlib import Path

FAI = "genome/thille.hifiv4.asm.hic.scaffold.masked.chr.only.fasta.fai"
TE_GFF = "annotation/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff"
SV_BED = "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/compare/v4/sv_bed/thille.sv.all.bed"
CLUSTERS = {"chr3": ("chr3", 73_340_678, 91_615_352), "chr5": ("chr5", 19_577_582, 30_337_161)}
FLANKS = [100_000, 1_000_000, 5_000_000]
PARTNERS = {"hass": "Hass (UoQ)", "westindian": "Unnamed Lowland", "puebla": "Puebla"}
OUT = Path("summaries/v4/cluster_te_sv.tsv")
SIMPLE = re.compile(r"^\(.*\)n$|-rich$")  # RepeatMasker simple-repeat "(AAACCCT)n" and low-complexity "GA-rich" targets


def chrom_lengths():
    return {f[0]: int(f[1]) for f in (l.split("\t") for l in open(FAI)) if f[0] in {f"chr{i}" for i in range(1, 13)}}


def merged_te(lengths):
    iv = defaultdict(list)
    for l in open(TE_GFF):
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] not in lengths:
            continue
        name = f[8].split('Target "Motif:', 1)[-1].split('"', 1)[0]
        if SIMPLE.search(name):                                 # simple repeats and low complexity, not TEs
            continue
        iv[f[0]].append((int(f[3]) - 1, int(f[4])))          # GFF 1-based closed -> half-open
    out = {}
    for c, v in iv.items():
        v.sort(); m = []
        for s, e in v:
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        starts = [s for s, _ in m]
        cum = [0]
        for s, e in m:
            cum.append(cum[-1] + e - s)
        out[c] = (starts, [e for _, e in m], cum)
    return out


def te_bp(te, c, a, b):
    """Union TE bp inside [a, b)."""
    starts, ends, cum = te[c]
    i = bisect.bisect_right(ends, a)          # first interval ending after a
    j = bisect.bisect_left(starts, b)         # intervals starting before b
    if i >= j:
        return 0
    total = cum[j] - cum[i]
    total -= max(0, a - starts[i])            # clip the first interval
    total -= max(0, ends[j - 1] - b)          # clip the last interval
    return total


def sv_positions(lengths):
    pos = {p: defaultdict(list) for p in list(PARTNERS) + ["pooled"]}
    with open(SV_BED) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p = r["partner_genome"]
            if p not in PARTNERS or r["#chrom"] not in lengths or int(r["size_bp"]) < 50:
                continue
            mid = (int(r["start"]) + int(r["end"])) // 2
            pos[p][r["#chrom"]].append(mid); pos["pooled"][r["#chrom"]].append(mid)
    for p in pos:
        for c in pos[p]:
            pos[p][c].sort()
    return pos


def sv_count(pos, c, a, b):
    v = pos.get(c, [])
    return bisect.bisect_left(v, b) - bisect.bisect_left(v, a)


def background(lengths, size, te, sv):
    te_d, sv_d = [], defaultdict(list)
    for c, L in lengths.items():
        for a in range(0, L - size + 1, size):
            te_d.append(te_bp(te, c, a, a + size) / size)
            for p in sv:
                sv_d[p].append(sv_count(sv[p], c, a, a + size) / size * 1e6)
    return te_d, sv_d


def pct(values, x):
    return 100 * sum(v <= x for v in values) / len(values) if values else float("nan")


def main():
    lengths = chrom_lengths()
    te = merged_te(lengths)
    sv = sv_positions(lengths)
    genome_te = sum(te_bp(te, c, 0, L) for c, L in lengths.items()) / sum(lengths.values())
    rows = []
    for name, (c, s, e) in CLUSTERS.items():
        regions = [("cluster interval", [(s, e)], e - s)]
        for f in FLANKS:
            parts = [(max(0, s - f), s), (e, min(lengths[c], e + f))]
            got = sum(b - a for a, b in parts)
            regions.append((f"flanks {f // 1000:,} kb", parts, got))
        for label, parts, L in regions:
            tb = sum(te_bp(te, c, a, b) for a, b in parts)
            bg_te, bg_sv = background(lengths, L, te, sv)
            row = {"cluster": f"{name} {s:,}-{e:,}", "region": label, "bp": L,
                   "flank_truncated": "yes" if label.startswith("flanks") and L < 2 * int(label.split()[1].replace(",", "")) * 1000 else "no",
                   "te_fraction": f"{tb / L:.4f}", "te_percentile": f"{pct(bg_te, tb / L):.1f}",
                   "genome_te_fraction": f"{genome_te:.4f}", "n_background_windows": len(bg_te)}
            for p in list(PARTNERS) + ["pooled"]:
                n = sum(sv_count(sv[p], c, a, b) for a, b in parts)
                d = n / L * 1e6
                row[f"sv_{p}_n"] = n; row[f"sv_{p}_per_mb"] = f"{d:.1f}"; row[f"sv_{p}_percentile"] = f"{pct(bg_sv[p], d):.1f}"
            rows.append(row)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("# TE and SV density at the flowering-gene clusters vs same-size windows genome-wide (scripts/cluster_te_sv.py)\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r["cluster"], r["region"], r["bp"], "TE", r["te_fraction"], f"(pct {r['te_percentile']})",
              "SV/Mb", r["sv_pooled_per_mb"], f"(pct {r['sv_pooled_percentile']})", "trunc", r["flank_truncated"])
    print(f"genome TE fraction {genome_te:.4f}")


if __name__ == "__main__":
    main()
