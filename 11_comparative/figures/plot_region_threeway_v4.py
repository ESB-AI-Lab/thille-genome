#!/usr/bin/env python3
# v4 copy (2026-10-05): Thille v4 genome, annotation, coverage, SyRI and centromere inputs; display name "Unnamed Lowland" (Edwin, 2026-10-05).
"""Three-genome view of one SV locus: West Indian (top), Thille (middle), Hass UoQ (bottom).

Thille window = the Thille SV block +/- --flank. The Hass and West Indian windows are the span of their SyRI alignments
to that Thille window on the same chromosome (Thille -> Hass and Thille -> West Indian runs in compare/syri/). Each genome
is drawn on its own scale, left to right in the same-strand orientation SyRI used (chromosomes reverse-complemented for
SyRI are drawn reverse-complemented; tick labels give original GenBank coordinates). Ribbons are SyRI alignment rows:
grey syntenic, orange inverted, green translocated, purple duplicated. The SV blocks of compare/regions/region_blocks.tsv
are shaded red. Under each chromosome: genes (Thille Tiberius, Hass UoQ GenBank genes, West Indian GenBank GFF), tandem
repeat arrays with period >= 100 bp (TRF .dat files), and Gwen CLR read depth in 20 kb bins relative to the chromosome
median (compare/regions/region_gwen_depth.tsv), one panel per genome.

Usage (from assemblies/thille.hic):
  plot_region_threeway.py <region_id> <out_prefix> --trf-hass H.dat --trf-wi W.dat [--flank 1500000] [--bin 20000]
The Hass/WI TRF files cover the plotted windows; their names end in _<start>.dat (0-based window start, same-strand).
"""
import argparse
import csv
import re
import statistics
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon, Rectangle  # noqa: E402

SAM = "/gpool/bin/samtools-1.24/samtools"
MB = "/epool/edwin/pamericana/avoflow/mappingback"
THILLE_GTF = "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf"
HASS_GENES = "/epool/edwin/pamericana/avoflow/z_GCA_029852735.1/ncbi_dataset/data/GCA_029852735.1/GCA_029852735.1_ASM2985273v1_genes.bed"
WI_GFF = "westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/genomic.gff"
RC = {"chr1", "chr4", "chr5", "chr6", "chr10", "chr11"}
IDS = {"hass": {f"chr{i + 1}": x for i, x in enumerate(
           "CM056809.1 CM056810.1 CM056811.1 CM056812.1 CM056813.1 CM056814.1 CM056815.1 CM056816.1 CM056817.1 "
           "CM056818.1 CM056819.1 CM056820.1".split())},
       "westindian": {f"chr{i + 1}": x for i, x in enumerate(
           "CP152278.1 CP152279.1 CP152280.1 CP152281.1 CP152282.1 CP152283.1 CP152284.1 CP152285.1 CP152286.1 "
           "CP152287.1 CP152288.1 CP152289.1".split())}}
FAI = {"thille": "thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai",
       "hass": "compare/v4/syri/hass/hass.chr.samestrand.fasta.fai",
       "westindian": "compare/v4/syri/westindian/westindian.chr.samestrand.fasta.fai"}
BAM = {"thille": "v4/gwen.clr.thille_v4", "hass": "gwen.clr.hass_uoq", "westindian": "gwen.clr.wi"}
NAME = {"westindian": "Unnamed Lowland", "thille": "Thille", "hass": "Hass UoQ"}
KIND = {"SYNAL": ("#b9b9b9", "syntenic"), "INVAL": ("#e08a00", "inverted"), "TRANSAL": ("#2e9e44", "translocated"),
        "INVTRAL": ("#2e9e44", "translocated"), "DUPAL": ("#8e44ad", "duplicated"), "INVDPAL": ("#8e44ad", "duplicated")}
INK, BLOCK, GENE, SAT = "#1f1f1f", "#c0392b", "#1f5fa8", "#6d4c41"
plt.rcParams.update({"font.size": 8, "text.color": INK, "pdf.fonttype": 42})


def length(genome, chrom):
    return {l.split()[0]: int(l.split()[1]) for l in open(FAI[genome])}[chrom]


def to_same(genome, chrom, s, e):
    """Original 0-based [s,e) on the GenBank sequence -> same-strand coordinates used by SyRI."""
    if genome != "thille" and chrom in RC:
        L = length(genome, chrom)
        return L - e, L - s
    return s, e


def to_orig(genome, chrom, x):
    return length(genome, chrom) - x if genome != "thille" and chrom in RC else x


def syri_rows(genome, chrom, ws, we):
    for line in open(f"compare/v4/syri/{genome}/{genome}syri.out"):
        f = line.split("\t")
        if f[0] != chrom or f[10] not in KIND or f[5] != chrom:
            continue
        rs, re_ = int(f[1]) - 1, int(f[2])
        if re_ > ws and rs < we:
            yield rs, re_, min(int(f[6]), int(f[7])) - 1, max(int(f[6]), int(f[7])), f[10], int(f[6]) > int(f[7])


def genes(genome, chrom, s, e):
    out = []
    if genome == "thille":
        for l in open(THILLE_GTF):
            f = l.split("\t")
            if len(f) > 4 and f[2] == "gene" and f[0] == chrom and int(f[3]) - 1 < e and int(f[4]) > s:
                out.append((int(f[3]) - 1, int(f[4])))
        return out
    oid = IDS[genome][chrom]
    if genome == "hass":
        for l in open(HASS_GENES):
            f = l.split("\t")
            if f[0] == oid:
                out.append(to_same(genome, chrom, int(f[1]), int(f[2])))
    else:
        for l in open(WI_GFF):
            f = l.split("\t")
            if len(f) > 4 and f[2] == "gene" and f[0] == oid:
                out.append(to_same(genome, chrom, int(f[3]) - 1, int(f[4])))
    return [g for g in out if g[0] < e and g[1] > s]


def trf_arrays(path, offset):
    out = []
    for l in open(path):
        if l.startswith("@"):
            continue
        p = l.split()
        if int(p[2]) >= 100:
            out.append((int(p[0]) - 1 + offset, int(p[1]) + offset))
    return out


def depth(genome, chrom, s, e, binsize, median):
    """Gwen CLR depth / chromosome median in bins; BAMs use original ids and coordinates."""
    oid = chrom if genome == "thille" else IDS[genome][chrom]
    bins = [(b, min(b + binsize, e)) for b in range(s, e, binsize)]
    bed = "".join(f"{oid}\t{a}\t{b}\n" for a, b in (sorted(to_same(genome, chrom, x, y)) for x, y in bins))
    r = subprocess.run([SAM, "bedcov", "/dev/stdin", f"{MB}/{BAM[genome]}.bam"], input=bed, capture_output=True,
                       text=True, check=True)
    vals = [int(l.split("\t")[3]) / (int(l.split("\t")[2]) - int(l.split("\t")[1])) / median
            for l in r.stdout.splitlines()]
    return [(a + b) / 2 for a, b in bins], vals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("region"); ap.add_argument("out")
    ap.add_argument("--trf-hass", required=True); ap.add_argument("--trf-wi", required=True)
    ap.add_argument("--flank", type=int, default=1500000); ap.add_argument("--bin", type=int, default=20000)
    a = ap.parse_args()
    blocks = [r for r in csv.DictReader(open("compare/regions/region_blocks.tsv"), delimiter="\t")
              if r["region"] == a.region]
    med = {(r["region"], r["genome"]): float(r["chrom_median_depth"])
           for r in csv.DictReader(open("compare/regions/region_gwen_depth.tsv"), delimiter="\t")}
    tb = next(b for b in blocks if b["genome"] == "thille")
    chrom = tb["chrom"]
    win = {"thille": (max(0, int(tb["start"]) - a.flank), min(length("thille", chrom), int(tb["end"]) + a.flank))}
    ribbons = {}
    for g in ("westindian", "hass"):
        rows = list(syri_rows(g, chrom, *win["thille"]))
        ribbons[g] = rows
        qs = [r[2] for r in rows if r[4] == "SYNAL"] + [r[3] for r in rows if r[4] == "SYNAL"]
        win[g] = (min(qs), max(qs))
    block = {b["genome"]: to_same(b["genome"], chrom, int(b["start"]), int(b["end"])) for b in blocks}
    trf = {"thille": trf_arrays(f"centromere/v4/trf/{chrom}.dat", 0),
           "hass": trf_arrays(a.trf_hass, int(re.search(r"_(\d+)\.dat$", a.trf_hass).group(1))),
           "westindian": trf_arrays(a.trf_wi, int(re.search(r"_(\d+)\.dat$", a.trf_wi).group(1)))}

    order = ["westindian", "thille", "hass"]
    ylev = {"westindian": 2.0, "thille": 1.0, "hass": 0.0}
    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(4, 1, height_ratios=[3.2, 1, 1, 1], hspace=0.55)
    ax = fig.add_subplot(gs[0])

    def X(g, x):
        s, e = win[g]
        return (x - s) / (e - s)

    for g in ("westindian", "hass"):
        for rs, re_, qs, qe, kind, rev in ribbons[g]:
            t0, t1 = X("thille", max(rs, win["thille"][0])), X("thille", min(re_, win["thille"][1]))
            q0, q1 = X(g, qs), X(g, qe)
            if kind in ("INVAL", "INVTRAL", "INVDPAL"):
                q0, q1 = q1, q0
            ya, yb = ylev["thille"] + (0.08 if g == "westindian" else -0.08), ylev[g] + (-0.08 if g == "westindian" else 0.08)
            ax.add_patch(Polygon([(t0, ya), (t1, ya), (q1, yb), (q0, yb)], closed=True, color=KIND[kind][0],
                                 alpha=0.55 if kind == "SYNAL" else 0.8, lw=0, zorder=1))
    for g in order:
        y = ylev[g]
        ax.add_patch(Rectangle((0, y - 0.07), 1, 0.14, color="#d9d9d9", zorder=2))
        if g in block:
            b0, b1 = X(g, block[g][0]), X(g, block[g][1])
            ax.add_patch(Rectangle((b0, y - 0.07), b1 - b0, 0.14, color=BLOCK, alpha=0.85, zorder=3))
        for s, e in trf[g]:
            if e > win[g][0] and s < win[g][1]:
                ax.add_patch(Rectangle((X(g, max(s, win[g][0])), y + 0.075), X(g, min(e, win[g][1])) - X(g, max(s, win[g][0])),
                                       0.05, color=SAT, lw=0, zorder=3))
        for s, e in genes(g, chrom, *win[g]):
            ax.add_patch(Rectangle((X(g, max(s, win[g][0])), y - 0.13), max(X(g, e) - X(g, s), 0.0015), 0.05,
                                   color=GENE, lw=0, zorder=3))
        s, e = win[g]
        lo, hi = sorted((to_orig(g, chrom, s), to_orig(g, chrom, e)))
        oid = chrom if g == "thille" else IDS[g][chrom]
        rc = " (reverse-complemented)" if g != "thille" and chrom in RC else ""
        ax.text(-0.01, y, f"{NAME[g]}\n{oid}", ha="right", va="center", fontsize=8, fontweight="bold")
        ax.text(0, y - 0.2, f"{to_orig(g, chrom, s) / 1e6:.2f} Mb", ha="left", va="top", fontsize=7)
        ax.text(1, y - 0.2, f"{to_orig(g, chrom, e) / 1e6:.2f} Mb{rc}", ha="right", va="top", fontsize=7)
    ax.set_xlim(0, 1); ax.set_ylim(-0.35, 2.3); ax.axis("off")
    kinds = {}
    for col, lab in KIND.values():
        kinds.setdefault(lab, col)
    handles = [Rectangle((0, 0), 1, 1, color=col, alpha=0.7, label=lab) for lab, col in kinds.items()]
    handles += [Rectangle((0, 0), 1, 1, color=BLOCK, label="SV block"),
                Rectangle((0, 0), 1, 1, color=SAT, label="tandem array, period >= 100 bp (above bar)"),
                Rectangle((0, 0), 1, 1, color=GENE, label="gene (below bar)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=4, fontsize=7, frameon=False)
    ax.set_title(f"{a.region}: {chrom} SV locus in Unnamed Lowland, Thille and Hass UoQ (SyRI alignments, Thille window "
                 f"{win['thille'][0] / 1e6:.2f} to {win['thille'][1] / 1e6:.2f} Mb)", fontsize=9)
    for i, g in enumerate(order):
        dx = fig.add_subplot(gs[i + 1])
        xs, vs = depth(g, chrom, *win[g], a.bin, med[(a.region, g)])
        dx.fill_between([X(g, x) for x in xs], vs, step="mid", color="#7f8c8d", lw=0)
        for yv in (0.25, 0.75, 1.25):
            dx.axhline(yv, color="#9e9e9e", lw=0.5, ls="--")
        dx.axhline(1, color=INK, lw=0.7)
        if g in block:
            dx.axvspan(X(g, block[g][0]), X(g, block[g][1]), color=BLOCK, alpha=0.12)
        dx.set_xlim(0, 1); dx.set_ylim(0, 2); dx.set_xticks([])
        dx.set_ylabel("depth /\nmedian", fontsize=7)
        dx.set_title(f"Gwen CLR reads on {NAME[g]} ({a.bin // 1000} kb bins; chromosome median "
                     f"{med[(a.region, g)]:.0f}x)", fontsize=8, loc="left")
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    print(f"wrote {a.out}.pdf/.png; windows: " + ", ".join(f"{g} {win[g][0]}-{win[g][1]}" for g in order))


if __name__ == "__main__":
    main()
