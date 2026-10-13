#!/usr/bin/env python3
"""SV map of a Thille region: Hass and West Indian SyRI events, Gwen read events, read-coverage classes, genes and
TEs, all in Thille coordinates, with the flowering-gene coverage table beside the map.

Tracks, top to bottom:
  Hass SVs / WI SVs    SyRI 1.7.0 events vs Thille (sv/thille.sv.all.bed) >= --min-event bp as boxes, one row per
                       type; below them, per-bin counts of SNPs (sv/<genome>.snps.chr3_chr10.tsv.gz) and of indels
                       (SyRI DEL + INS < --min-event bp)
  Gwen                 Gwen read-coverage events >= 100 kb (sv/gwen.sv.ge100kb.thille_coords.bed), by read class
  Coverage classes     Gwen, Thille and WI reads on Thille: per bin, the fraction of positions in each class
                       (figure classes, M = int(chromosome median))
  Genes                Tiberius gene spans; flowering genes (flowering_beds/igv/flowering_genes_thille.bed) labeled
  TEs                  per-bin fraction of positions covered by RepeatMasker TEs
Shaded: the Thille centromere call (grey) and the regions passed with --highlight (yellow).

Usage (from the repo root):
    python3 scripts/sv_region_map.py --region chr3:7100000-12200000 --highlight chr3:9098588-10196806 --name chr3_deletion
    python3 scripts/sv_region_map.py --region chr3 --name chr3_whole
Writes Results_sv/<name>.pdf and .png.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coverage_class_summary import class_index, read_bed, split_name  # noqa: E402
import plot_gene_stacks as pgs  # noqa: E402

SV_BED = "sv/thille.sv.all.bed"
GWEN_BED = "sv/gwen.sv.ge100kb.thille_coords.bed"
CEN_BED = "sv/thille_centromeres.satellitefirst.bed"
GTF = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
TE_GFF = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.TE.gff"
FLOWER_BED = "flowering_beds/igv/flowering_genes_thille.bed"
FAI = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta.fai"
COVER = [("Gwen (A) reads", "gwen.clr.thille"), ("Thille (B) reads", "thille.hifi.thille"),
         ("West Indian reads", "westindian.hifi.thille")]
SV_TYPES = ["NOTAL", "HDR", "DEL", "INS", "CPL", "CPG", "DUP", "INVDP", "INV", "TRANS", "INVTR", "TDM"]
SV_COLORS = dict(zip(SV_TYPES, ["#c0392b", "#7f7f7f", "#e08a00", "#4a78b5", "#8c564b", "#17becf",
                                "#9467bd", "#bcbd22", "#2e9e44", "#d62728", "#1f77b4", "#ff7f0e"]))
GWEN_COLORS = {"absent": "#c0392b", "half": "#029e73", "two_copies": "#0173b2", "multi": "#949494"}
TABLE_COLS = [("Gwen", "Hass"), ("Thille", "Hass"), ("WI", "Hass"), ("Gwen", "Thille"), ("Thille", "Thille"),
              ("WI", "Thille"), ("Gwen", "WI"), ("Thille", "WI")]
TRACK_OF = {("Gwen", "Hass"): ("hass", "gwen.clr.hass_uoq"), ("Thille", "Hass"): ("hass", "thille.hifi.hass_uoq"),
            ("WI", "Hass"): ("hass", "westindian.hifi.hass_uoq"), ("Gwen", "Thille"): ("thille", "gwen.clr.thille"),
            ("Thille", "Thille"): ("thille", "thille.hifi.thille"), ("WI", "Thille"): ("thille", "westindian.hifi.thille"),
            ("Gwen", "WI"): ("wi", "gwen.clr.wi"), ("Thille", "WI"): ("wi", "thille.hifi.wi")}
CLASS_NAMES = ["low", "half", "medium", "high", "outlier"]


def parse_region(text):
    lengths = dict(line.split("\t")[:2] for line in open(FAI))
    m = re.fullmatch(r"(\w+)(?::(\d+)-(\d+))?", text.replace(",", ""))
    if not m:
        sys.exit(f"bad region {text}")
    chrom = m.group(1)
    if m.group(2):
        return chrom, int(m.group(2)), int(m.group(3))
    return chrom, 1, int(lengths[chrom])


def read_sv(chrom, lo, hi):
    df = pd.read_csv(SV_BED, sep="\t", dtype={"partner_chrom": str})
    df = df[(df["#chrom"] == chrom) & (df["end"] >= lo) & (df["start"] <= hi)]
    return df


def gene_table(genes_here):
    """Rows: flowering genes in the region; cells: dominant class, fraction, depth ratio for each read set x genome."""
    tabs = {k: pd.read_csv(f"summaries/flowering_gene_coverage_classes_{k}.tsv", sep="\t") for k in ("thille", "hass", "wi")}
    rows = []
    for fam, gwen_gene in genes_here:
        cells = []
        for col in TABLE_COLS:
            key, track = TRACK_OF[col]
            t = tabs[key]
            r = t[(t.gene == fam) & (t.gwen_gene == gwen_gene)]
            if r.empty or f"{track}_low" not in t.columns:
                cells.append("n/a")
                continue
            r = r.iloc[0]
            fr = [r[f"{track}_{c}"] for c in CLASS_NAMES]
            i = int(np.argmax(fr))
            cells.append(f"{CLASS_NAMES[i][:3]} {fr[i]:.2f} ({r[f'{track}_median_depth'] / r[f'{track}_chrom_median']:.2f}x)")
        rows.append([fam] + cells)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region", required=True, help="chrN or chrN:start-end (1-based, Thille)")
    p.add_argument("--highlight", nargs="*", default=[], help="regions to shade, chrN:start-end")
    p.add_argument("--name", required=True)
    p.add_argument("--min-event", type=int, help="smallest SyRI event drawn as a box (default 50 kb whole chromosome, 1 kb zoom)")
    p.add_argument("--bins", type=int, default=400, help="bins across the region for densities and class fractions")
    p.add_argument("--outdir", default="Results_sv")
    a = p.parse_args()

    chrom, lo, hi = parse_region(a.region)
    span = hi - lo + 1
    whole = lo == 1 and a.region == chrom
    min_event = a.min_event or (50_000 if whole else 1_000)
    edges = np.linspace(lo, hi + 1, a.bins + 1)
    mids = (edges[:-1] + edges[1:]) / 2
    width = edges[1] - edges[0]

    sv = read_sv(chrom, lo, hi)
    snps = {}
    for g in ("hass", "westindian"):
        pos = []
        with gzip.open(f"sv/{g}.snps.chr3_chr10.tsv.gz", "rt") as fh:
            for line in fh:
                c, x = line.split("\t")
                if c == chrom:
                    pos.append(int(x))
        snps[g] = np.histogram(pos, bins=edges)[0]
    gwen = pd.read_csv(GWEN_BED, sep="\t")
    gwen = gwen[(gwen["#chrom"] == chrom) & (gwen["end"] >= lo) & (gwen["start"] <= hi)]
    cen = [(int(f[1]), int(f[2])) for f in (l.split("\t") for l in open(CEN_BED)) if f[0] == chrom]

    cov = {}
    for label, track in COVER:
        depth, meta = pgs.load_track("cache", track, chrom)
        cls = class_index(np.asarray(depth[lo - 1:hi]), int(meta["median"]))
        idx = np.minimum(((np.arange(lo, hi + 1) - lo) / width).astype(int), a.bins - 1)
        counts = np.zeros((a.bins, 5))
        np.add.at(counts, (idx, cls), 1)
        cov[label] = counts / counts.sum(1, keepdims=True)

    genes = []
    for line in open(GTF):
        f = line.split("\t")
        if len(f) > 8 and f[0] == chrom and f[2] == "gene" and int(f[4]) >= lo and int(f[3]) <= hi:
            genes.append((int(f[3]), int(f[4])))
    flower = [(n, s, e) for n, c, s, e in read_bed(FLOWER_BED) if c == chrom and e >= lo and s <= hi]
    te_cov = np.zeros(a.bins)
    te_mask = np.zeros(span, dtype=bool)
    for line in open(TE_GFF):
        if line.startswith("#"):
            continue
        f = line.split("\t", 5)
        if f[0] == chrom and int(f[4]) >= lo and int(f[3]) <= hi:
            te_mask[max(int(f[3]), lo) - lo:min(int(f[4]), hi) - lo + 1] = True
    idx = np.minimum(((np.arange(span)) / width).astype(int), a.bins - 1)
    te_cov = np.bincount(idx, weights=te_mask, minlength=a.bins) / np.bincount(idx, minlength=a.bins)

    colors = pgs.palette()
    heights = [2.2, 0.7, 2.2, 0.7, 0.6] + [0.7] * len(COVER) + [1.0, 0.5]
    fig = plt.figure(figsize=(26, sum(heights) + 3.2))
    gs = fig.add_gridspec(len(heights), 2, width_ratios=[3.2, 2.2], height_ratios=heights, hspace=0.35, wspace=0.03)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(len(heights))]
    for ax in axes[1:]:
        ax.sharex(axes[0])
    mb = lambda x, _p: f"{x / 1e6:.2f}" if span < 20e6 else f"{x / 1e6:.0f}"

    def shade(ax):
        for s, e in cen:
            ax.axvspan(s, e, color="#bdbdbd", alpha=0.35, lw=0, zorder=0)
        for h in a.highlight:
            hc, hs, he = parse_region(h)
            if hc == chrom:
                ax.axvspan(hs, he, color="#f0e442", alpha=0.35, lw=0, zorder=0)

    row = 0
    for g, label in (("hass", "Hass vs Thille (SyRI)"), ("westindian", "West Indian vs Thille (SyRI)")):
        ax = axes[row]
        ev = sv[(sv.partner_genome == g) & (sv.size_bp >= min_event)]
        shown = [t for t in SV_TYPES if (ev.type == t).any()]
        for yi, t in enumerate(shown):
            for s, e in ev[ev.type == t][["start", "end"]].itertuples(index=False):
                ax.add_patch(Rectangle((s, yi - 0.4), max(e - s, span / 2000), 0.8, color=SV_COLORS[t], lw=0))
        ax.set_yticks(range(len(shown)), shown, fontsize=7)
        ax.set_ylim(-0.6, max(len(shown), 1) - 0.4)
        ax.set_title(f"{label}: events >= {min_event:,} bp", loc="left", fontsize=9)
        shade(ax)
        ax2 = axes[row + 1]
        small = sv[(sv.partner_genome == g) & (sv.size_bp < min_event) & sv.type.isin(["DEL", "INS"])]
        ax2.bar(mids, snps[g], width=width, color="#7f7f7f", lw=0, label="SNPs")
        ax2.bar(mids, np.histogram((small.start + small.end) / 2, bins=edges)[0], width=width, color="#e08a00", lw=0,
                label=f"indels < {min_event:,} bp")
        ax2.set_ylabel("per bin", fontsize=7)
        ax2.legend(fontsize=6, loc="upper right", frameon=False, ncol=2)
        shade(ax2)
        row += 2
    ax = axes[row]
    for s, e, det in gwen[["start", "end", "detail"]].itertuples(index=False):
        rc = re.search(r"read_class=(\w+)", str(det))
        k = rc.group(1) if rc else "multi"
        ax.add_patch(Rectangle((s, 0), e - s, 1, color=GWEN_COLORS.get(k, "#949494"), lw=0))
    ax.set_yticks([])
    ax.set_title("Gwen read events >= 100 kb (absent red, half green, two copies blue, multi-contig grey)", loc="left", fontsize=9)
    shade(ax)
    row += 1
    for label, _t in COVER:
        ax = axes[row]
        ax.stackplot(mids, cov[label].T, colors=[colors[c] for c in CLASS_NAMES], lw=0)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label.replace(" reads", ""), fontsize=7, rotation=0, ha="right", va="center")
        row += 1
    axes[row - len(COVER)].set_title("Coverage class per bin on Thille (low / half / medium / high / outlier)", loc="left", fontsize=9)
    ax = axes[row]
    for s, e in genes:
        ax.add_patch(Rectangle((s, 0.1), max(e - s, span / 3000), 0.3, color="#b8bdb5", lw=0))
    fams = []
    loci = {}
    for n, s, e in flower:  # several Gwen genes can share one Thille model (e.g. PRP8A x4): one row per locus
        loci.setdefault((s, e), n)
    for i, ((s, e), n) in enumerate(sorted(loci.items())):
        fam, gg = split_name(n)[:2]
        fams.append((fam, gg))
        ax.add_patch(Rectangle((s, 0.5), max(e - s, span / 1500), 0.4, color=pgs.PRIMARY_FILL, lw=0))
        ax.text((s + e) / 2, 0.95 + 0.35 * (i % 3), fam, fontsize=6, ha="center", va="bottom", rotation=0)
    ax.set_ylim(0, 2.1)
    ax.set_yticks([])
    ax.set_title("Genes (grey) and flowering genes (purple)", loc="left", fontsize=9)
    shade(ax)
    row += 1
    ax = axes[row]
    ax.bar(mids, te_cov, width=width, color="#8c564b", lw=0)
    ax.set_ylim(0, 1)
    ax.set_ylabel("TE", fontsize=7)
    shade(ax)
    axes[0].set_xlim(lo, hi)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    axes[-1].xaxis.set_major_formatter(plt.FuncFormatter(mb))
    axes[-1].set_xlabel(f"Thille {chrom} (Mb)")

    tax = fig.add_subplot(gs[:, 1])
    tax.axis("off")
    rows = gene_table(list(dict.fromkeys(fams)))
    if rows:
        header = ["Gene"] + [f"{r}\n→ {g}" for r, g in TABLE_COLS]
        tbl = tax.table(cellText=rows, colLabels=header, loc="upper center", cellLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6.5)
        tbl.scale(1, 1.6)
        tax.set_title("Flowering genes in the region: dominant class, fraction of the gene, gene depth / chromosome median",
                      fontsize=9, loc="left")
    else:
        tax.text(0, 0.95, "No flowering genes in this region", fontsize=9)
    legend = [Patch(color=SV_COLORS[t], label=t) for t in SV_TYPES] + \
             [Patch(color=colors[c], label=f"coverage {c}") for c in CLASS_NAMES] + \
             [Patch(color="#bdbdbd", alpha=0.35, label="centromere call"), Patch(color="#f0e442", alpha=0.35, label="highlighted region")]
    tax.legend(handles=legend, loc="lower left", ncol=3, fontsize=7, frameon=False)
    fig.suptitle(f"Structural variation on Thille {chrom}:{lo:,}-{hi:,}", fontsize=13, y=0.995)
    Path(a.outdir).mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(Path(a.outdir) / f"{a.name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    ev_counts = sv.groupby(["partner_genome", "type"]).size()
    print(f"{a.name}: {len(sv)} SyRI events in region; flowering genes {len(fams)}", file=sys.stderr)
    print(ev_counts.to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
