#!/usr/bin/env python3
"""Read-level breakpoint figures for the SyRI v4 inversions (one figure per inversion).

Columns: the inversion's 100 kb window(s) on Thille, then on the compared genome (Hass UoQ, Unnamed Lowland or Puebla).
Rows, top to bottom: one row per read set aligned to that genome (depth / chromosome median: filled area = all reads,
line = reads with MAPQ >= 20; ticks at the top = reads clipped >= 500 bp within 1 kb of a breakpoint, filled when the
clipped part aligns, on the opposite strand, near the other breakpoint of the inversion); flowering genes; genes; TEs.
The inversion is the shaded span, its breakpoints the dashed lines. Each read row reports, for the breakpoint(s) in
the window, reads spanning the breakpoint by >= 5 kb on both sides and split reads ("span" / "split"; start, end).
Inputs are the outputs of make_windows.py, inv_pileup.py and te/sbatch_repeatmasker_invwin.sh. Colors: the dataviz
reference palette (read sets = categorical slots 1, 2, 3 and 7, validated adjacent; flowering = slot 8; TEs green and
yellow, validated).

Usage (inside the asmbenchmarkers container, from assemblies/thille.hic):
  plot_inversions.py [--only INV_ID ...] [--dpi 200]
Writes compare/v4/inversions/<outdir>/<pair>/<inv_id>.{pdf,png} (default figures_v2; figures/ holds the three-read-set
version drawn before the Puebla alignments existed).
"""
import argparse
import collections
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = "compare/v4/inversions"
INK, INK2, MUTED, GRID, BASE, SURF, SPAN = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#f0efec"
READSETS = [("thille_hifi", "Thille HiFi", "#2a78d6"), ("lowland_hifi", "Unnamed Lowland HiFi", "#eb6834"),
            ("gwen_clr", "Gwen CLR", "#1baf7a"), ("puebla_hifi", "Puebla HiFi", "#4a3aa7")]
FLOWER = "#e34948"
TE_COL = {"LTR": "#008300", "non-LTR / DNA": "#eda100", "Unknown": BASE, "Simple / low complexity": INK2}
NAME = {"thille": "Thille", "hass": "Hass UoQ", "westindian": "Unnamed Lowland", "puebla": "Puebla"}
LIB = "/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/consensi.fa.classified"
FL = "/gpool/edwin/avocados/thille_paper/flowering_beds/v4/igv"
FLOWER_BED = {"thille": f"{FL}/flowering_genes_thille.bed", "hass": f"{FL}/flowering_genes_hass.bed",
              "westindian": f"{FL}/flowering_genes_wi.bed"}
GENES = {"thille": "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf",
         "hass": "/epool/edwin/pamericana/avoflow/z_GCA_029852735.1/ncbi_dataset/data/GCA_029852735.1/"
                 "GCA_029852735.1_ASM2985273v1_genes.bed",
         "westindian": "westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/genomic.gff"}
THILLE_TE = "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "pdf.fonttype": 42, "axes.edgecolor": BASE,
                     "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
                     "ytick.labelcolor": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.major.size": 2.5, "ytick.major.size": 2.5, "text.color": INK})


def te_group(cls):
    c = cls or "Unknown"
    if c.startswith("LTR"):
        return "LTR"
    if c.startswith(("LINE", "SINE", "DNA", "RC")):
        return "non-LTR / DNA"
    if c.startswith(("Simple", "Low_complexity", "Satellite", "rRNA", "snRNA", "tRNA", "srpRNA", "scRNA")):
        return "Simple / low complexity"
    return "Unknown"


def family_classes():
    cls = {}
    for line in open(LIB):
        if line.startswith(">"):
            fam, _, c = line[1:].split()[0].partition("#")
            cls[fam] = c
    return cls


def motif_group(attr, fams):
    m = attr.split('Motif:', 1)[1].split('"', 1)[0] if "Motif:" in attr else ""
    if m.startswith("(") or m.endswith("-rich"):
        return "Simple / low complexity"
    return te_group(fams.get(m))


def load_intervals(windows):
    """Annotation intervals for every window, keyed by (genome, chrom, win_start)."""
    want = collections.defaultdict(list)
    for w in windows:
        want[(w["genome"], w["chrom"])].append((w["ws"], w["we"]))

    def hit(g, c, s, e):
        return [ws for ws, we in want.get((g, c), []) if s < we and e > ws]

    flow, genes, tes = (collections.defaultdict(list) for _ in range(3))
    for g, path in FLOWER_BED.items():
        for line in open(path):
            if line.startswith("track"):
                continue
            f = line.rstrip("\n").split("\t")
            for ws in hit(g, f[0], int(f[1]), int(f[2])):
                flow[(g, f[0], ws)].append((int(f[1]), int(f[2]), f[3].split("|")[0]))
    for line in open(GENES["thille"]):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "gene":
            for ws in hit("thille", f[0], int(f[3]) - 1, int(f[4])):
                genes[("thille", f[0], ws)].append((int(f[3]) - 1, int(f[4]), f[6]))
    for line in open(GENES["hass"]):
        f = line.rstrip("\n").split("\t")
        for ws in hit("hass", f[0], int(f[1]), int(f[2])):
            genes[("hass", f[0], ws)].append((int(f[1]), int(f[2]), f[5]))
    for line in open(GENES["westindian"]):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "gene":
            for ws in hit("westindian", f[0], int(f[3]) - 1, int(f[4])):
                genes[("westindian", f[0], ws)].append((int(f[3]) - 1, int(f[4]), f[6]))
    fams = family_classes()
    for line in open(THILLE_TE):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        for ws in hit("thille", f[0], int(f[3]) - 1, int(f[4])):
            tes[("thille", f[0], ws)].append((int(f[3]) - 1, int(f[4]), motif_group(f[8], fams)))
    for path in glob.glob(f"{D}/te/windows_part*/windows_part*.fa.out.gff"):
        for line in open(path):
            if line.startswith("#"):
                continue
            f = line.split("\t")
            g, c, off = f[0].split("__")
            off = int(off)
            tes[(g, c, off)].append((off + int(f[3]) - 1, off + int(f[4]), motif_group(f[8], fams)))
    return flow, genes, tes


def read_bins(path):
    rows = [l.rstrip("\n").split("\t") for l in list(open(path))[1:]]
    return ([(int(r[0]) + int(r[1])) / 2 for r in rows], [float(r[2]) for r in rows], [float(r[3]) for r in rows])


def read_clips(path):
    return [(int(l.split("\t")[0]), l.split("\t")[1].strip()) for l in list(open(path))[1:]]


def lanes(iv):
    """Greedy lane assignment so overlapping intervals do not draw on top of each other."""
    ends, out = [], []
    for s, e, *rest in sorted(iv):
        for i, last in enumerate(ends):
            if s > last:
                ends[i] = e
                out.append((i, s, e, *rest))
                break
        else:
            ends.append(e)
            out.append((len(ends) - 1, s, e, *rest))
    return out, max(1, len(ends))


def style(ax, ws, we, show_x=False):
    ax.set_xlim(ws, we)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if show_x:
        ticks = [ws, ws + 50000, we]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t / 1e6:.3f}" for t in ticks])
        for lab, ha in zip(ax.get_xticklabels(), ("left", "center", "right")):
            lab.set_horizontalalignment(ha)
    else:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)


def shade(ax, w):
    a, b = max(w["s"], w["ws"]), min(w["e"], w["we"])
    if b > a:
        ax.axvspan(a, b, color=SPAN, lw=0, zorder=0)
    for x in (w["s"], w["e"]):
        if w["ws"] <= x <= w["we"]:
            ax.axvline(x, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--no-text", action="store_true", help="omit coordinates, notes and read counts drawn in the panels (they go in the legend and inversion tables); keep only the inversion id; default off")
    ap.add_argument("--outdir", default="figures_v2", help="output folder under compare/v4/inversions/")
    a = ap.parse_args()
    inv = {r[0]: r for r in (l.rstrip("\n").split("\t") for l in list(open(f"{D}/inversions.tsv"))[1:])}
    windows = []
    for l in list(open(f"{D}/windows.tsv"))[1:]:
        f = l.rstrip("\n").split("\t")
        windows.append(dict(inv=f[0], pair=f[1], genome=f[2], side=f[3], chrom=f[4], ws=int(f[5]), we=int(f[6]),
                            s=int(f[7]), e=int(f[8])))
    if a.only:
        windows = [w for w in windows if w["inv"] in a.only]
    sup = collections.defaultdict(dict)
    medians = {}
    for p in glob.glob(f"{D}/support/task*.tsv"):
        for l in list(open(p))[1:]:
            f = l.rstrip("\n").split("\t")
            sup[(f[0], f[2], f[3], f[7])][f[5]] = (int(f[9]), int(f[11]))
            medians[(f[2], f[4], f[7])] = float(f[8])  # chromosome median recorded with each support row
    flow, genes, tes = load_intervals(windows)
    by_inv = collections.OrderedDict()
    for w in windows:
        by_inv.setdefault(w["inv"], []).append(w)
    nr = len(READSETS)
    rows = ["reads"] * nr + ["flowering", "genes", "tes"]
    hr = [1.0] * nr + [0.42, 0.55, 0.42]
    for iid, ws_ in by_inv.items():
        r = inv[iid]
        pair = r[1]
        cols = sorted(ws_, key=lambda w: (w["genome"] != "thille", {"whole": 0, "left": 1, "right": 2}[w["side"]]))
        W, H = 1.55 + 2.7 * len(cols), 6.5
        fig, axes = plt.subplots(len(rows), len(cols), figsize=(W, H),
                                 gridspec_kw=dict(height_ratios=hr, hspace=0.22, wspace=0.1), squeeze=False)
        fig.patch.set_facecolor("white")
        for j, w in enumerate(cols):
            g = w["genome"]
            key = (g, w["chrom"], w["ws"])
            size = w["e"] - w["s"]
            sub = {"whole": f"inversion {w['s'] + 1:,}-{w['e']:,} ({size:,} bp)",
                   "left": f"start breakpoint at {w['s'] + 1:,} (inversion {size / 1e6:.2f} Mb)",
                   "right": f"end breakpoint at {w['e']:,} (inversion {size / 1e6:.2f} Mb)"}[w["side"]]
            axes[0][j].set_title(f"{NAME[g]}  {w['chrom']}" if a.no_text else f"{NAME[g]}  {w['chrom']}\n{sub}",
                                 fontsize=7, color=INK, loc="left", pad=4)
            for i, (rs, label, col) in enumerate(READSETS):
                ax = axes[i][j]
                style(ax, w["ws"], w["we"])
                bins = f"{D}/bins/{iid}_{g}_{w['side']}.{rs}.bins.tsv"
                if not os.path.exists(bins):
                    ax.set_ylim(0, 2)
                    ax.set_yticks([])
                    ax.spines["left"].set_visible(False)
                    if i == nr // 2:
                        ax.text(0.5, 0.5, "no reads" if a.no_text else f"no read alignments to {NAME[g]}", transform=ax.transAxes, ha="center",
                                va="center", color=MUTED, fontsize=7)
                    continue
                med = medians[(g, w["chrom"], rs)]
                x, dep, depq = read_bins(bins)
                med = med or 1.0
                shade(ax, w)
                ax.axhline(1.0, color=GRID, lw=0.6, zorder=1)
                ax.fill_between(x, [min(d / med, 2.0) for d in dep], color=col, alpha=0.28, lw=0, zorder=2)
                ax.plot(x, [min(d / med, 2.0) for d in depq], color=col, lw=0.9, zorder=3)
                for pos, kind in read_clips(f"{D}/clips/{iid}_{g}_{w['side']}.{rs}.clips.tsv"):
                    ax.plot(pos, 2.18, marker="v", ms=3.2, mew=0.5, color=INK if kind == "split" else MUTED,
                            mfc=INK if kind == "split" else "none", zorder=4)
                ax.set_ylim(0, 2.3)
                ax.set_yticks([0, 1, 2])
                if j:
                    ax.set_yticklabels([])
                s = sup.get((iid, g, w["side"], rs), {})
                txt = "  ".join(f"{k}: span {s[k][0]}, split {s[k][1]}" for k in ("start", "end") if k in s)
                if txt and not a.no_text:
                    ax.text(0.01, 0.97, txt, transform=ax.transAxes, ha="left", va="top", fontsize=5.6,
                            color=INK2, bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85),
                            zorder=5)
                if j == 0:
                    ax.set_ylabel(label, fontsize=6.5, color=INK, rotation=0, ha="right",
                                  va="center", labelpad=6)
            # flowering genes
            ax = axes[nr][j]
            style(ax, w["ws"], w["we"])
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
            shade(ax, w)
            fl = flow.get(key, [])
            if g not in FLOWER_BED:
                ax.text(0.5, 0.5, "no annotation", transform=ax.transAxes, ha="center", va="center", color=MUTED,
                        fontsize=6)
            for s_, e_, nm in fl:
                ax.add_patch(plt.Rectangle((s_, 0.15), e_ - s_, 0.4, color=FLOWER, lw=0, zorder=3))
                ax.text(max(s_, w["ws"]) + 300, 0.62, nm, fontsize=6, color=INK, va="bottom", zorder=4)
            # genes
            ax = axes[nr + 1][j]
            style(ax, w["ws"], w["we"])
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
            shade(ax, w)
            if g not in GENES:
                ax.set_ylim(0, 1)
                ax.text(0.5, 0.5, "no gene annotation", transform=ax.transAxes, ha="center", va="center",
                        color=MUTED, fontsize=6)
            else:
                lv, n = lanes(genes.get(key, []))
                for lane, s_, e_, strand in lv:
                    ax.add_patch(plt.Rectangle((s_, lane + 0.2), e_ - s_, 0.6, color=INK2, lw=0, zorder=3))
                ax.set_ylim(n + 0.1, -0.1)
            # TEs
            ax = axes[nr + 2][j]
            style(ax, w["ws"], w["we"], show_x=True)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
            shade(ax, w)
            for s_, e_, grp in tes.get(key, []):
                ax.add_patch(plt.Rectangle((s_, 0.15), e_ - s_, 0.7, color=TE_COL[grp], lw=0, zorder=3))
            ax.set_xlabel(f"{w['chrom']} position (Mb)", fontsize=6.5, color=INK2, labelpad=2)
            if j == 0:
                for i, lab in ((nr, "flowering genes"), (nr + 1, "genes"), (nr + 2, "TEs")):
                    axes[i][0].set_ylabel(lab, fontsize=6.5, color=INK, rotation=0, ha="right", va="center",
                                          labelpad=6)
        tl, ol = int(r[5]), int(r[9])
        fig.text(0.01, 1 - 0.12 / H, iid if a.no_text else f"{iid.split('_', 1)[1]}: Thille vs {NAME[pair]}", ha="left", va="top",
                 fontsize=9, fontweight="bold", color=INK)
        if not a.no_text:
          fig.text(0.01, 1 - 0.33 / H, f"Thille {r[2]}:{int(r[3]) + 1:,}-{int(r[4]):,} ({tl:,} bp);  {NAME[pair]} "
                 f"{r[6]}:{int(r[7]) + 1:,}-{int(r[8]):,} ({ol:,} bp).\nDepth = read depth / chromosome median; 100 kb "
                 f"windows.", ha="left", va="top", fontsize=6.5, color=INK2, linespacing=1.4)
        handles = [Patch(color=c, alpha=0.5, label=l) for _, l, c in READSETS]
        handles += [Line2D([], [], color=INK2, lw=0.9, label="line: reads with MAPQ >= 20; fill: all reads"),
                    Line2D([], [], marker="v", ls="", color=INK, ms=4, label="split read to partner breakpoint"),
                    Line2D([], [], marker="v", ls="", color=MUTED, mfc="none", ms=4, label="clipped read"),
                    Patch(color=SPAN, label="inversion"), Patch(color=FLOWER, label="flowering gene")]
        handles += [Patch(color=c, label=f"TE: {k}") for k, c in TE_COL.items()]
        fig.legend(handles=handles, loc="lower left", ncol=4 if len(cols) <= 2 else 7, fontsize=5.8, frameon=False,
                   bbox_to_anchor=(0.01, 0.0), handlelength=1.4, columnspacing=1.0)
        fig.subplots_adjust(left=1.3 / W, right=1 - 0.1 / W, top=1 - 1.1 / H,
                            bottom=(0.95 if len(cols) <= 2 else 0.7) / H)
        out = f"{D}/{a.outdir}/{pair}"
        os.makedirs(out, exist_ok=True)
        for ext in ("pdf", "png"):
            fig.savefig(f"{out}/{iid}.{ext}", dpi=a.dpi, facecolor="white")
        plt.close(fig)
        print(f"wrote {out}/{iid}.pdf/.png", flush=True)


if __name__ == "__main__":
    main()
