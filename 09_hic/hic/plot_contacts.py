#!/usr/bin/env python3
"""Hi-C contact map figures from a balanced multi-resolution cooler (.mcool).

For each scope it writes (PNG 300 dpi + PDF):
  <scope>_genome_<wg>kb.{png,pdf}          whole genome, log10 raw counts | log10 ICE-balanced
  <scope>_perchrom_<pc>kb_{raw,balanced}  3x4 grid of chromosomes (only named chrN sequences)

Color: single-hue sequential (white -> dark red); zeros/NaN (unbalanced bins) are left white. Limits are
computed per panel from that panel's own values (--scale):
  pct     log, 1st-99.5th percentile of non-zero values (original default)
  log     log, floor = median of same-sequence pixels >= 10 Mb apart, ceiling = 95th percentile of
          same-sequence pixels within 1 Mb of the diagonal (background drops to white)
  linear  linear, 0 to the 99.5th percentile of non-zero values (Juicebox/PretextView-like saturation)
In the "all" scope only chr* and scaffold_13 get boundary lines/labels; the remaining small scaffolds
are bracketed as one "unplaced" block.

Usage: plot_contacts.py <scope.mcool> <scope> <outdir> <pairstats.tsv> [--wg-res 500000] [--chrom-res 100000]
"""
import argparse
import csv
import os

import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize  # noqa: E402

CMAP = LinearSegmentedColormap.from_list("hic_red", ["#ffffff", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"])
CMAP.set_bad("#ffffff")
INK, MUTED, RULE = "#1f1f1f", "#5f5f5f", "#9a9a9a"
plt.rcParams.update({"font.size": 8, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": RULE,
                     "pdf.fonttype": 42})


SCALE_TEXT = {"pct": "log color, 1st–99.5th percentile of non-zero pixels",
              "log": "log color, floor = median far-cis (≥10 Mb), ceiling = 95th pct near-diagonal (≤1 Mb)",
              "linear": "linear color, 0 to 99.5th percentile of non-zero pixels"}


def norm_for(mat, scale="pct", seq_ids=None, res=None):
    ok = np.isfinite(mat) & (mat > 0)
    v = mat[ok]
    if v.size == 0:
        raise ValueError("matrix has no positive finite values")
    if scale == "pct":
        lo, hi = np.percentile(v, [1, 99.5])
        return LogNorm(vmin=lo, vmax=max(hi, lo * 10))
    if scale == "linear":
        return Normalize(vmin=0, vmax=np.percentile(v, 99.5))
    if scale == "log":
        n = mat.shape[0]
        ids = np.zeros(n, dtype=int) if seq_ids is None else np.asarray(seq_ids)
        same = ids[:, None] == ids[None, :]
        d = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        far = mat[ok & same & (d >= int(np.ceil(10e6 / res)))]
        near = mat[ok & same & (d <= max(1, int(1e6 // res)))]
        if far.size == 0 or near.size == 0:
            raise ValueError("log scale needs same-sequence pixels both >= 10 Mb apart and within 1 Mb")
        lo, hi = float(np.median(far)), float(np.percentile(near, 95))
        return LogNorm(vmin=lo, vmax=max(hi, lo * 10))
    raise ValueError(f"unknown scale {scale}")


def value_label(raw, scale):
    what = "raw contacts" if raw else "ICE-balanced"
    return f"{what} (linear)" if scale == "linear" else f"log10 {what}"


def labelled(names):
    return [n for n in names if n.startswith("chr") or n == "scaffold_13"]


NO_TEXT = False


def genome_figure(clr, scope, outdir, res_kb, weight="weight", raw_only=False, label="", suffix="", scale="pct"):
    names = list(clr.chromnames)
    raw = clr.matrix(balance=False)[:].astype(float)
    bal = None if raw_only else clr.matrix(balance=weight)[:]
    n = raw.shape[0]
    marks = [(name, *clr.extent(name)) for name in labelled(names)]
    unplaced = [nm for nm in names if nm not in {m[0] for m in marks}]

    seq_ids = np.concatenate([[i] * len(clr.bins().fetch(nm)) for i, nm in enumerate(names)])
    panels = [(raw, value_label(True, scale))]
    if not raw_only:
        panels.append((bal, value_label(False, scale) + (f" ({weight})" if weight != "weight" else "")))
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2 * len(panels) - 0.4, 7.4), constrained_layout=True, squeeze=False)
    for ax, (mat, title) in zip(axes[0], panels):
        im = ax.imshow(np.ma.masked_invalid(np.where(mat > 0, mat, np.nan)), cmap=CMAP, norm=norm_for(mat, scale, seq_ids, clr.binsize),
                       interpolation="none", origin="upper", extent=(0, n, n, 0))
        for _, s, e in marks:
            for b in (s, e):
                ax.axhline(b, color=RULE, lw=0.4)
                ax.axvline(b, color=RULE, lw=0.4)
        ticks = [(s + e) / 2 for _, s, e in marks]
        labels = [nm.replace("scaffold_", "scf") for nm, _, _ in marks]
        if unplaced:
            s = clr.extent(unplaced[0])[0]
            ticks.append((s + n) / 2)
            labels.append(f"unplaced\n({len(unplaced)})")
        ax.set_xticks(ticks, labels, rotation=90)
        ax.set_yticks(ticks, labels)
        ax.tick_params(length=0)
        if NO_TEXT:
            ax.set_title("AB"[list(axes[0]).index(ax)] if len(panels) > 1 else "", loc="left", fontweight="bold", color=INK)
        else:
            ax.set_title(title, color=INK)
        cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.01)
        cb.set_label("contacts per bin" if mat is raw else "balanced contact frequency", color=MUTED)
        cb.outline.set_edgecolor(RULE)
    what = {"chr12": "chr1–12", "chr12scf13": "chr1–12 + scaffold_13"}.get(
        scope, f"all {len(names)} scaffolds (chr1–12, then by size)")
    if not NO_TEXT:
      fig.suptitle(f"Thille Hi-C contact map — {what}, {res_kb} kb bins; {SCALE_TEXT[scale]}"
                 + (f"\n{label}" if label else ""), color=INK)
    base = os.path.join(outdir, f"{scope}_genome_{res_kb}kb{suffix}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", dpi=300)
    plt.close(fig)
    return base


def perchrom_figure(clr, scope, outdir, res_kb, balance, scale="pct", suffix=""):
    chroms = [c for c in clr.chromnames if c.startswith("chr")]
    rows, cols = 3, 4
    fig, axes = plt.subplots(rows, cols, figsize=(13, 10), constrained_layout=True)
    for ax in axes.flat[len(chroms):]:
        ax.set_visible(False)
    for ax, ch in zip(axes.flat, chroms):
        mat = clr.matrix(balance=balance).fetch(ch)
        mat = mat.astype(float) if not balance else mat
        length_mb = clr.chromsizes[ch] / 1e6
        im = ax.imshow(np.ma.masked_invalid(np.where(mat > 0, mat, np.nan)), cmap=CMAP, norm=norm_for(mat, scale, None, clr.binsize),
                       interpolation="none", extent=(0, length_mb, length_mb, 0))
        ax.set_title(ch, color=INK)
        ax.tick_params(length=2)
        fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02).outline.set_edgecolor(RULE)
    kind = "balanced" if balance else "raw"
    if not NO_TEXT:
      fig.suptitle(f"Thille Hi-C per chromosome — {value_label(not balance, scale)}, {res_kb} kb bins "
                 f"(axes in Mb; limits per panel: {SCALE_TEXT[scale]})", color=INK)
    base = os.path.join(outdir, f"{scope}_perchrom_{res_kb}kb_{kind}{suffix}")
    for ext in ("png", "pdf"):
        fig.savefig(f"{base}.{ext}", dpi=300)
    plt.close(fig)
    return base


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mcool")
    ap.add_argument("scope")
    ap.add_argument("outdir")
    ap.add_argument("pairstats")
    ap.add_argument("--wg-res", type=int, default=500000, help="whole-genome bin size (default 500000)")
    ap.add_argument("--chrom-res", type=int, default=100000, help="per-chromosome bin size (default 100000)")
    ap.add_argument("--weight-name", default="weight", help="balancing weight column to plot (default weight)")
    ap.add_argument("--raw-only", action="store_true", help="plot only the raw panel")
    ap.add_argument("--variant-label", default="", help="extra title line naming the balancing variant")
    ap.add_argument("--suffix", default="", help="filename tag appended to the genome figure name")
    ap.add_argument("--genome-only", action="store_true", help="skip per-chromosome panels")
    ap.add_argument("--scale", choices=sorted(SCALE_TEXT), default="pct", help="color scaling (default pct)")
    ap.add_argument("--no-text", action="store_true", help="omit titles; genome panels get only the letters A/B (text goes in the figure legend); default off")
    a = ap.parse_args()
    global NO_TEXT
    NO_TEXT = a.no_text
    os.makedirs(a.outdir, exist_ok=True)

    wg = cooler.Cooler(f"{a.mcool}::resolutions/{a.wg_res}")
    pc = cooler.Cooler(f"{a.mcool}::resolutions/{a.chrom_res}")
    with open(a.pairstats) as fh:
        stats = next(csv.DictReader(fh, delimiter="\t"))
    pixel_sum = int(wg.pixels()[:]["count"].sum())
    if pixel_sum != int(stats["kept"]):
        raise SystemExit(f"ERROR: cooler count sum {pixel_sum} != kept pairs {stats['kept']}")
    print(f"pairs kept {stats['kept']} == cooler sum {pixel_sum}; cis fraction {stats['cis_frac_of_kept']}")

    print("wrote", genome_figure(wg, a.scope, a.outdir, a.wg_res // 1000, a.weight_name, a.raw_only,
                                 a.variant_label, a.suffix, a.scale))
    # Per-chromosome panels only for the chr1-12 scope; in "all" the raw intra-chromosome blocks are identical.
    if a.scope == "chr12" and not a.genome_only:
        for bal in (False, True):
            print("wrote", perchrom_figure(pc, a.scope, a.outdir, a.chrom_res // 1000, bal, a.scale, a.suffix))


if __name__ == "__main__":
    main()
