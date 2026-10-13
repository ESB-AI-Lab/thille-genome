#!/usr/bin/env python3
"""Genome-wide read coverage on the Thille assembly: supplementary table and figure (no UoQ reference, Edwin 2026-10-02).

Table: for Gwen CLR, Thille HiFi and Lowland (WI) HiFi reads on Thille, per chromosome and genome-wide, the chromosome
median depth and the bp and % of positions in each coverage class (low <= 0.25 M, half <= 0.75 M, medium <= 1.25 M,
high <= 1.5 M, outlier > 1.5 M; M = chromosome median), from summaries/coverage_class_counts_thille.tsv.
Figure: the 12 chromosomes end to end; for each read set, the median depth of non-overlapping 100 kb bins divided by
the chromosome median (per-position depth from cache/<track>.<chr>.depth.npy), points colored by the class of the
bin value (display only; bins are not used for classification) and drawn at 2.0 when above it.

Usage (from the repo root): python3 scripts/genome_wide_coverage.py [--bin 100000] [--no-text]
Writes summaries/genome_wide_coverage_thille.tsv and Results_genome_coverage/genome_wide_coverage_thille.{pdf,png}.
With --no-text (Edwin, 2026-10-09: no explanatory text inside figures) only genome_wide_coverage_thille_notext.{pdf,png},
without the footnote (it goes to the legend md); the TSV is not rewritten.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

TRACKS = [("gwen.clr.thille", "Gwen CLR"), ("thille.hifi.thille", "Thille HiFi"), ("westindian.hifi.thille", "Unnamed Lowland HiFi")]
CLASSES = ["low", "half", "medium", "high", "outlier"]
COLORS = {"low": "#E69F00", "half": "#009E73", "medium": "#0072B2", "high": "#CC79A7", "outlier": "#D55E00"}
CHROMS = [f"chr{i}" for i in range(1, 13)]


def cls(r):
    return np.select([r <= 0.25, r <= 0.75, r <= 1.25, r <= 1.5], CLASSES[:4], "outlier")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bin", type=int, default=100000)
    p.add_argument("--no-text", action="store_true", help="no footnote; writes *_notext only (default off; output "
                   "without it is unchanged)")
    a = p.parse_args()
    rows = [r for r in csv.DictReader(open("summaries/coverage_class_counts_thille.tsv"), delimiter="\t")]
    out = []
    for track, label in TRACKS:
        sub = [r for r in rows if r["alignment"] == track and r["chromosome"] != "genome"]  # skip the summary script's total row
        if len(sub) != 12:
            raise SystemExit(f"{track}: {len(sub)} chromosomes in coverage_class_counts_thille.tsv, expected 12")
        tot = {c: 0 for c in CLASSES}; tl = 0
        for r in sorted(sub, key=lambda r: int(r["chromosome"][3:])):
            L = int(r["length"]); tl += L
            rec = {"read_set": label, "chromosome": r["chromosome"], "median_depth": r["median"], "length_bp": L}
            for c in CLASSES:
                rec[f"{c}_bp"] = int(r[c]); rec[f"{c}_pct"] = f"{100 * int(r[c]) / L:.1f}"; tot[c] += int(r[c])
            out.append(rec)
        meds = [int(r["median"]) for r in sub]
        rec = {"read_set": label, "chromosome": "genome", "median_depth": f"{min(meds)} to {max(meds)}", "length_bp": tl}
        for c in CLASSES:
            rec[f"{c}_bp"] = tot[c]; rec[f"{c}_pct"] = f"{100 * tot[c] / tl:.1f}"
        out.append(rec)
    if not a.no_text:
        with open("summaries/genome_wide_coverage_thille.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
    for r in out:
        if r["chromosome"] == "genome":
            print(r["read_set"], "medians", r["median_depth"], " ".join(f"{c} {r[c + '_pct']}%" for c in CLASSES))

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7})
    fig, axes = plt.subplots(len(TRACKS), 1, figsize=(7.2, 4.6), sharex=True)
    offsets, x0 = {}, 0
    for c in CHROMS:
        offsets[c] = x0
        x0 += json.load(open(f"cache/thille.hifi.thille.{c}.json"))["length"]
    for ax, (track, label) in zip(axes, TRACKS):
        for c in CHROMS:
            meta = json.load(open(f"cache/{track}.{c}.json"))
            dep = np.load(f"cache/{track}.{c}.depth.npy", mmap_mode="r")
            n = len(dep) // a.bin
            med = np.median(np.asarray(dep[:n * a.bin]).reshape(n, a.bin), axis=1) / meta["median"]
            xs = offsets[c] + (np.arange(n) + 0.5) * a.bin
            k = cls(med)
            for name in CLASSES:
                m = k == name
                ax.scatter(xs[m] / 1e6, np.minimum(med[m], 2.0), s=1.2, color=COLORS[name], lw=0, rasterized=True)
        for y in (0.25, 0.75, 1.25, 1.5):
            ax.axhline(y, color="#9aa0a6", lw=0.4, ls="--")
        ax.axhline(1.0, color="#333333", lw=0.5)
        for c in CHROMS[1:]:
            ax.axvline(offsets[c] / 1e6, color="#d0d0d0", lw=0.4)
        ax.set_ylim(-0.05, 2.1); ax.set_ylabel(f"{label}\ndepth / median", fontsize=6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks([(offsets[c] + json.load(open(f"cache/thille.hifi.thille.{c}.json"))["length"] / 2) / 1e6
                         for c in CHROMS], [c.replace("chr", "") for c in CHROMS])
    axes[-1].set_xlabel("Thille chromosome")
    handles = [plt.Line2D([], [], marker="o", ls="", color=COLORS[c], label=c, markersize=4) for c in CLASSES]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=6, bbox_to_anchor=(0.5, -0.04))
    if not a.no_text:
        fig.text(0.01, -0.09, f"Median depth of {a.bin // 1000} kb bins divided by the chromosome median, Thille reference; "
                 "dashed lines: class boundaries 0.25, 0.75, 1.25, 1.5; values above 2.0 drawn at 2.0.", fontsize=5.5,
                 color="#555555")
    Path("Results_genome_coverage").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"Results_genome_coverage/genome_wide_coverage_thille{'_notext' if a.no_text else ''}.{ext}", dpi=600,
                    bbox_inches="tight")


if __name__ == "__main__":
    main()
