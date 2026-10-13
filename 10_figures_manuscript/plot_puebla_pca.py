#!/usr/bin/env python3
"""Supplementary Figure S16: Puebla projected onto the Gwen-study PCA, PC1 vs PC2 only (Edwin, 2026-10-11).

Points are the 30 panel accessions at their published coordinates and the samples projected by the thille-hic-bf
session (assemblies/thille.hic/puebla_pca/<proj>/projection.tsv: PLINK 2 --score with the panel allele weights,
variance-standardized, calibrated per PC on the panel): Puebla HiFi, and Thille HiFi as a control that lands on its
own panel point. Panel accessions are colored by the population-genetic group of Solares et al. 2023, Table 1
(Guatemalan, Mexican, Lowland; accessions without a group are hybrids). No title or notes in the figure (they are in manuscript/figure_legends/Figure_S16.md).

Usage (from the repo root):
  python3 scripts/plot_puebla_pca.py [--projection <projection.tsv>] [-o Results_v4_notext/figureS16/puebla_pca_pc1_pc2]
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PROJ = "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/puebla_pca/proj_contig_noOG_v3/projection.tsv"
# Solares et al. 2023 (G3 13:jkac323), Table 1: groups G, M, L only (Edwin, 2026-10-11). The three wild Chiapas
# accessions take their historical ecotype (G). Accessions without a PopGen group are assigned from the K = 3
# admixture of Table 1: at least 80% of one group puts them in that group (Gwen and Mendez 100% G, Pinkerton 81% G;
# Edwin, 2026-10-11), otherwise they are hybrids (069-02, Bacon, Fuerte, Zutano, Velvick, CH-CR-25). The 80% cut
# matches the paper's own G calls (Lyon 84%, Anaheim 93%).
GROUP = {
    "Guatemalan": ["Anaheim", "Carlsbad", "Hass-R", "Linda", "Lyon", "Nabal", "Nimlioh", "Reed", "Taft", "Thille",
                   "CH-G-07", "CH-G-10", "CH-G-11", "Gwen", "Mendez", "Pinkerton"],
    "Mexican": ["Ganter", "PequenoCharly", "TopaTopa"],
    "Lowland": ["263-C", "Fairchild", "Simmonds", "VC26", "Waldin"],
    "Hybrid": ["CH-CR-25", "069-02", "Bacon", "Fuerte_1", "Velvick", "Zutano"],
}
COLOR = {"Guatemalan": "#009E73", "Mexican": "#0072B2", "Lowland": "#D55E00",   # Edwin: G green, M blue, L orange
         "Hybrid": "#999999"}
LABEL = {"Fuerte_1": "Fuerte", "Hass-R": "Hass, Mendez", "Mendez": "", "PequenoCharly": "Pequeno Charly"}
# no leader lines (Edwin); Hass and Mendez (a Hass sport) sit on one point and share a label
COLUMN = {}
# label offsets (points) for crowded accessions, set by eye on the PC1-PC2 plot; others go up and right
OFFSET = {"Taft": (4, 3), "Nabal": (-30, 3), "Linda": (-28, -2), "Nimlioh": (-38, -6), "Carlsbad": (-8, -11),
          "Reed": (4, -2), "Anaheim": (7, -13), "CH-G-07": (-44, -3), "CH-G-10": (-16, -15), "CH-G-11": (-36, 2),
          "Lyon": (5, -4), "Gwen": (-4, 7), "Hass-R": (5, 4), "Thille": (7, -7),
          "069-02": (-14, -10), "Zutano": (4, 3), "Bacon": (4, -9), "Fairchild": (-40, -3), "263-C": (4, -4),
          "Simmonds": (4, 2), "Waldin": (4, -8)}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--projection", default=PROJ)
    p.add_argument("-o", "--out", default="Results_v4_notext/figureS16/puebla_pca_pc1_pc2")
    a = p.parse_args()
    rows = list(csv.DictReader(open(a.projection), delimiter="\t"))
    group_of = {s: g for g, ss in GROUP.items() for s in ss}
    panel = [r for r in rows if r["source"].startswith("panel")]
    missing = [r["sample"] for r in panel if r["sample"] not in group_of]
    if missing:
        raise SystemExit(f"no group for panel samples {missing}")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for g in GROUP:
        pts = [r for r in panel if group_of[r["sample"]] == g]
        ax.scatter([float(r["PC1"]) for r in pts], [float(r["PC2"]) for r in pts], s=36, c=COLOR[g], label=g,
                   edgecolors="white", linewidths=0.5, zorder=2)
        for r in pts:
            xy = (float(r["PC1"]), float(r["PC2"]))
            if r["sample"] in COLUMN:
                ax.annotate(LABEL.get(r["sample"], r["sample"]), xy, fontsize=6.5, xytext=COLUMN[r["sample"]],
                            textcoords="data", va="center", color="#333333",
                            arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5, shrinkA=0, shrinkB=3))
            else:
                ax.annotate(LABEL.get(r["sample"], r["sample"]), xy, fontsize=6.5,
                            xytext=OFFSET.get(r["sample"], (3, 2)), textcoords="offset points", color="#333333")
    proj = {r["sample"]: r for r in rows if r["source"] == "projected" and r["PC1"] != "nan"}
    for name, key, marker, color in (("Puebla (projected)", "puebla.hifi", "*", "black"),
                                      ("Thille HiFi (projected control)", "thille.hifi", "X", "#009E73")):
        if key in proj:
            r = proj[key]
            ax.scatter(float(r["PC1"]), float(r["PC2"]), s=170 if marker == "*" else 70, marker=marker, c=color,
                       edgecolors="black", linewidths=0.6, label=name, zorder=4)
            if key == "puebla.hifi":
                ax.annotate("Puebla", (float(r["PC1"]), float(r["PC2"])), fontsize=8, fontweight="bold",
                            xytext=(-38, -3), textcoords="offset points")
            print(f"{key}: PC1 {float(r['PC1']):.3f} PC2 {float(r['PC2']):.3f}")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="lower left", frameon=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=600 if ext == "pdf" else 300, bbox_inches="tight")


if __name__ == "__main__":
    main()
