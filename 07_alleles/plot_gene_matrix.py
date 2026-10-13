#!/usr/bin/env python3
"""Figure: variant x haplotype matrix ("3b overview") for one Fst flowering gene, same design as
scripts/plot_sdmyb_matrix.py, from summaries/fst_alleles/<label>/variants_all.tsv (scripts/gene_allele_summary.py).

a  the locus to scale (Thille window, gene drawn 5' -> 3', CDS exons as boxes), regions linked to the column groups
b  matrix: columns = sites carried by at least one haplotype, grouped 5' flank | E1 | I1 | E2 | ... | 3' flank
   (exons and introns numbered in transcript order; "gene" when there is no Thille CDS model); rows = Gwen h1/h2,
   Thille h1/h2, Hass pri/alt (UoQ pseudohaplotypes), Puebla hap1/hap2, WI h1/h2; a cell is filled when that haplotype carries the variant (reads:
   >= --carry of its reads covering the site; assemblies: present), colored by variant type; the allele letter
   (GENE^a = Thille reference CDS, summaries/fst_alleles/allele_summary.tsv) is printed at the right of each row;
   under each column, the number of haplotypes carrying it
c  carried sites per region and haplotype (more than 9 regions: 5' flank, exons, introns, 3' flank)
Sites are merged across haplotypes when position, type, length and (SNP) alternative base match.

Usage (from the repo root): python3 scripts/plot_gene_matrix.py <label> [--carry 0.5] [--no-text]
Writes Results_fst_alleles/<label>_variant_matrix.{pdf,png}; with --no-text (Edwin, 2026-10-09: no title, region counts or
footnote inside the figure) <label>_variant_matrix_notext.{pdf,png} and <label>_variant_matrix_notext_notes.tsv (the
title, the sites per region and the footnote, for the legend md)
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import ConnectionPatch, Patch, Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fst_common import OUT, load_models  # noqa: E402

ROWS = [("Gwen", 1), ("Gwen", 2), ("Thille", 1), ("Thille", 2), ("Hass", 1), ("Hass", 2), ("Puebla", 1), ("Puebla", 2),
        ("WI", 1), ("WI", 2)]  # display order (Edwin, 2026-10-01); allele letters keep the order in gene_allele_summary.py
ROW_GROUP = {"Thille": "Thille  (B · HiFi)", "WI": "Unnamed Lowland  (? · HiFi)", "Gwen": "Gwen  (A · CLR)",
             "Puebla": "Puebla  (A · assembly)", "Hass": "Hass  (A · assembly)"}
HLAB = {"Puebla": "hap"}
TABLE_NAME = {"WI": "Unnamed Lowland"}  # display name (Edwin, 2026-10-05); "WI" stays the internal key
COL = {"nonsyn": "#C8501E", "syn": "#4F9FD6", "noncoding": "#8E959F", "DEL": "#1D3557", "INS": "#E0A100"}
GREY = "#B9BEC5"


def rowname(c, h):
    return {1: "pri", 2: "alt"}[h] if c == "Hass" else f"{HLAB.get(c, 'h')}{h}"


def region_label(r):
    return {"5prime_flank": "5′ flank", "3prime_flank": "3′ flank", "gene": "gene"}.get(r, r)


def vtype(r):
    if r["type"] in ("DEL", "INS"):
        return r["type"]
    if isinstance(r["aa_change"], str) and r["aa_change"]:
        return "syn" if r["aa_change"][0] == r["aa_change"][-1] else "nonsyn"
    return "noncoding"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("label")
    p.add_argument("--carry", type=float, default=0.5)
    p.add_argument("--outdir", default="Results_fst_alleles")
    p.add_argument("--no-text", action="store_true", help="no title, region counts or footnote in the figure; they go to "
                   "<out>_notes.tsv (default off; output without it is unchanged)")
    a = p.parse_args()
    g = load_models()[a.label]
    REGIONS = g.regions()
    v = pd.read_csv(OUT / a.label / "variants_all.tsv", sep="\t")
    v["vtype"] = v.apply(vtype, axis=1)
    letters, afp = {}, {}
    for r in csv.DictReader(open(OUT / "allele_summary.tsv"), delimiter="\t"):
        if r["label"] == a.label:
            cv, h = r["haplotype"].split(" h")
            letters[(cv, int(h))] = r["allele"].split("^")[1] if "^" in r["allele"] else ""
            afp[(cv, int(h))] = r.get("allele_freq_pct", "")
    # covered intervals per row: reads (>= 3 reads of that haplotype) and assemblies (aligned blocks)
    covered = {}
    for rs in ("Thille", "WI", "Gwen"):
        f = OUT / a.label / f"{rs}.covered.tsv"
        if f.exists():
            for r in csv.DictReader(open(f), delimiter="\t"):
                covered.setdefault((rs, int(r["hap"])), []).append((int(r["start"]), int(r["end"])))
    for name, row in (("Puebla_hap1", ("Puebla", 1)), ("Puebla_hap2", ("Puebla", 2)), ("Hass_pri", ("Hass", 1)),
                      ("Hass_alt", ("Hass", 2))):
        f = OUT / a.label / f"{name}.paf"
        covered[row] = []
        if f.exists():
            for line in open(f):
                x = line.split("\t")
                if "tp:A:P" in x[12:]:
                    covered[row].append((g.win_start + int(x[7]), g.win_start + int(x[8]) - 1))

    def is_covered(row, pos):
        return any(s <= pos <= e for s, e in covered.get(row, []))

    carried = {}
    for r in v.itertuples(index=False):
        key = (r.pos, r.type, r.length if r.type != "SNP" else 0, r.alt if r.type == "SNP" else "")
        for h in (1, 2):
            w, c = getattr(r, f"hap{h}_with"), getattr(r, f"hap{h}_covering")
            flag = getattr(r, f"hap{h}_carried", float("nan"))
            hit = (flag == 1) if flag == flag else (c and w / c >= a.carry)  # read sets carry an explicit call
            if (r.read_set, h) in ROWS and hit:
                carried.setdefault(key, {"region": r.region, "vtype": r.vtype, "aa": r.aa_change, "rows": set()})
                carried[key]["rows"].add((r.read_set, h))
    rev = g.strand == "-"
    sites = sorted(carried.items(), key=lambda kv: (REGIONS.index(kv[1]["region"]), -kv[0][0] if rev else kv[0][0]))
    n_cds = sum(1 for _, d in sites if d["region"][0] == "E")
    n_nc = len(sites) - n_cds
    w_cds, w_nc = 6.0, min(1.0, 300.0 / max(n_nc, 1))  # coding columns wide enough for their labels

    xs, widths, groups, x = [], [], [], 0.0
    for reg in REGIONS:
        cols = [i for i, (_, d) in enumerate(sites) if d["region"] == reg]
        if not cols:
            continue
        start = x
        w = w_cds if reg[0] == "E" else w_nc
        for _ in cols:
            xs.append(x); widths.append(w); x += w
        groups.append((reg, start, x, len(cols)))
        x += 4.0
    total = max(x, 1.0)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.linewidth": 0.5})
    fig = plt.figure(figsize=(7.2, 7.2))
    gs = fig.add_gridspec(4, 1, height_ratios=[0.8, 4.2, 0.5, 1.8], hspace=0.45)
    ax_loc, ax_m, ax_bar, ax_tab = (fig.add_subplot(gs[i]) for i in range(4))

    # a: locus to scale, 5' -> 3'
    span = (g.win_end, g.win_start) if rev else (g.win_start, g.win_end)
    ax_loc.set_xlim(span); ax_loc.set_ylim(-1, 1.6); ax_loc.axis("off")
    ax_loc.plot(span, [0, 0], color="#4a4a4a", lw=0.6)
    if g.exons:
        for s, e in g.exons:
            ax_loc.add_patch(Rectangle((s, -0.35), e - s + 1, 0.7, color="#2b2b2b", lw=0))
    else:
        ax_loc.add_patch(Rectangle((g.gene_start, -0.2), g.gene_end - g.gene_start + 1, 0.4, color="#8a8a8a", lw=0))
    g5, g3 = (g.gene_end, g.gene_start) if rev else (g.gene_start, g.gene_end)
    ax_loc.annotate("", xy=(g3, 0.9), xytext=(g5, 0.9), arrowprops=dict(arrowstyle="->", lw=0.6, color="#4a4a4a"))
    model = g.thille_model if g.exons else "no Thille CDS model"
    title = f"{g.gene} ({model}) on Thille {g.chrom}, 5′ → 3′"
    if not a.no_text:
        ax_loc.text((g5 + g3) / 2, 1.05, title, ha="center", va="bottom", fontsize=6.5)
    span_bp = g.win_end - g.win_start
    step = next(s * 10 ** k for k in range(2, 8) for s in (1, 2, 5) if span_bp / (s * 10 ** k) <= 5)
    for t in range(int(np.ceil(g.win_start / step) * step), g.win_end + 1, step):
        ax_loc.text(t, -0.95, f"{t / 1e6:.3f} Mb", ha="center", va="bottom", fontsize=5.5, color="#555555")
    ax_loc.text(span[0], 1.55, "a", fontsize=9, fontweight="bold", va="top")

    # b: matrix
    ax_m.set_xlim(-0.5, total); ax_m.set_ylim(len(ROWS) - 0.5 + 0.6, -0.9)
    for j, (key, d) in enumerate(sites):
        for i, row in enumerate(ROWS):
            if row in d["rows"]:
                ax_m.add_patch(Rectangle((xs[j] + 0.04 * widths[j], i - 0.38), 0.92 * widths[j], 0.76, color=COL[d["vtype"]], lw=0))
            elif is_covered(row, key[0]):
                ax_m.add_patch(Rectangle((xs[j] + 0.04 * widths[j], i - 0.38), 0.92 * widths[j], 0.76, color="#F1F2F4", lw=0))
            else:  # no reads / no alignment for this haplotype here: missing, not reference
                ax_m.add_patch(Rectangle((xs[j] + 0.04 * widths[j], i - 0.38), 0.92 * widths[j], 0.76, facecolor="white",
                                         edgecolor="#B9BEC5", hatch="//////", lw=0))
    for reg, s, e, k in groups:
        bold = reg[0] == "E"
        ax_m.plot([s, e], [-0.65, -0.65], color="#2b2b2b" if bold else "#8a8a8a", lw=1.0 if bold else 0.5)
        small = (e - s) < 0.12 * total
        ax_m.text((s + e) / 2, -0.75, region_label(reg) if small or a.no_text else f"{region_label(reg)} (n={k})", ha="center",
                  va="bottom", fontsize=6 if len(groups) < 12 else 5, fontweight="bold" if bold else "normal")
        lo, hi = g.region_span(reg)
        g0, g1 = (hi, lo) if rev else (lo, hi)
        for gx, mx in ((g0, s), (g1, e)):
            fig.add_artist(ConnectionPatch(xyA=(gx, -0.45), coordsA=ax_loc.transData, xyB=(mx, -1.15),
                                           coordsB=ax_m.transData, color="#c4c8cd", lw=0.4))
    label_all = n_cds <= 60
    for j, (key, d) in enumerate(sites):
        if d["region"][0] == "E" and (label_all or d["vtype"] == "nonsyn"):
            lab = d["aa"] if isinstance(d["aa"], str) and d["aa"] else f"{key[1]}{key[2]}"
            ax_m.text(xs[j] + widths[j] / 2, len(ROWS) - 0.35, lab, rotation=90, ha="center", va="top", fontsize=5)
    ax_m.set_yticks(range(len(ROWS)), [rowname(c, h) for c, h in ROWS], fontsize=6)
    ax_m.tick_params(length=0); ax_m.set_xticks([])
    for side in ax_m.spines.values():
        side.set_visible(False)
    for grp, label in ROW_GROUP.items():
        idx = [i for i, (c, _) in enumerate(ROWS) if c == grp]
        ax_m.text(-0.035 * total, np.mean(idx), label, ha="right", va="center", fontsize=6.5)
        ax_m.plot([-0.4, total - 4.0], [max(idx) + 0.5] * 2, color="#d0d0d0", lw=0.4)
    xr = total - 3.0
    for i, row in enumerate(ROWS):
        cnt = sum(row in d["rows"] for _, d in sites)
        ax_m.text(xr, i, str(cnt), ha="left", va="center", fontsize=6, color="#444444")
        if letters.get(row):
            ax_m.text(xr + 0.06 * total, i, f"{g.gene}$^{{{letters[row]}}}$", ha="left", va="center", fontsize=6.5,
                      style="italic")
            if afp.get(row):
                ax_m.text(xr + 0.06 * total + 0.04 * total + 0.012 * total * len(g.gene), i, f"{afp[row]}%" if afp[row] != "n.d." else "n.d.",
                          ha="left", va="center", fontsize=6, color="#444444")
    ax_m.text(xr, -0.2, "sites", ha="left", va="bottom", fontsize=5.5, color="#444444")
    if letters:
        ax_m.text(xr + 0.06 * total, -0.2, "allele  (AF)", ha="left", va="bottom", fontsize=5.5, color="#444444")
    ax_m.text(-0.5, -1.25, "b", fontsize=9, fontweight="bold", va="bottom")

    # carrier bar
    ax_bar.set_xlim(-0.5, total); ax_bar.set_ylim(0, len(ROWS))
    for j, (key, d) in enumerate(sites):
        ax_bar.add_patch(Rectangle((xs[j] + 0.04 * widths[j], 0), 0.92 * widths[j], len(d["rows"]), color=GREY, lw=0))
    ax_bar.set_yticks([0, 5, 10]); ax_bar.set_ylabel("haplotypes\ncarrying", fontsize=5.5)
    ax_bar.set_xticks([])
    for s in ("top", "right", "bottom"):
        ax_bar.spines[s].set_visible(False)

    # c: carried sites per region and haplotype
    ax_tab.axis("off")
    regs = [gr[0] for gr in groups]
    if len(regs) <= 9:
        cols = [(region_label(r), {r}) for r in regs]
    else:
        cols = [("5′ flank", {"5prime_flank"}), ("exons", {r for r in regs if r[0] == "E"}),
                ("introns", {r for r in regs if r[0] == "I"}), ("3′ flank", {"3prime_flank"})]
    cols.append(("total", set(regs)))
    cells = [[str(sum(1 for _, d in sites if d["region"] in rs and row in d["rows"])) for _, rs in cols] for row in ROWS]
    tbl = ax_tab.table(cellText=cells, rowLabels=[f"{TABLE_NAME.get(c, c)} {rowname(c, h)}" for c, h in ROWS],
                       colLabels=[c for c, _ in cols], loc="center", cellLoc="center", edges="horizontal")
    tbl.auto_set_font_size(False); tbl.set_fontsize(6); tbl.scale(1, 0.85)
    for (_, _), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3); cell.set_edgecolor("#9aa0a6")
    ax_tab.text(-0.07, 1.15, "c", fontsize=9, fontweight="bold", va="bottom", transform=ax_tab.transAxes)

    handles = [Patch(color=COL[k], label=l) for k, l in (("nonsyn", "SNP, amino-acid change"), ("syn", "SNP, synonymous"),
                                                          ("noncoding", "SNP, noncoding"), ("DEL", "deletion ≥ 2 bp"),
                                                          ("INS", "insertion ≥ 2 bp"))]
    handles += [Patch(facecolor="#F1F2F4", label="reference allele"),
                Patch(facecolor="white", edgecolor="#B9BEC5", hatch="//////", label="no data")]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False, fontsize=6, bbox_to_anchor=(0.5, -0.01))
    homo = sorted({r.read_set for r in v.itertuples() if str(getattr(r, "unphased", "")) in ("1", "1.0")})
    note = (f"Carried = the haplotype's allele: in ≥ {int(a.carry * 100)}% of its reads, else its most frequent allele; "
            "unphased read sets: most frequent allele on h1, second (≥ 25% of reads) on h2 (assemblies: present in the aligned sequence). Variants are relative to the Thille assembly (one Thille "
            "haplotype).\nAllele letters: identical CDS = same letter; a = Thille reference CDS; ? = CDS not fully "
            "covered by that haplotype's reads. AF = share of the variety's reads carrying that exon allele, pooled over "
            "CDS sites where its two haplotypes differ\n(top two alleles per site; intronic and flanking variation not "
            "used); 100% = both haplotypes share the exon allele; unphased read sets with < 25% second allele: majority allele on both rows, AF = its real frequency; assemblies: 50% or 100% of the two assembled haplotypes.")
    if homo:
        note += (f" {', '.join(homo)}: not phased (no two reads share enough heterozygous sites); h1/h2 = first and "
                 "second most frequent allele.")
    if not g.exons:
        note += " No Thille CDS model at this locus: regions are 5′ flank, placed gene span, 3′ flank; no allele letters."
    Path(a.outdir).mkdir(parents=True, exist_ok=True)
    out = f"{a.outdir}/{a.label}_variant_matrix"
    if a.no_text:
        out += "_notext"
        with open(f"{out}_notes.tsv", "w") as fh:
            fh.write(f"title\t{title}\n")
            fh.write("regions\t" + "; ".join(f"{region_label(reg)} {k}" for reg, _, _, k in groups) + "\n")
            fh.write(f"note\t{' '.join(note.split())}\n")
    else:
        fig.text(0.01, -0.045, note, fontsize=5.5, color="#555555", linespacing=1.4)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"{a.label}: {len(sites)} sites ({n_cds} coding); " +
          ", ".join(f"{c}{h} {sum((c, h) in d['rows'] for _, d in sites)}" for c, h in ROWS))


if __name__ == "__main__":
    main()
