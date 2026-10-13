#!/usr/bin/env python3
"""Summary of the alleles carried by all three A-type cultivars and absent from Thille (observation, not a test).

Genes (Edwin, 2026-10-02): SDMYB, AGAMOUS-like, CSU2, NPR5-BOP2, SGRL, LCYB. For each gene, the A-shared allele is the
CDS allele (summaries/fst_alleles/allele_summary.tsv) present in Gwen, Hass (MH2022) and Puebla and in neither Thille
haplotype. The figure shows the allele letter and AF of every haplotype, the A-shared allele highlighted, its
amino-acid changes, and the gene's Fst status (contig and scaffold window ranks from manuscript Table 1; the top 1% is
ranks 1 to 517 of 51,701 contig windows and 1 to 508 of 50,775 scaffold windows). Fst values not available are shown
as "n/a" and listed in summaries/fst_alleles/a_shared_alleles.tsv as [needs updating].

Usage (from the repo root): python3 scripts/plot_a_shared_alleles.py [--no-text]
Writes summaries/fst_alleles/a_shared_alleles.tsv and Results_fst_alleles/a_shared_alleles.{pdf,png}. With --no-text
(Edwin, 2026-10-09: no explanatory text inside figures) only Results_fst_alleles/a_shared_alleles_notext.{pdf,png}: the
footnote is left out (it belongs in the legend md) and the column header is shortened; the TSV is not rewritten.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

GENES = ["SDMYB", "AGAMOUS-like", "CSU2", "NPR5-BOP2", "SGRL", "LCYB"]
COLS = [("Gwen", 1, "h1"), ("Gwen", 2, "h2"), ("Thille", 1, "h1"), ("Thille", 2, "h2"), ("Hass", 1, "h1"),
        ("Hass", 2, "h2"), ("Puebla", 1, "h1"), ("Puebla", 2, "h2"), ("WI", 1, "h1"), ("WI", 2, "h2")]
# Display order (Edwin, 2026-10-07): haplotype phase is not shared between loci, so within each A-type cultivar the
# haplotype carrying the A-shared allele is drawn second for every gene (Hass pri/alt and Puebla hap1/hap2 per gene are
# in the carriers column of a_shared_alleles.tsv); Thille and Unnamed Lowland keep the h1/h2 order of the calls.
A_TYPE = ("Gwen", "Hass", "Puebla")
GROUP = {"Gwen": "Gwen\ntype A", "Thille": "Thille\ntype B", "Hass": "Hass\ntype A", "Puebla": "Puebla\ntype A",
         "WI": "Unnamed\nLowland\ntype ?"}
# Fst status: manuscript v14 Table 1 (contig / scaffold rank, mean Fst); LCYB contig and scaffold ranks from
# summaries/fst_alleles/pigment_fst.tsv (no mean Fst given there)
FST = {"SDMYB": ("138 (0.283)", "157 (0.293)", "both"), "AGAMOUS-like": ("274 (0.222)", "765 (0.173)", "contig"),
       "CSU2": ("719 (0.168)", "326 (0.224)", "scaffold"), "NPR5-BOP2": ("1,925 (0.123)", "124 (0.316)", "scaffold"),
       "LCYB": ("5,297", "5,869", "none"),
       "SGRL": ("11,968 (0.056)", "22,445 (0.033)", "none")}
# SGRL: tblastn of the Thille protein g32057 against the Gwen contigs and scaffolds (the Arabidopsis SGRL protein has no
# hit at 1e-20), by the thille-hic-bf session 2026-10-02; best rank of the two Gwen copies (placed copy: contig 33,842)
ACCENT, FAINT = "#C8501E", "#EEF0F3"
# frameshift alleles: the residue-by-residue list after the shift is replaced by its description
AA_TEXT = {"SDMYB": "N7S, D65G, S139G, D144N, T158I, T163M, A167T, S170N; frameshift from codon 183",
           "AGAMOUS-like": "1 bp deletion at CDS position 5; frameshift, stop at codon 21"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-text", action="store_true", help="no footnote, short column header; writes *_notext only "
                    "(default off; output without it is unchanged)")
    no_text = ap.parse_args().no_text
    t = defaultdict(dict)
    for r in csv.DictReader(open("summaries/fst_alleles/allele_summary.tsv"), delimiter="\t"):
        cv, h = r["haplotype"].split(" h")
        t[r["label"]][(cv, int(h))] = r
    rows = []
    for g in GENES:
        let = {k: (v["allele"].split("^")[1] if "^" in v["allele"] else "") for k, v in t[g].items()}
        per = lambda cv: {let[(cv, 1)], let[(cv, 2)]}
        shared = sorted((per("Gwen") & per("Hass") & per("Puebla")) - per("Thille") - {"?", "n/a", ""})
        assert len(shared) == 1, (g, shared)
        s = shared[0]
        carrier = next(k for k in t[g] if let[k] == s)
        rows.append({"gene": g, "a_shared_allele": f"{g}^{s}", "aa_changes": AA_TEXT.get(g) or (t[g][carrier]["aa_diffs_vs_a"].replace(",", ", ") or "synonymous only"),
                     "in_lowland_wi": "no" if s not in per("WI") else "yes",
                     "carriers": ";".join(f"{cv} h{h}" for (cv, h) in sorted(t[g]) if let[(cv, h)] == s),
                     "fst_contig_rank": FST[g][0], "fst_scaffold_rank": FST[g][1], "top1pct": FST[g][2],
                     "note": "[needs updating] Fst" if "n/a" in FST[g] else ""})
        rows[-1]["_let"], rows[-1]["_af"], rows[-1]["_s"] = let, {k: v["allele_freq_pct"] for k, v in t[g].items()}, s
    cols = [c for c in rows[0] if not c.startswith("_")]
    if not no_text:
        with open("summaries/fst_alleles/a_shared_alleles.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7})
    fig, ax = plt.subplots(figsize=(7.2, 3.3))
    ax.set_xlim(-3.2, len(COLS) + 9.5); ax.set_ylim(len(GENES) + 0.2, -1.9); ax.axis("off")
    prev = None
    for j, (cv, h, lab) in enumerate(COLS):
        ax.text(j + 0.5, -0.15, lab, ha="center", va="bottom", fontsize=6)
        if cv != prev:
            ax.text(j + 1.0, -0.6, GROUP[cv], ha="center", va="bottom", fontsize=6, fontweight="bold", linespacing=1.1)
            if prev is not None:
                ax.plot([j, j], [-0.9, len(GENES)], color="#9aa0a6", lw=0.5)
        prev = cv
    for i, r in enumerate(rows):
        ax.text(-0.15, i + 0.5, r["gene"], ha="right", va="center", fontsize=7, style="italic")
        flip = {cv for cv in A_TYPE if r["_let"].get((cv, 1)) == r["_s"] and r["_let"].get((cv, 2)) != r["_s"]}
        for j, (cv, h0, lab) in enumerate(COLS):
            h = 3 - h0 if cv in flip else h0
            l = r["_let"].get((cv, h), "")
            hit = l == r["_s"]
            ax.add_patch(Rectangle((j + 0.05, i + 0.08), 0.9, 0.84, color=ACCENT if hit else FAINT, lw=0))
            ax.text(j + 0.5, i + 0.42, l, ha="center", va="center", fontsize=7.5, color="white" if hit else "#333333",
                    fontweight="bold" if hit else "normal", style="italic")
            af = r["_af"].get((cv, h), "")
            if af:
                ax.text(j + 0.5, i + 0.78, af if af == "n.d." else f"{af}%", ha="center", va="center", fontsize=4.8,
                        color="white" if hit else "#555555")
        x0 = len(COLS) + 0.3
        g, s = r["gene"], r["_s"]
        ax.text(x0, i + 0.32, f"{g}$^{{{s}}}$ vs {g}$^{{a}}$: " + r["aa_changes"], ha="left",
                va="center", fontsize=5.5)
        ax.text(x0, i + 0.72, f"Fst rank contig {r['fst_contig_rank']}, scaffold {r['fst_scaffold_rank']}", ha="left",
                va="center", fontsize=5.3, color="#444444")
    ax.text(len(COLS) + 0.3, -0.15, "Amino-acid changes; Fst rank" if no_text else "A-shared allele (filled) vs Thille reference allele a:\namino-acid changes; type A vs B Fst rank (mean Fst)", ha="left", va="bottom", fontsize=6,
            fontweight="bold")
    if not no_text:
        fig.text(0.02, 0.01, "Filled cells: the allele carried by all three A-type cultivars and by neither Thille haplotype "
                 "(observation), drawn on the second haplotype of each A-type cultivar (haplotype order is not phased between genes). Letters: CDS alleles, a = Thille reference CDS; numbers: AF. Amino-acid changes: reference residue, position, A-shared residue (N7S: Asn7 to Ser).\nUnnamed Lowland carries "
                 "none of the six CDS alleles (at the protein level it encodes the SGRL$^{d}$ protein on both haplotypes). Fst ranks: Gwen 20 kb windows (top 1% = ranks 1 to 517 contig, 1 to 508 scaffold).",
                 fontsize=5.3, color="#555555", linespacing=1.4)
    Path("Results_fst_alleles").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"Results_fst_alleles/a_shared_alleles{'_notext' if no_text else ''}.{ext}", dpi=600, bbox_inches="tight")
    for r in rows:
        print({c: r[c] for c in cols})


if __name__ == "__main__":
    main()
