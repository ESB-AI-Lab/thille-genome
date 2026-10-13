#!/usr/bin/env python3
# v4 copy (2026-10-05): Thille v4 genome, annotation, coverage, SyRI and centromere inputs; display name "Unnamed Lowland" (Edwin, 2026-10-05).
# 2026-10-07 (audit A3): Puebla HiFi reads on Thille added as a fourth Thille coverage panel; Puebla reads are not
# aligned to Hass UoQ or the Unnamed Lowland, so those columns keep their read sets.
"""Karyotype, gene, TE and read-coverage tracks for one SV locus in West Indian T2T, Thille and Hass UoQ.

One column per genome (West Indian, Thille, Hass UoQ). Top row: the whole chromosome as an ideogram, with the plotted
window boxed and the SV block (compare/regions/region_blocks.tsv) in red; Thille also shows its centromere call. Below,
the window: genes (Thille Tiberius, Hass UoQ GenBank genes, West Indian GenBank GFF), TEs colored by class, tandem arrays
with period >= 100 bp (TRF), then one read-depth panel per raw read alignment to that genome (bedtools genomecov -d
output reduced to 10 kb bins by compare/regions/cov/bin_genomecov.sh), as depth / chromosome median with the coverage
class boundaries 0.25, 0.75, 1.25 and 1.5 dashed. Windows are those used by plot_region_threeway.py (the Thille block
+/- 1.5 Mb and the span of the SyRI alignments to it). Chromosomes that SyRI used reverse-complemented (Hass and West
Indian chr1, 4, 5, 6, 10, 11) are drawn in that same-strand orientation; tick labels give original GenBank coordinates.

TE classes come from the family name in the RepeatMasker Target, looked up in the Thille RepeatModeler library headers
(consensi.fa.classified), for all three genomes; the Hass and West Indian windows were masked with that same library
(compare/regions/sbatch_repeatmasker_windows.sh).

Usage (from assemblies/thille.hic): plot_region_tracks.py <region_id> <out_prefix>
"""
import argparse
import csv
import glob
import importlib.util
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle  # noqa: E402

spec = importlib.util.spec_from_file_location("tw", "compare/plot_region_threeway_v4.py")
tw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tw)

LIB = "/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/consensi.fa.classified"
THILLE_TE = "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff"
COV = "compare/v4/regions/cov"
CEN = "centromere/v4/thille_centromeres.satellitefirst.bed"
READS = {"thille": [("thille.hifi.thille_v4", "Thille HiFi"), ("westindian.hifi.thille_v4", "Unnamed Lowland HiFi"),
                    ("gwen.clr.thille_v4", "Gwen CLR"), ("puebla.hifi.thille_v4", "Puebla HiFi")],
         "hass": [("thille.hifi.hass_uoq", "Thille HiFi"), ("westindian.hifi.hass_uoq", "Unnamed Lowland HiFi"),
                  ("gwen.clr.hass_uoq", "Gwen CLR")],
         "westindian": [("thille.hifi.wi", "Thille HiFi"), ("gwen.clr.wi", "Gwen CLR")]}
TE_COL = {"LTR/Gypsy": "#b2182b", "LTR/Copia": "#ef8a62", "LTR other": "#fddbc7", "LINE/SINE": "#67a9cf",
          "DNA/RC": "#2166ac", "Unknown": "#bdbdbd", "Simple/low complexity": "#4d4d4d"}
INK, BLOCK, GENE, SAT, IDEO = "#1f1f1f", "#c0392b", "#1f5fa8", "#6d4c41", "#d9d9d9"
CAND = "#d4a017"
FLOWER = "/gpool/edwin/avocados/thille_paper/flowering_beds/igv"
CAND_BEDS = {"thille": ["flowering_genes_thille.bed"], "hass": ["flowering_genes_hass.bed"],
             "westindian": ["flowering_genes_wi.bed", "flowering_genes_wi_fst.bed"]}
SV_COL = {"NOTAL": ("#c0392b", "not aligned"), "INV": ("#e08a00", "inversion"), "TRANS": ("#2e9e44", "translocation"),
          "INVTR": ("#2e9e44", "translocation"), "DUP": ("#8e44ad", "duplication"), "INVDP": ("#8e44ad", "duplication"),
          "DEL": ("#1f1f1f", "deletion"), "INS": ("#17a2b8", "insertion"), "CPG": ("#8c564b", "copy gain/loss"),
          "CPL": ("#8c564b", "copy gain/loss"), "HDR": ("#bdbdbd", "highly diverged")}
SV_MIN = 50000
SHORT = {"westindian": "Unnamed Lowland", "hass": "Hass", "thille": "Thille"}
plt.rcParams.update({"font.size": 7, "text.color": INK, "pdf.fonttype": 42})


def family_classes():
    cls = {}
    for line in open(LIB):
        if line.startswith(">"):
            fam, _, c = line[1:].split()[0].partition("#")
            cls[fam] = c
    return cls


def te_group(c):
    if c.startswith("LTR/Gypsy"):
        return "LTR/Gypsy"
    if c.startswith("LTR/Copia"):
        return "LTR/Copia"
    if c.startswith("LTR"):
        return "LTR other"
    if c.startswith(("LINE", "SINE")):
        return "LINE/SINE"
    if c.startswith(("DNA", "RC")):
        return "DNA/RC"
    if c in ("Simple_repeat", "Low_complexity", "Satellite"):
        return "Simple/low complexity"
    return "Unknown"


def te_rows(gff, cls, chrom_filter=None, offset=0):
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or (chrom_filter and f[0] != chrom_filter):
            continue
        m = re.search(r'Motif:([^"]+)"', f[8])
        fam = m.group(1) if m else ""
        c = "Simple_repeat" if fam.startswith("(") else ("Low_complexity" if fam.endswith("-rich") else cls.get(fam, "Unknown"))
        yield int(f[3]) - 1 + offset, int(f[4]) + offset, te_group(c)


def coverage(name, genome, chrom):
    oid = chrom if genome == "thille" else tw.IDS[genome][chrom]
    path = f"{COV}/{name}.{oid}.10kb.tsv"
    head = open(path).readline()
    med = int(re.search(r"median_depth=(\d+)", head).group(1))
    rows = []
    for line in open(path):
        if line[0] in "#b":
            continue
        s, e, d = line.split("\t")
        s, e = tw.to_same(genome, chrom, int(s), int(e))
        rows.append(((s + e) / 2, float(d) / med if med else 0.0))
    return sorted(rows), med


def candidates(genome, chrom):
    """Flowering candidate genes (top 1% Fst windows) from the flowering_beds BEDs, same-strand coordinates."""
    oid = chrom if genome == "thille" else tw.IDS[genome][chrom]
    seen = {}
    for bed in CAND_BEDS[genome]:
        for line in open(f"{FLOWER}/{bed}"):
            f = line.split("\t")
            if line.startswith("track") or len(f) < 4 or f[0] != oid:
                continue
            name = f[3].split("|")[0]
            seen[(name, f[1])] = (*tw.to_same(genome, chrom, int(f[1]), int(f[2])), name)
    return sorted(seen.values())


def sv_calls(genome, chrom, s0, e0):
    """SyRI annotations >= SV_MIN in this genome's window: lanes of (label, [(start, end, type)])."""
    lanes = []
    if genome == "thille":
        for other in ("westindian", "hass"):
            ev = []
            for line in open(f"compare/v4/syri/{other}/{other}syri.out"):
                f = line.split("\t")
                if f[0] == chrom and f[10] in SV_COL and f[1] != "-":
                    s, e = int(f[1]) - 1, int(f[2])
                    if e - s >= SV_MIN and e > s0 and s < e0:
                        ev.append((s, e, f[10]))
            lanes.append((f"vs {SHORT[other]}", ev))
    else:
        ev = []
        for line in open(f"compare/v4/syri/{genome}/{genome}syri.out"):
            f = line.split("\t")
            if f[5] == chrom and f[10] in SV_COL and f[6] != "-":
                s, e = sorted((int(f[6]) - 1, int(f[7])))
                if e - s >= SV_MIN and e > s0 and s < e0:
                    ev.append((s, e, f[10]))
        lanes.append(("vs Thille", ev))
    return lanes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("region"); ap.add_argument("out")
    ap.add_argument("--no-text", action="store_true", help="omit the title, description and read medians (they go in the figure legend); default off")
    a = ap.parse_args()
    blocks = [r for r in csv.DictReader(open("compare/regions/region_blocks.tsv"), delimiter="\t") if r["region"] == a.region]
    tb = next(b for b in blocks if b["genome"] == "thille")
    chrom = tb["chrom"]
    win = {"thille": (max(0, int(tb["start"]) - 1500000), min(tw.length("thille", chrom), int(tb["end"]) + 1500000))}
    for g in ("westindian", "hass"):
        rows = [r for r in tw.syri_rows(g, chrom, *win["thille"]) if r[4] == "SYNAL"]
        win[g] = (min(r[2] for r in rows), max(r[3] for r in rows))
    block = {b["genome"]: tw.to_same(b["genome"], chrom, int(b["start"]), int(b["end"])) for b in blocks}
    cls = family_classes()
    te = {"thille": [r for r in te_rows(THILLE_TE, cls, chrom) if r[1] > win["thille"][0] and r[0] < win["thille"][1]]}
    trf = {"thille": tw.trf_arrays(f"centromere/v4/trf/{chrom}.dat", 0)}
    for g in ("westindian", "hass"):
        start = win[g][0]
        base = f"win_{a.region}_{g}_{start}"
        gffs = glob.glob(f"compare/regions/repeatmasker/{base}/{base}.fa.out.gff")
        te[g] = list(te_rows(gffs[0], cls, offset=start)) if gffs else []
        trf[g] = tw.trf_arrays(f"compare/regions/trf/{base}.dat", start)
    cen = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in open(CEN)}
    desc = next(r for r in csv.DictReader(open("compare/v4/regions/region_descriptions_relabel2.tsv"), delimiter="\t")
                if r["region"] == a.region)

    order = ["westindian", "thille", "hass"]
    nrows = 5 + max(len(v) for v in READS.values())
    fig = plt.figure(figsize=(13, 2.6 + 1.05 * nrows))
    gs = fig.add_gridspec(nrows, 3, height_ratios=[0.9, 0.7, 0.8, 0.55, 0.3] + [1] * (nrows - 5), hspace=0.5, wspace=0.12)
    for col, g in enumerate(order):
        s0, e0 = win[g]
        L = tw.length(g, chrom)
        oid = chrom if g == "thille" else tw.IDS[g][chrom]
        rc = g != "thille" and chrom in tw.RC
        # ideogram of the whole chromosome
        ax = fig.add_subplot(gs[0, col])
        ax.add_patch(FancyBboxPatch((0, 0.35), L, 0.3, boxstyle="round,pad=0,rounding_size=0.15",
                                    mutation_aspect=L / 8, color=IDEO))
        if g == "thille":
            cs, ce = cen[chrom]
            ax.add_patch(Rectangle((cs, 0.3), ce - cs, 0.4, color="#4a4a4a"))
        if g in block:
            ax.add_patch(Rectangle((block[g][0], 0.3), block[g][1] - block[g][0], 0.4, color=BLOCK))
        ax.add_patch(Rectangle((s0, 0.22), e0 - s0, 0.56, fill=False, ec=INK, lw=0.9))
        ax.set_xlim(-L * 0.01, L * 1.01); ax.set_ylim(0, 1); ax.axis("off")
        ax.set_title(f"{tw.NAME[g]}  {oid}{' (rc)' if rc else ''}" if a.no_text else
                     f"{tw.NAME[g]}  {oid}{' (reverse-complemented)' if rc else ''}  {L / 1e6:.1f} Mb",
                     fontsize=8, fontweight="bold")
        # SyRI structural annotations >= SV_MIN, labelled by type
        sx = fig.add_subplot(gs[1, col])
        lanes = sv_calls(g, chrom, s0, e0)
        for k, (lab, ev) in enumerate(lanes):
            y = 1 - (k + 1) / len(lanes)
            for s, e, t in ev:
                c, name = SV_COL[t]
                sx.add_patch(Rectangle((max(s, s0), y + 0.08), min(e, e0) - max(s, s0), 0.84 / len(lanes), color=c, lw=0))
                if min(e, e0) - max(s, s0) > (e0 - s0) * 0.06 and t != "HDR":
                    sx.text((max(s, s0) + min(e, e0)) / 2, y + 0.5 / len(lanes), name, ha="center", va="center",
                            fontsize=5.5, color="white" if t in ("NOTAL", "DEL", "DUP", "INVDP") else INK)
            sx.text(s0, y + 0.92 / len(lanes) + 0.02, lab, ha="left", va="bottom", fontsize=6, color="#5f5f5f")
        sx.set_xlim(s0, e0); sx.set_ylim(0, 1); sx.axis("off")
        if col == 0:
            sx.text(-0.16, 0.5, "SyRI calls\n(>= 50 kb)", transform=sx.transAxes, ha="right", va="center", fontsize=7)
        # genes (flowering candidates highlighted and named), TEs, satellites in the window
        tracks = [("genes", [(s, e, GENE) for s, e in tw.genes(g, chrom, s0, e0)]),
                  ("TEs", [(s, e, TE_COL[c]) for s, e, c in te[g] if e > s0 and s < e0]),
                  ("tandem\narrays", [(s, e, SAT) for s, e in trf[g] if e > s0 and s < e0])]
        for i, (label, feats) in enumerate(tracks):
            tx = fig.add_subplot(gs[2 + i, col])
            for s, e, c in feats:
                tx.add_patch(Rectangle((max(s, s0), 0.1), max(min(e, e0) - max(s, s0), (e0 - s0) * 0.0008), 0.8,
                                       color=c, lw=0))
            if g in block:
                tx.axvspan(*block[g], color=BLOCK, alpha=0.12, lw=0)
            if label == "genes":
                for cs_, ce_, name in candidates(g, chrom):
                    if ce_ > s0 and cs_ < e0:
                        mid = (cs_ + ce_) / 2
                        tx.add_patch(Rectangle((mid - (e0 - s0) * 0.003, 0.0), (e0 - s0) * 0.006, 1.0, color=CAND, lw=0,
                                               zorder=5))
                        tx.text(mid, 1.05, name, ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="#8a6500",
                                clip_on=False)
            tx.set_xlim(s0, e0); tx.set_ylim(0, 1); tx.set_yticks([]); tx.set_xticks([])
            for sp in tx.spines.values():
                sp.set_visible(False)
            if col == 0:
                tx.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=7)
        # coverage panels
        for j, (name, label) in enumerate(READS[g]):
            cx = fig.add_subplot(gs[5 + j, col])
            pts, med = coverage(name, g, chrom)
            pts = [(x, y) for x, y in pts if s0 <= x <= e0]
            cx.fill_between([x for x, _ in pts], [min(y, 2.0) for _, y in pts], step="mid", color="#7f8c8d", lw=0)
            for yv in (0.25, 0.75, 1.25, 1.5):
                cx.axhline(yv, color="#9e9e9e", lw=0.5, ls="--")
            cx.axhline(1, color=INK, lw=0.7)
            if g in block:
                cx.axvspan(*block[g], color=BLOCK, alpha=0.12, lw=0)
            cx.set_xlim(s0, e0); cx.set_ylim(0, 2)
            cx.set_title(label if a.no_text else f"{label} reads (median {med}x)", fontsize=7, loc="left", pad=2)
            if col == 0:
                cx.set_ylabel("depth / median", fontsize=6)
            last = j == len(READS[g]) - 1
            ticks = [s0 + (e0 - s0) * k / 4 for k in range(5)]
            cx.set_xticks(ticks)
            cx.set_xticklabels([f"{tw.to_orig(g, chrom, t) / 1e6:.2f}" for t in ticks] if last else [])
            if last:
                cx.set_xlabel(f"{oid} position (Mb)", fontsize=7)
    handles = [Patch(color=BLOCK, label="SV block"), Patch(color=GENE, label="gene"),
               Patch(color=CAND, label="flowering candidate (top 1% Fst windows)"),
               Patch(color=SAT, label="tandem array, period >= 100 bp")]
    seen = {}
    for c, name in SV_COL.values():
        seen.setdefault(name, c)
    handles += [Patch(color=c, label=f"SyRI: {n}") for n, c in seen.items()]
    handles += [Patch(color=c, label=f"TE: {k}") for k, c in TE_COL.items()]
    handles += [Patch(color="#4a4a4a", label="Thille centromere call")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.01))
    if not a.no_text:
        fig.subplots_adjust(top=0.93)
        fig.suptitle(f"{a.region}: {chrom} {desc['sv_type']}", fontsize=10, y=0.985, fontweight="bold")
        fig.text(0.5, 0.968, desc["description"], ha="center", va="top", fontsize=7.5, wrap=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300, bbox_inches="tight")
    print(f"wrote {a.out}.pdf/.png; windows " + ", ".join(f"{g} {win[g][0]}-{win[g][1]}" for g in order))


if __name__ == "__main__":
    main()
