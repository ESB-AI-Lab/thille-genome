#!/usr/bin/env python3
"""Legend md files for the figures drawn by this repo (Edwin, 2026-10-09: no figure text inside the PDFs/PNGs; titles,
subtitles, notes and explanatory numbers go to the legend). Figures 1, 2, 6, 7 and Supplementary Figures S5 to S8,
S11, S13 and S14; the thille-hic-bf session writes 3, 4, 5, S2, S9, S10, S12 and S15 in the same folder.

Format (agreed with thille-hic-bf): manuscript/figure_legends/Figure_<id>.md
  # Figure <id>                      (or "# Supplementary Figure S<n>")
  **<title>.**
  (A-C) <shared text>                shared details first, then one line per panel ID (Edwin, 2026-10-10: a title,
  (A) <panel A text>                 then a sentence or two per panel); a figure without panels has one body line
Plain text; *italic* only for gene and species names.

The text moved out of the figures comes from the *_notes.tsv sidecars written with --no-text by plot_gene_stacks.py
(chromosome medians) and plot_gene_matrix.py (title, sites per region, footnote), and from the composition TSV of
Figure 2; the rest is the v37 legend text.
Usage (from the repo root, after the --no-text figures exist): python3 scripts/write_figure_legends.py
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_paths_v4 import B6, MOVED_COVERAGE, S11, S5, V3_TO_V4  # noqa: E402

OUT = Path("manuscript/figure_legends")
NT = Path("Results_v4_notext")
FA = Path("v4work/Results_fst_alleles")
COMP = Path("/epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/thille_genome_composition.tsv")
TRACK_NAMES = {"Gwen": "Gwen CLR", "Thille": "Thille HiFi", "Unnamed Lowland": "Unnamed Lowland HiFi", "Puebla": "Puebla HiFi"}


def v(g):
    return V3_TO_V4.get(g, g)


def matrix(g):
    """title, regions and footnote of a --no-text allele matrix."""
    d = dict(line.rstrip("\n").split("\t", 1) for line in open(FA / f"{v(g)}_variant_matrix_notext_notes.tsv"))
    title = re.sub(r" on Thille chr(\w+), 5′ → 3′", r", chromosome \1", d["title"])
    return title, d["regions"].replace(";", ","), d["note"]


def medians(path):
    """'Gwen CLR 76, Thille HiFi 80, ...' and the figure title from a plot_gene_stacks.py notes TSV."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    parts = []
    for r in rows:
        reads = r["track"].split(" → ")[0]
        parts.append(f"{TRACK_NAMES.get(reads, reads)} {r['chrom_median']}")
    return rows[0]["title"], rows[0]["chrom"], ", ".join(parts)


def gene_id(g):
    """'SDMYB (g34985, chromosome 10)' from the matrix title."""
    t = matrix(g)[0]                       # "SDMYB (g34985), chromosome 10"
    m = re.match(r"(.+?) \((.+?)\), (chromosome \w+)", t)
    return f"{m.group(1)} ({m.group(2)}, {m.group(3)})" if m else t


def sites(g):
    return sum(int(x.split()[-1]) for x in matrix(g)[1].split(","))


def cov(path):
    """'SDMYB (g34985, chromosome 10; medians 76x, 80x, 81x, 30x)': the figure title and panel medians that were drawn
    in the figure, in track order."""
    title, chrom, _ = medians(path)
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    gene = re.split(r" on (?:Thille|Hass|Unnamed Lowland)", title)[0]
    gid = re.match(r"(.+?) \((.+?)\)$", gene)
    meds = ", ".join(f"{r['chrom_median']}x" for r in rows)
    chrom = chrom.replace("Chromosome", "chromosome")
    return f"{gid.group(1)} ({gid.group(2)}, {chrom}; medians {meds})" if gid else f"{gene} ({chrom}; medians {meds})"


COV_SHARED = ("Read depth per position in a 100 kb window centered on the gene (shaded), colored by coverage class "
              "relative to the chromosome median (solid line); below, the class strip, gene models, and TEs. Genes on "
              "the minus strand are drawn right to left.")
MATRIX_SHARED = ("(a) the Thille locus, gene drawn 5′ to 3′; (b) variant sites (columns) carried by each haplotype "
                 "(rows), with the allele letter and allele frequency (AF) of each haplotype; (c) carried sites per "
                 "region and haplotype. Allele a is the Thille reference CDS.")


def write(name, lines):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {OUT / name}.md")


def main():
    comp = {(r["class"], r["subclass"]): float(r["percent_of_assembly_exclN"])
            for r in csv.DictReader(open(COMP), delimiter="\t")}
    rep = 100 - comp[("Non-repetitive", "")]
    ft, pt, wt = NT / "fst_coverage/thille_reference", NT / "pigment_coverage/thille_reference", NT / "fst_coverage/wi_reference"
    cpath = {g: (pt if g in ("SGRL", "LCYB") else ft) / f"{v(g)}_100kb_stacked_notes.tsv" for g in B6}
    write("Figure_1", [
        "# Figure 1",
        "**Cumulative length of the Thille genome assembly.**",
        "Cumulative length of the 12 chromosomes (885.4 Mb) and of the 2,021 hifiasm contigs (989.3 Mb), largest "
        "first. The dashed lines mark the expected genome size (896 Mb) and 12 chromosomes."])
    write("Figure_2", [
        "# Figure 2",
        "**Genome annotation composition.**",
        f"Repeat classes (inner ring), their subfamilies (outer ring), and gene space (outer track) as percentages of "
        f"the 885.4 Mb assembly (GC 39.7%), of which {rep:.1f}% is repetitive."])
    write("Figure_6", [
        "# Figure 6",
        "**Alleles of the genes that separate the B-type and A-type genomes.**",
        "(A) CDS allele and AF of each haplotype at the six genes where Gwen, Hass, and Puebla each carry an allele "
        "absent from Thille and Unnamed Lowland (filled), with its amino-acid changes and type A versus type B Fst "
        "ranks (mean Fst). Haplotypes are not phased between genes, so the A-shared allele is drawn second.",
        f"(B-D) Allele matrices on Thille: {MATRIX_SHARED}",
        f"(B) {gene_id('SDMYB')} and {gene_id('AGAMOUS-like')}.",
        f"(C) {gene_id('CSU2')} and {gene_id('NPR5-BOP2')}.",
        f"(D) {gene_id('SGRL')} and {gene_id('LCYB')}."])
    write("Figure_7", [
        "# Figure 7",
        "**Coverage at the genes with an A-shared allele.**",
        "(A-C) Gwen CLR, Thille HiFi, Unnamed Lowland HiFi, and Puebla HiFi reads on Thille (top to bottom). "
        + COV_SHARED + " Chromosome medians are given in that order.",
        f"(A) {cov(cpath['SDMYB'])} and {cov(cpath['AGAMOUS-like'])}.",
        f"(B) {cov(cpath['CSU2'])} and {cov(cpath['NPR5-BOP2'])}.",
        f"(C) {cov(cpath['SGRL'])} and {cov(cpath['LCYB'])}."])
    write("Figure_S5", [
        "# Supplementary Figure S5",
        "**Allele matrices of the ten genes recovered by the combined Fst screen.**",
        f"Drawn as in Figure 6B, left to right: " + "; ".join(gene_id(g) for g in S5) + "."])
    write("Figure_S6", [
        "# Supplementary Figure S6",
        "**Coverage at the ten genes recovered by the combined Fst screen.**",
        "Drawn as in Figure 7, on Thille (first ten panels; Gwen, Thille, Unnamed Lowland, Puebla) and on Unnamed "
        "Lowland (last ten; Gwen, Thille, Unnamed Lowland). On Thille: "
        + "; ".join(cov(ft / f"{v(g)}_100kb_stacked_notes.tsv") for g in S5) + ". On Unnamed Lowland: "
        + "; ".join(cov(wt / f"{v(g)}_100kb_stacked_notes.tsv") for g in S5) + "."])
    write("Figure_S7", [
        "# Supplementary Figure S7",
        "**Read coverage at SDMYB on the Hass (UoQ) assembly.**",
        f"Gwen CLR, Thille HiFi, and Unnamed Lowland HiFi reads, drawn as in Figure 7; "
        f"{cov(NT / 'fst_coverage/hass_reference/SDMYB_100kb_stacked_notes.tsv')}. Hass reads are not shown on their "
        "own assembly because the HiFi reads available from NCBI are a fraction of those reported (Nath et al. 2022)."])
    write("Figure_S8", [
        "# Supplementary Figure S8",
        "**Genome-wide read coverage on Thille.**",
        "Median depth of 100 kb bins divided by the chromosome median for Gwen CLR, Thille HiFi, and Unnamed Lowland "
        "HiFi reads; dashed lines mark the class boundaries (0.25, 0.75, 1.25, and 1.5) and values above 2.0 are drawn "
        "at 2.0."])
    write("Figure_S11", [
        "# Supplementary Figure S11",
        "**Allele matrices of the other candidate flowering and pigment genes.**",
        "Drawn as in Figure 6B, left to right: " + "; ".join(gene_id(g) for g in S11) + "."])
    write("Figure_S13", [
        "# Supplementary Figure S13",
        "**Coverage at the other candidate flowering genes.**",
        "Drawn as in Figure 7, left to right: "
        + "; ".join(cov(ft / f"{v(g)}_100kb_stacked_notes.tsv") for g in MOVED_COVERAGE) + "."])
    tn = sorted((NT / "tandem_arrays").glob("*_notes.tsv"))
    s14 = {p.name.split("_")[0].split("-")[0]: p for p in tn}

    def arr(key):   # copy count, span and medians that were in the tandem figure title and panel titles
        t, _, _ = medians(s14[key])
        m = re.search(r"\((\d+) copies; span B ([\d,]+) bp\)", t)
        meds = ", ".join(f"{r['chrom_median']}x" for r in csv.DictReader(open(s14[key]), delimiter="\t"))
        return f"{m.group(1)} copies over {m.group(2)} bp; medians {meds}"
    write("Figure_S14", [
        "# Supplementary Figure S14",
        "**Selected tandem arrays.**",
        "(A-D) Coverage on Thille drawn as in Figure 7, in the smallest fixed window that holds the array, with the "
        "array copies shaded and drawn as a separate track.",
        f"(A) CCR1, chromosome 2, 100 kb ({arr('CCR1')}).",
        f"(B) KCS11, chromosome 6, 500 kb ({arr('KCS11')}).",
        f"(C) The leucine-rich repeat receptor kinase array containing MIK2, chromosome 12, 500 kb ({arr('MIK2')}).",
        f"(D) The interleaved SCPDL (pink) and HST (fuchsia) arrays, chromosome 3, 100 kb ({arr('SCPDL')})."])


if __name__ == "__main__":
    main()
