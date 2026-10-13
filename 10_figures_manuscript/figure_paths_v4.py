#!/usr/bin/env python3
"""Absolute paths of every manuscript figure panel and data table on the Thille v4 assembly, in manuscript order.

Panels follow manuscript/build_scripts/build_latex.py (FIGURES, F6, F7, S5, S11). Assembly-level figures come from the
thille-hic-bf session; --relabel points at its folder of figures relabelled "Unnamed Lowland".
Usage (from the repo root): python3 scripts/figure_paths_v4.py [--relabel <dir>]
Writes manuscript/figure_paths_v4.tsv (item, panel, path, exists) and prints missing paths to stderr.
"""
import argparse
import sys
from pathlib import Path

REPO = Path("/gpool/edwin/avocados/thille_paper")
TH = Path("/epool/edwin/pamericana/avoflow/assemblies/thille.hic")
AV = Path("/epool/edwin/pamericana/avoflow")
R4 = REPO / "Results_v4"
FA = REPO / "v4work/Results_fst_alleles"
FT, FW, FH = R4 / "fst_coverage/thille_reference", R4 / "fst_coverage/wi_reference", R4 / "fst_coverage/hass_reference"
PT = R4 / "pigment_coverage/thille_reference"
# No figure text inside the panels (Edwin, 2026-10-09): every panel points at its --no-text render; titles, notes and
# explanatory numbers are in manuscript/figure_legends/Figure_<id>.md. The earlier renders are kept where they were.
R4N = REPO / "Results_v4_notext"
FTN, FWN, FHN = R4N / "fst_coverage/thille_reference", R4N / "fst_coverage/wi_reference", R4N / "fst_coverage/hass_reference"
PTN = R4N / "pigment_coverage/thille_reference"
NT = TH / "compare/v4/notext"   # thille-hic-bf no-text panels (plan v18)
V3_TO_V4 = {"VOZ1_g24979": "VOZ1_g22280", "VOZ1_g26946": "VOZ1_g26844"}  # labels that carry a gene ID
# Layout (Edwin, 2026-10-07): Figure 6 is the allele figure of the six genes that separate the B-type genomes (Thille,
# Unnamed Lowland) from the A-type ones (Gwen, Hass, Puebla): A = summary table of alleles per haplotype with AF, then
# one matrix per gene; Figure 7 shows coverage of the same six genes. The other matrices join S11, the other coverage
# figures form S13, the tandem arrays (former Figure 8) S14, and the ModDotPlot comparisons S15.
B6 = ["SDMYB", "AGAMOUS-like", "CSU2", "NPR5-BOP2", "SGRL", "LCYB"]
MOVED_MATRICES = ["SPA1", "VOZ1_g24979", "PLA2-LHD2", "PAR1", "NRX1", "CCR1_g514", "CCR1_g512", "PSY", "CCS-like"]
MOVED_COVERAGE = ["SPA1", "VOZ1_g24979", "PLA2-LHD2", "PAR1", "NRX1", "CCR1_g514", "CCR1_g512"]
S11 = MOVED_MATRICES + ["ALP1", "AP2", "ARF2A-like", "EMF1", "HY5", "ICR1", "KCS11", "MIK2", "MIK2_g15851", "NAC25", "SGR1",
                        "TIC", "VOZ1_g26946", "ZFP2"]
S5 = ["FBD-RNI-like", "RER3", "RPT3", "PBL34", "Ole-e-10-like", "LIP1", "ZEP", "CIGR2", "CCCH-zinc-finger", "CTR1"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--relabel", default=str(TH / "compare/v4/relabel2"))
    a = p.parse_args()
    RL = Path(a.relabel)
    v = lambda g: V3_TO_V4.get(g, g)
    rows = [
        ("Figure 1", "", R4N / "figure1/thille_cumulative_length_notext.pdf"),
        ("Figure 2", "", R4N / "figure2/thille_genome_composition_notext.pdf"),
        ("Figure 3", "A", TH / "merqury/thille_hifi_v4_k21/thille_hifi_v4.thille.hifiv4.asm.hic.scaffold.masked.chr.only.spectra-cn.st.pdf"),
        ("Figure 3", "B", NT / "Fig3B_busco.pdf"),
        ("Figure 4", "A raw and B balanced", NT / "all_genome_500kb_V2_strict-filters_notext.pdf"),
        ("Figure 5", "A Unnamed Lowland", NT / "Fig5A_Unnamed_Lowland.pdf"),
        ("Figure 5", "A Hass", NT / "Fig5A_Hass.pdf"),
        ("Figure 5", "A Puebla", NT / "Fig5A_Puebla.pdf"),
        ("Figure 5", "B", RL / "westindian_thille_hass_puebla_plotsr.pdf"),
        ("Figure 5", "C", NT / "karyotype/thille_karyotype_satellitefirst_split-idiogram_lavender_flowering.pdf"),
    ]
    rows += [("Figure 6", "A summary of alleles", FA / "a_shared_alleles_notext.pdf")]
    rows += [("Figure 6", g, FA / f"{g}_variant_matrix_notext.pdf") for g in B6]
    rows += [("Figure 7", g, (PTN if g in ("SGRL", "LCYB") else FTN) / f"{g}_100kb_stacked.pdf") for g in B6]
    rows += [("Figure S2", "", NT / "S2_Gwen.pdf")]
    rows += [("Figure S5", v(g), FA / f"{v(g)}_variant_matrix_notext.pdf") for g in S5]
    rows += [("Figure S6", f"Thille {g}", FTN / f"{v(g)}_100kb_stacked.pdf") for g in S5]
    rows += [("Figure S6", f"Unnamed Lowland {g}", FWN / f"{v(g)}_100kb_stacked.pdf") for g in S5]
    rows += [("Figure S7", "SDMYB on Hass UoQ", FHN / "SDMYB_100kb_stacked.pdf")]
    rows += [("Figure S8", "", REPO / "v4work/Results_genome_coverage/genome_wide_coverage_thille_notext.pdf")]
    rows += [("Figure S9", "A", NT / "S9_R3_tracks.pdf"), ("Figure S9", "B", NT / "S9_R10_tracks.pdf")]
    rows += [("Figure S10", "A", RL / "R3_fourway_ribbon.pdf"), ("Figure S10", "B", RL / "R10_fourway_ribbon.pdf")]
    rows += [("Figure S11", v(g), FA / f"{v(g)}_variant_matrix_notext.pdf") for g in S11]
    rows += [("Figure S13", v(g), FTN / f"{v(g)}_100kb_stacked.pdf") for g in MOVED_COVERAGE]
    rows += [("Figure S14", g, R4N / f"tandem_arrays/{g}_tandem.pdf")
             for g in ("CCR1_100kb", "KCS11_500kb", "MIK2_500kb", "SCPDL-HST_100kb")]
    rows += [("Figure S15", "A Hass", NT / "S15A_Hass.pdf"), ("Figure S15", "B Unnamed Lowland", NT / "S15B_Unnamed_Lowland.pdf"),
             ("Figure S15", "C Puebla", NT / "S15C_Puebla.pdf")]
    S = REPO / "v4work/summaries"
    rows += [
        ("Table S1", "SVs >= 100 kb (SyRI, incl. Puebla; Unnamed Lowland labels)", REPO / "manuscript/supplementary/Table_S1_sv_events_ge100kb.tsv"),
        ("Table 1", "assembly and Hi-C", REPO / "summaries/v4/table1_v4.tsv"),
        ("Table 2", "annotation", REPO / "summaries/tables/annotation_v4.tsv"),
        ("Table 3", "mapping, Thille HiFi", REPO / "summaries/v4/mapping_stats/thille.hifi.thille_v4.full_stats.tsv"),
        ("Table 3", "mapping, Thille HiFi identity", REPO / "summaries/v4/mapping_stats/thille.hifi.thille_v4.identity.tsv"),
        ("Table 3", "mapping, Gwen CLR", REPO / "summaries/v4/mapping_stats/gwen.clr.thille_v4.full_stats.tsv"),
        ("Table 3", "mapping, Gwen CLR identity", REPO / "summaries/v4/mapping_stats/gwen.clr.thille_v4.identity.tsv"),
        ("Table 3", "mapping, Unnamed Lowland HiFi", REPO / "summaries/v4/mapping_stats/westindian.hifi.thille_v4.full_stats.tsv"),
        ("Table 3", "mapping, Unnamed Lowland HiFi identity", REPO / "summaries/v4/mapping_stats/westindian.hifi.thille_v4.identity.tsv"),
        ("Table 3", "mapping, Thille HiFi to full v4 scaffolds", REPO / "summaries/v4/mapping_stats/thille.hifi.thille_v4_full.full_stats.tsv"),
        ("Table 3", "mapping, Hass UoQ 12 chr / full (unchanged)", REPO / "summaries/mapping_stats/rerun_20261002/mapping_summary.tsv"),
        ("Tables 4-6", "gene tables (docx source)", S / "tables/gene_tables.json"),
        ("Table 4", "Fst flowering genes", S / "tables/fst_flowering.tsv"),
        ("Table 5", "flowering genes without Fst", S / "tables/no_fst_flowering.tsv"),
        ("Table 6", "coverage at Gwen genes", S / "tables/coverage_gwen.tsv"),
        ("Table 7", "tandem arrays", S / "tables/tandem_arrays.tsv"),
        ("Table S3", "genome-wide coverage", S / "genome_wide_coverage_thille.tsv"),
        ("Table S4", "A-shared alleles", S / "fst_alleles/a_shared_alleles.tsv"),
        ("Table S5", "AF skew", S / "fst_alleles/af_skew_check.tsv"),
        ("Table S6", "pigment Fst", S / "fst_alleles/pigment_fst.tsv"),
        ("Text", "allele letters and AF, all 39 genes", S / "fst_alleles/allele_summary.tsv"),
        ("Text", "allele class pattern", S / "fst_alleles/allele_class_pattern.tsv"),
        ("Text", "per-gene coverage classes, Thille reference", S / "flowering_gene_coverage_classes_thille.tsv"),
        ("Text", "coverage SV hits, Thille reference", S / "coverage_sv_hits_thille.tsv"),
        ("Text", "coverage class counts per chromosome", S / "coverage_class_counts_thille.tsv"),
        ("Text", "AB window classes", S / "ab_window_classes.tsv"),
        ("Text", "Unnamed Lowland segment check", S / "wi_segment_check.tsv"),
        ("Text", "SDMYB alleles", S / "sdmyb_alleles/variants.tsv"),
        ("Text", "gene model functions (Swiss-Prot)", S / "gene_model_functions.tsv"),
        ("Text", "tandem arrays, members, TE", S / "tandem_array_members.tsv"),
        ("Reference", "v3 to v4 gene ID map", REPO / "summaries/gene_id_map_v3_v4.tsv"),
        ("Reference", "v3 to v4 change list", REPO / "summaries/v3_vs_v4_changes.md"),
    ]
    out = REPO / "manuscript/figure_paths_v4.tsv"
    with open(out, "w") as fh:
        fh.write("item\tpanel\tpath\texists\n")
        for item, panel, path in rows:
            ok = path.is_file()
            fh.write(f"{item}\t{panel}\t{path}\t{'yes' if ok else 'NO'}\n")
            if not ok:
                print(f"missing: {item} {panel} {path}", file=sys.stderr)
    print(f"{out}: {len(rows)} rows")


if __name__ == "__main__":
    main()
