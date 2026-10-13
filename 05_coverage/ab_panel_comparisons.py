#!/usr/bin/env python3
"""Every coverage comparison for the A/B panel genes in one long table: gene x reference x read set.

Reads the per-gene class tables and hit tables for the three references (Thille, Hass UoQ, West Indian).
Columns: class fractions of the gene span, median depth over the gene, chromosome median, and the deletion /
duplication hit (gene bp inside >= 1 kb runs, longest run). Writes summaries/ab_panel_all_comparisons.tsv.
Usage (from the repo root): python3 scripts/ab_panel_comparisons.py
"""
import pandas as pd

GENES = ["SDMYB", "AP2", "NAC25", "MIK2", "SUD1", "Y5713", "OTU7", "APT3", "SPA1", "KCS11", "PAR1", "CCR1"]
REFS = [("Thille", "thille", {"gwen.clr.thille": "Gwen (A)", "thille.hifi.thille": "Thille (B)", "westindian.hifi.thille": "West Indian"}),
        ("Hass UoQ", "hass", {"gwen.clr.hass_uoq": "Gwen (A)", "thille.hifi.hass_uoq": "Thille (B)", "westindian.hifi.hass_uoq": "West Indian"}),
        ("West Indian", "wi", {"gwen.clr.wi": "Gwen (A)", "thille.hifi.wi": "Thille (B)", "westindian.hifi.wi": "West Indian"})]
rows = []
for ref_name, key, tracks in REFS:
    cls = pd.read_csv(f"summaries/flowering_gene_coverage_classes_{key}.tsv", sep="\t")
    hits = pd.read_csv(f"summaries/coverage_sv_hits_{key}.tsv", sep="\t")
    for i, r in cls.iterrows():
        if r.gene not in GENES or (r.gene == "MIK2" and not str(r.gwen_gene).startswith("Chr12")):
            continue
        h = hits.iloc[i]
        assert h.locus == r.locus
        for t, reads in tracks.items():
            hit = lambda s: (f"{h[f'{t}_{s}_gene_bp_in_runs']} bp (run {h[f'{t}_{s}_longest_run_bp']:,})"
                             if h[f"{t}_{s}_hit"] == "yes" else "no")
            rows.append({"gene": r.gene, "gwen_gene": r.gwen_gene, "reference": ref_name, "model": r.model,
                         "locus": r.locus, "gene_bp": r.gene_bp, "reads": reads,
                         **{c: r[f"{t}_{c}"] for c in ("low", "half", "medium", "high", "outlier")},
                         "median_depth": r[f"{t}_median_depth"], "chrom_median": r[f"{t}_chrom_median"],
                         "deletion_hit": hit("deletion"), "duplication_hit": hit("duplication")})
df = pd.DataFrame(rows)
df["gene"] = pd.Categorical(df.gene, GENES, ordered=True)
df = df.sort_values(["gene", "gwen_gene", "reference", "reads"], kind="stable")
df.to_csv("summaries/ab_panel_all_comparisons.tsv", sep="\t", index=False)
