#!/usr/bin/env python3
"""Build the four gene tables (flowering genes with Fst peaks, flowering genes without Fst peaks, genes with coverage
differences from the Gwen paper, tandem arrays) as TSV and JSON for the Word document.

Sources:
  flowering_beds/genes.tsv                 curated genes, categories, Gwen IDs, Thille models, curated best hits
  flowering_beds/igv/flowering_genes_thille.bed   Thille locus for genes without a Thille gene model
  summaries/tandem_arrays.tsv, tandem_array_members.tsv   supported Thille tandem arrays
  summaries/gene_model_functions.tsv       Swiss-Prot best hit per Thille model (scripts/swissprot_homology.py)
  db/swissprot/swissprot_meta.tsv.gz       UniProtKB/Swiss-Prot 2026_03 names
Function text is UniProt "Function [CC]" (first sentence, evidence tags removed).

Usage (from the repo root):
    python3 scripts/make_gene_tables.py     # writes summaries/tables/*.tsv and summaries/tables/gene_tables.json
"""

import collections
import csv
import gzip
import json
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import swissprot_homology as sh  # noqa: E402  (reuses the UniProt function fetch and its cache)

TABLES = [
    ("fst_flowering", "Candidate flowering genes in top 1% Fst windows",
     {"strong", "floral/reproductive/scent/light"}),
    ("no_fst_flowering", "Flowering genes without Fst peaks (selected by coverage)",
     {"annex* promoted (coverage, not Fst)"}),
    ("coverage_gwen", "Genes with coverage differences from the Gwen study",
     {"paper control"}),
]
COLUMNS = ["SwissProt gene", "Thille location", "Thille gene ID", "Gwen gene (ver1 / jg)", "Function"]


def entry_name(best_hit):
    token = best_hit.split(" ")[0] if best_hit else ""
    return token if re.fullmatch(r"[A-Z0-9]+_[A-Z0-9]+", token) else ""


def location(span_0b):
    """'chr3:81996880-81999600(-)' (0-based start) -> 'chr3:81,996,881-81,999,600 (-)'."""
    m = re.fullmatch(r"(\S+):(\d+)-(\d+)\(([+-])\)", span_0b or "")
    if not m:
        return ""
    return f"{m.group(1)}:{int(m.group(2)) + 1:,}-{int(m.group(3)):,} ({m.group(4)})"


def short_function(protein_name, cc_function):
    text = re.sub(r"\s*\{ECO:[^}]*\}", "", cc_function or "").removeprefix("FUNCTION: ").strip()
    text = re.sub(r"\s*\(PubMed:[^)]*\)", "", text)
    first = re.split(r"(?<=[a-z0-9)\]])\.\s", text, maxsplit=1)[0].rstrip(".")
    name = re.split(r" \(", protein_name or "", maxsplit=1)[0]
    if first:
        return f"{name}. {first}." if name else f"{first}."
    return f"{name}. No UniProt function annotation." if name else "No SwissProt hit."


def main():
    genes = list(csv.DictReader(open("flowering_beds/genes.tsv"), delimiter="\t"))
    bed = {}
    for line in open("flowering_beds/igv/flowering_genes_thille.bed"):
        if line.startswith("chr"):
            f = line.rstrip("\n").split("\t")
            bed["|".join(f[3].split("|")[:2])] = f"{f[0]}:{int(f[1])}-{f[2]}({f[5]})"
    functions = {r["model"]: r for r in csv.DictReader(open("summaries/gene_model_functions.tsv"), delimiter="\t")}

    by_entry, by_accession = {}, {}
    with gzip.open(sh.SPROT_META, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_entry[row["Entry Name"]] = row
            by_accession[row["Entry"]] = row

    curated = {entry_name(g["thille_model_best_hit"]) or entry_name(g["gwen_model_best_hit"]) for g in genes}
    curated.discard("")
    accessions = {by_entry[e]["Entry"] for e in curated if e in by_entry}
    accessions |= {r["swissprot_accession"] for r in functions.values() if r["swissprot_accession"]}
    cc = sh.fetch_cc_function(accessions)

    out = {}
    for key, title, categories in TABLES:
        rows = []
        for g in genes:
            if g["category"] not in categories:
                continue
            entry = entry_name(g["thille_model_best_hit"])
            source = "curated (Thille model)"
            if not entry:
                entry, source = entry_name(g["gwen_model_best_hit"]), "curated (Gwen model)"
            if not entry and g["thille_model"] in functions and functions[g["thille_model"]]["swissprot_entry"]:
                entry, source = functions[g["thille_model"]]["swissprot_entry"], "homology"
            meta = by_entry.get(entry, {})
            gene = meta.get("Gene Names (primary)", "") or (entry.split("_")[0] if entry else "")
            label = "" if g["label"] in (entry, gene) else f"{g['label']}: "
            suffix = " (homology)" if source == "homology" else ""
            model = g["thille_model"]
            span = g["thille_model_span_0b"] or bed.get(f"{g['family']}|{g['ver1_id']}", "")
            note = ""
            if not re.fullmatch(r"g\d+", model):
                note = "; aligned locus" if model == "unannotated" else ""
                model = ("none (unannotated)" if model == "unannotated"
                         else "none (" + model.replace("absent (closest paralog", "ortholog absent; closest paralog").rstrip(")") + ")")
            rows.append({
                "SwissProt gene": f"{label}{gene} ({entry}){suffix}" if entry else f"{g['label']}: no SwissProt hit",
                "Thille location": location(span) + note,
                "Thille gene ID": model,
                "Gwen gene (ver1 / jg)": f"{g['ver1_id']} / {g['contig_genes'].replace(';', ', ')}" if g["contig_genes"] else g["ver1_id"],
                "Function": short_function(meta.get("Protein names", ""), cc.get(meta.get("Entry", ""), "")),
                "_label": g["label"], "_category": g["category"], "_hit_from": source,
            })
        out[key] = {"title": title, "columns": COLUMNS, "rows": rows}

    arrays = list(csv.DictReader(open("summaries/tandem_arrays.tsv"), delimiter="\t"))
    members = list(csv.DictReader(open("summaries/tandem_array_members.tsv"), delimiter="\t"))
    rows = []
    for a in arrays:
        ms = [m for m in members if m["array"] == a["array"] and m["in_ge_70pct_pair"] == "True"]
        genes_hit = collections.Counter(functions[m["model"]]["swissprot_gene"] for m in ms if functions[m["model"]]["swissprot_gene"])
        top_gene, top_n = genes_hit.most_common(1)[0] if genes_hit else ("", 0)
        top_row = next((functions[m["model"]] for m in ms if functions[m["model"]]["swissprot_gene"] == top_gene), None) if top_gene else None
        others = ", ".join(f"{g} ({n})" for g, n in genes_hit.most_common()[1:])
        name = (f"{top_gene} ({top_row['swissprot_entry']}; {top_n} of {len(ms)} models)" if top_row
                else "no SwissProt hit")
        if others:
            name += f"; also {others}"
        curated_here = sorted({f"{g['label']} ({g['thille_model']})" for g in genes
                               for m in ms if g["thille_model"] == m["model"]})
        if curated_here:
            name += "; contains curated " + ", ".join(curated_here)
        contig = ms[0]["contig"]
        start, end = min(int(m["start"]) for m in ms), max(int(m["end"]) for m in ms)
        groups = [g.split("|") for g in a["groups"].split(",")]
        gwen = ", ".join(sorted({p[1] if p[0] != p[1] else p[0] for p in groups}))
        rows.append({
            "SwissProt gene": name,
            "Thille location": f"{contig}:{start:,}-{end:,} ({'/'.join(sorted({m['strand'] for m in ms}))})",
            "Thille gene ID": ", ".join(m["model"] for m in ms),
            "Gwen gene (ver1 / jg)": gwen,
            "Function": short_function(top_row["protein_name"], top_row["function"]) if top_row else "No SwissProt hit.",
            "Copies": len(ms),
            "Mean copy size (bp)": f"{statistics.mean(int(m['gene_bp']) for m in ms):,.0f}",
            "_groups": a["groups"],
        })
    out["tandem_arrays"] = {"title": "Supported tandem arrays in Thille (gene models with >= 70% pairwise exon identity)",
                            "columns": COLUMNS + ["Copies", "Mean copy size (bp)"], "rows": rows}

    folder = Path("summaries/tables")
    folder.mkdir(parents=True, exist_ok=True)
    for key, table in out.items():
        with open(folder / f"{key}.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table["rows"][0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(table["rows"])
    (folder / "gene_tables.json").write_text(json.dumps(out, indent=1))
    for key, table in out.items():
        print(f"{key}: {len(table['rows'])} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
