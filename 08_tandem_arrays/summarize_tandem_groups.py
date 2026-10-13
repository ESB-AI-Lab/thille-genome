#!/usr/bin/env python3
"""Summarise the tandem-copy groups per assembly: copy counts and sizes, array span, gaps, figure window,
copy status (BED itemRgb), id/r from the copy names, overlaps between groups, and function where known.

Function comes only from flowering_beds/genes.tsv (SwissProt best hit of the chosen model); the jg / ver1-only
groups have no function data in the local files.

Copies shorter than DUBIOUS_BELOW bp get status "dubious (<1 kb)" (Edwin, 2026-09-26: most genes are over 1 kb,
so a sub-kb copy is unlikely to be a gene copy); the BED's own colour status is kept in bed_status.
Each copy is also checked against the gene annotation of its assembly (Thille: Tiberius GTF; Hass: GenBank
GCA_029852735.1 GTF, gene spans only) and, for Thille, against the RepeatMasker TE annotation.

Usage (from the repo root):
    python3 scripts/summarize_tandem_groups.py    # writes summaries/tandem_groups.tsv, summaries/tandem_copies.tsv
"""

import collections
import csv
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_gene_stacks as p  # noqa: E402  (window and tier math are shared with the figures)

BEDS = "flowering_beds/igv/"
DUBIOUS_BELOW = 1000
GENE_GTF = {"b": "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf",
            "a": "annotation/GCA_029852735.1_ASM2985273v1_genes.gtf"}
TE_GFF = {"b": "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.TE.gff"}


def covered_bp(start, stop, intervals):
    """Number of bases in [start, stop] covered by the union of (s, e) intervals (1-based inclusive)."""
    clipped = sorted((max(s, start), min(e, stop)) for s, e in intervals if s <= stop and e >= start)
    total, cur_s, cur_e = 0, None, None
    for s, e in clipped:
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def gene_model_median(key, copies, models):
    lengths = [m.stop - m.start + 1 for c in copies for m in models[key].get(c.contig, [])
               if m.start <= c.stop and m.stop >= c.start]
    return int(statistics.median(lengths)) if lengths else ""


def te_cover_median(key, copies, tes):
    """Median over copies of the fraction of the copy inside annotated TEs; empty when the assembly has no TE file."""
    if key not in tes:
        return ""
    fractions = []
    for c in copies:
        part = tes[key].get(c.contig)
        near = part[(part["start"] <= c.stop) & (part["stop"] >= c.start)] if part is not None else None
        bp = covered_bp(c.start, c.stop, zip(near["start"], near["stop"])) if near is not None else 0
        fractions.append(bp / (c.stop - c.start + 1))
    return round(statistics.median(fractions), 3)


def copy_status(copy, size):
    return f"dubious (<{DUBIOUS_BELOW // 1000} kb)" if size < DUBIOUS_BELOW else STATUS.get(copy.rgb, str(copy.rgb))
STATUS = {(90, 110, 160): "unresolved", (230, 120, 20): "flowering", (200, 0, 0): "curated/paralog"}
ASSEMBLY = {"b": "Thille", "a": "Hass"}


def main():
    names = {"b": p.ContigNames(), "a": p.ContigNames("hass_cm_chromosomes.tsv")}
    refs = []
    for key, asm in (("b", "thille"), ("a", "hass")):
        ref = p.Reference(key, p.read_bed(BEDS + f"flowering_genes_{asm}.bed", 2), [], [], names[key])
        ref.tandem = p.read_tandem(BEDS + f"tandem_copies_{asm}.bed")
        refs.append(ref)
    display = p.display_names(set(refs[0].genes) | set(refs[1].genes), 2)
    figures = {f.key: f for f in p.tandem_figures(refs, display, [])}
    genes = {f"{row['family']}|{row['ver1_id']}": row
             for row in csv.DictReader(open("flowering_beds/genes.tsv"), delimiter="\t")}

    contigs = {ref.key: {c for c in ref.tandem} for ref in refs}
    models = {key: p.read_gene_models(GENE_GTF[key], contigs[key]) for key in GENE_GTF}
    tes = {}
    for key, path in TE_GFF.items():
        tes[key], _counts = p.read_te(path, contigs[key])

    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for ref in refs:
        for copies in ref.tandem.values():
            for copy in copies:
                groups[copy.group][ref.key].append(copy)
    groups_at = collections.defaultdict(set)  # identical copy coordinates -> groups using them
    for group, by_ref in groups.items():
        for key, copies in by_ref.items():
            for c in copies:
                groups_at[(key, c.contig, c.start, c.stop)].add(group)

    rows, copy_rows = [], []
    for group in sorted(groups):
        figure = figures[group]
        for key in ("b", "a"):
            copies = sorted(groups[group].get(key, []), key=lambda c: c.start)
            if not copies:
                continue
            sizes = [c.stop - c.start + 1 for c in copies]
            gaps = [b.start - a.stop - 1 for a, b in zip(copies, copies[1:])]
            ids = [float(m.group(1)) for c in copies for m in [re.search(r"\|id([\d.]+)", c.name)] if m]
            rs = [float(m.group(1)) for c in copies for m in [re.search(r"\|r([\d.]+)", c.name)] if m]
            status = collections.Counter(copy_status(c, size) for c, size in zip(copies, sizes))
            bed_status = collections.Counter(STATUS.get(c.rgb, str(c.rgb)) for c in copies)
            with_model = sum(any(m.start <= c.stop and m.stop >= c.start for m in models[key].get(c.contig, []))
                             for c in copies)
            coords = sorted((c.contig, c.start, c.stop) for c in copies)
            shared = sorted({o for c in copies for o in groups_at[(key, c.contig, c.start, c.stop)]} - {group})
            identical = [o for o in shared
                         if sorted((c.contig, c.start, c.stop) for c in groups[o].get(key, [])) == coords]
            lo, hi = p.window(figure.regions[key], figure.size)
            gene = genes.get(group, {})
            hit = gene.get("thille_model_best_hit" if key == "b" else "hass_model_best_hit", "")
            end = max(c.stop for c in copies)
            rows.append({
                "group": group, "display": figure.display, "flowering_family": gene.get("family", ""),
                "category": gene.get("category", ""), "assembly": ASSEMBLY[key],
                "chromosome": names[key].label(copies[0].contig), "contig": copies[0].contig,
                "n_copies": len(copies), "copy_size_min": min(sizes), "copy_size_median": int(statistics.median(sizes)),
                "copy_size_max": max(sizes), "array_start": copies[0].start, "array_end": end,
                "array_span": end - copies[0].start + 1,
                "gap_median": int(statistics.median(gaps)) if gaps else "", "gap_max": max(gaps) if gaps else "",
                "window_tier_kb": figure.size // 1000, "window": f"{lo}-{hi}",
                "copies_ge_1kb": sum(size >= DUBIOUS_BELOW for size in sizes),
                "copies_dubious_lt_1kb": status.get(f"dubious (<{DUBIOUS_BELOW // 1000} kb)", 0),
                "copies_unresolved": bed_status.get("unresolved", 0), "copies_flowering": bed_status.get("flowering", 0),
                "copies_curated_or_paralog": bed_status.get("curated/paralog", 0),
                "copies_overlapping_gene_model": with_model,
                "median_gene_model_bp": gene_model_median(key, copies, models),
                "median_copy_TE_cover": te_cover_median(key, copies, tes),
                "identity_min": min(ids) if ids else "", "identity_max": max(ids) if ids else "",
                "r_min": min(rs) if rs else "", "r_max": max(rs) if rs else "",
                "identical_copy_set_as": ",".join(identical),
                "shares_copies_with": ",".join(o for o in shared if o not in identical),
                "function_best_hit": hit or ("(no function data in local files)" if not gene else ""),
                "copy_sizes": ",".join(map(str, sizes))})
            for c, size in zip(copies, sizes):
                hits = [m for m in models[key].get(c.contig, []) if m.start <= c.stop and m.stop >= c.start]
                te_part = tes.get(key, {}).get(c.contig)
                te_bp = ""
                if te_part is not None:
                    near = te_part[(te_part["start"] <= c.stop) & (te_part["stop"] >= c.start)]
                    te_bp = covered_bp(c.start, c.stop, zip(near["start"], near["stop"]))
                copy_rows.append({"group": group, "assembly": ASSEMBLY[key], "contig": c.contig, "start_1b": c.start,
                                  "end": c.stop, "size": size, "strand": c.strand,
                                  "status": copy_status(c, size), "bed_status": STATUS.get(c.rgb, str(c.rgb)),
                                  "gene_models": ",".join(f"{m.gene_id}({m.stop - m.start + 1}bp,{m.strand})" for m in hits),
                                  "copy_bp_in_gene_models": covered_bp(c.start, c.stop, [(m.start, m.stop) for m in hits]),
                                  "copy_bp_in_TEs": te_bp,
                                  "name": c.name})

    Path("summaries").mkdir(exist_ok=True)
    for path, data in (("summaries/tandem_groups.tsv", rows), ("summaries/tandem_copies.tsv", copy_rows)):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(data)
        print(f"{path}: {len(data)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
