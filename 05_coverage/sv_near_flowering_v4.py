#!/usr/bin/env python3
"""Inversions and deletions near the flowering genes on Thille v4, by distance class (Edwin, 2026-10-12).

SVs: SyRI events on Thille v4 coordinates (compare/v4/sv_bed/thille.sv.all.bed; comparisons with Hass (UoQ), Unnamed
Lowland and Puebla), types INV and DEL, at least 50 bp. Genes: flowering_beds/v4/igv/flowering_genes_thille.bed
(name field = label|Gwen ID|Thille model|...), with categories from flowering_beds/v4/genes.tsv.
Distance classes, per gene and SV (nearest edges): body = overlaps the gene span; cis = 1 bp to 10 kb from the gene;
distal = 10 to 100 kb. An SV is counted once per gene, in its closest class. Positional "trans" effects cannot be
assigned from coordinates and are not counted.
Also writes the length and count of all 152 inversions (compare/v4/inversions/inversions.tsv) per comparison.
Usage (repo root): python3 scripts/sv_near_flowering_v4.py   # -> summaries/v4/sv_near_flowering/
"""
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

SV = "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/compare/v4/sv_bed/thille.sv.all.bed"
INV = "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/compare/v4/inversions/inversions.tsv"
BED = "flowering_beds/v4/igv/flowering_genes_thille.bed"
GENES = "flowering_beds/v4/genes.tsv"
OUT = Path("summaries/v4/sv_near_flowering")
PAIRS = {"hass": "Hass (UoQ)", "westindian": "Unnamed Lowland", "puebla": "Puebla"}
CIS, DISTAL = 10_000, 100_000


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cat = {}
    for r in csv.DictReader(open(GENES), delimiter="\t"):
        cat.setdefault(r["label"], r["category"])
    genes = []
    for line in open(BED):
        if line.startswith("track"):
            continue
        f = line.rstrip("\n").split("\t")
        label = f[3].split("|")[0]
        genes.append((label, f[3].split("|")[2], f[0], int(f[1]), int(f[2])))
    svs = defaultdict(list)
    for r in csv.DictReader(open(SV), delimiter="\t"):
        if r["type"] in ("INV", "DEL") and r["partner_genome"] in PAIRS and int(r["size_bp"]) >= 50:
            svs[r["#chrom"]].append((int(r["start"]), int(r["end"]), r["type"], r["partner_genome"], int(r["size_bp"])))
    rows, tot = [], Counter()
    for label, model, c, s, e in genes:
        counts = Counter()
        for a, b, typ, pair, size in svs.get(c, []):
            dist = 0 if (a < e and b > s) else min(abs(a - e), abs(s - b))
            cls = "body" if dist == 0 else "cis" if dist <= CIS else "distal" if dist <= DISTAL else None
            if cls:
                counts[(typ, cls, pair)] += 1
                tot[(cat.get(label, ""), typ, cls)] += 1
        row = {"gene": label, "thille_model": model, "category": cat.get(label, ""), "locus": f"{c}:{s + 1}-{e}"}
        for typ in ("INV", "DEL"):
            for cls in ("body", "cis", "distal"):
                row[f"{typ}_{cls}"] = sum(counts[(typ, cls, p)] for p in PAIRS)
                row[f"{typ}_{cls}_by_pair"] = ";".join(f"{p}:{counts[(typ, cls, p)]}" for p in PAIRS if counts[(typ, cls, p)])
        rows.append(row)
    with open(OUT / "sv_near_flowering_genes.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    # inversions overall
    inv = list(csv.DictReader(open(INV), delimiter="\t"))
    with open(OUT / "inversion_lengths.tsv", "w") as fh:
        fh.write("comparison\tcount\ttotal_thille_bp\tmin_bp\tmedian_bp\tmax_bp\tge_100kb\tge_1Mb\n")
        for pair in list(PAIRS) + ["all"]:
            L = [int(r["thille_len"]) for r in inv if pair == "all" or r["pair"] == pair]
            fh.write(f"{PAIRS.get(pair, 'all')}\t{len(L)}\t{sum(L)}\t{min(L)}\t{int(statistics.median(L))}\t{max(L)}\t"
                     f"{sum(x >= 100_000 for x in L)}\t{sum(x >= 1_000_000 for x in L)}\n")
    print(open(OUT / "inversion_lengths.tsv").read())
    for k in sorted(tot):
        print(*k, tot[k], sep="\t")


if __name__ == "__main__":
    main()
