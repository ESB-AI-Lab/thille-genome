#!/usr/bin/env python3
"""Write the 152 SyRI inversions (> 100 bp; Thille vs Hass UoQ, Unnamed Lowland and Puebla) as BED files with their
raw-read support calls.

  inversions.thille_coords.bed            all 152 on Thille v4 coordinates
  inversions.<pair>_coords.bed            the inversions of one pair on that genome's original GenBank coordinates
                                          (pair = hass, westindian [Unnamed Lowland], puebla)

Inputs: inversions.tsv (coordinates, already converted back to the original accession and strand and checked by
alignment) and support_summary_v2.tsv (calls from summarize_support.py). Coordinates are BED (0-based start, end).
Columns (BED6+): chrom start end name score strand(-) pair partner_chrom partner_start partner_end thille_len
other_len, then one column per read set and assembly: <reads>_on_<thille|other> = call (spanning/split read counts),
"." where that read set was not aligned to that assembly.

Usage (from compare/v4/inversions): make_inversion_beds.py [outdir]
"""
import csv
import os
import sys

READS = ["thille_hifi", "lowland_hifi", "gwen_clr", "puebla_hifi"]
SUPPORT_COLS = [f"{r}_on_{side}" for side in ("thille", "other") for r in READS]
HEADER = ["#chrom", "start", "end", "name", "score", "strand", "pair", "partner_chrom", "partner_start",
          "partner_end", "thille_len", "other_len", *SUPPORT_COLS]
ORDER = {f"chr{i}": i for i in range(1, 13)}


def call(s, col):
    c = s.get(f"{col}_call", "")
    if not c:
        return "."
    return f"{c} ({s[f'{col}_span']}/{s[f'{col}_split']})"


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    inv = list(csv.DictReader(open("inversions.tsv"), delimiter="\t"))
    sup = {r["inv_id"]: r for r in csv.DictReader(open("support_summary_v2.tsv"), delimiter="\t")}
    missing = [r["inv_id"] for r in inv if r["inv_id"] not in sup]
    if missing:
        sys.exit(f"no support row for {len(missing)} inversions, e.g. {missing[:3]}")
    thille, per_pair = [], {}
    for r in inv:
        s = sup[r["inv_id"]]
        calls = [call(s, c) for c in SUPPORT_COLS]
        common = [r["inv_id"], "0", "-", r["pair"]]
        thille.append([r["thille_chrom"], r["thille_start"], r["thille_end"], *common,
                       r["other_chrom"], r["other_start"], r["other_end"], r["thille_len"], r["other_len"], *calls])
        per_pair.setdefault(r["pair"], []).append(
            [r["other_chrom"], r["other_start"], r["other_end"], *common,
             r["thille_chrom"], r["thille_start"], r["thille_end"], r["thille_len"], r["other_len"], *calls])
    files = {"inversions.thille_coords.bed": sorted(thille, key=lambda x: (ORDER[x[0]], int(x[1])))}
    for p, rows in per_pair.items():
        files[f"inversions.{p}_coords.bed"] = sorted(rows, key=lambda x: (x[0], int(x[1])))
    for name, rows in files.items():
        path = os.path.join(out, name)
        with open(path, "w") as fh:
            fh.write("\t".join(HEADER) + "\n")
            fh.writelines("\t".join(x) + "\n" for x in rows)
        print(f"{path}: {len(rows)} inversions")


if __name__ == "__main__":
    main()
