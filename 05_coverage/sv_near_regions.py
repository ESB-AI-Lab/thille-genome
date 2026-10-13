#!/usr/bin/env python3
"""List every SV within +-flank of the chr3/chr10 deletions and of each flowering gene, in Thille coordinates.

Events: SyRI 1.7.0 Hass-vs-Thille and West-Indian-vs-Thille calls (sv/thille.sv.all.bed; all sizes and types:
NOTAL, HDR, DEL, INS, CPL, CPG, DUP, INVDP, INV, TRANS, INVTR, TDM) and the Gwen read-coverage events >= 100 kb
(also in sv/thille.sv.all.bed, partner_genome gwen). SNPs are not included.
Distance is 0 when the event overlaps the target, otherwise the gap in bp to the nearer target edge.

Usage (from the repo root):
    python3 scripts/sv_near_regions.py [--flank 100000]
Writes summaries/sv_near_regions.tsv (one row per target x event) and summaries/sv_near_regions_summary.tsv
(counts and bp per target x genome x type).
"""

import argparse

import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from coverage_class_summary import read_bed, split_name  # noqa: E402

DELETIONS = [("chr3_deletion", "chr3", 9098588, 10196806), ("chr10_deletion", "chr10", 12268718, 13566546)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--flank", type=int, default=100_000)
    a = p.parse_args()

    # thille.sv.all.bed already holds the Gwen read events (partner_genome gwen), so it is the only input
    ev = pd.read_csv("sv/thille.sv.all.bed", sep="\t", dtype={"partner_chrom": str})
    targets = list(DELETIONS)
    seen = set()
    for name, contig, start, stop in read_bed("flowering_beds/igv/flowering_genes_thille.bed"):
        if (contig, start, stop) in seen:  # Gwen genes sharing one Thille model: one target per locus
            continue
        seen.add((contig, start, stop))
        f = split_name(name)
        targets.append((f"{f[0]}|{f[1]}|{f[2]}", contig, start, stop))

    rows = []
    for tname, contig, start, stop in targets:
        lo, hi = start - a.flank, stop + a.flank
        near = ev[(ev["#chrom"] == contig) & (ev["end"] >= lo) & (ev["start"] + 1 <= hi)]
        for r in near.itertuples(index=False):
            s1, e1 = r.start + 1, r.end
            dist = 0 if (e1 >= start and s1 <= stop) else (start - e1 if e1 < start else s1 - stop)
            rows.append({"target": tname, "target_locus": f"{contig}:{start}-{stop}", "genome": r.partner_genome,
                         "event": r.name, "type": r.type, "size_bp": r.size_bp, "event_locus": f"{contig}:{s1}-{e1}",
                         "distance_bp": dist, "partner": f"{r.partner_chrom}:{r.partner_start}-{r.partner_end}",
                         "detail": r.detail})
    df = pd.DataFrame(rows)
    df.to_csv("summaries/sv_near_regions.tsv", sep="\t", index=False)
    summ = (df.assign(overlaps=df.distance_bp == 0)
              .groupby(["target", "target_locus", "genome", "type"])
              .agg(events=("event", "size"), events_overlapping=("overlaps", "sum"), total_bp=("size_bp", "sum"),
                   largest_bp=("size_bp", "max"))
              .reset_index())
    summ.to_csv("summaries/sv_near_regions_summary.tsv", sep="\t", index=False)
    print(f"{len(targets)} targets, {len(df)} target-event rows", file=sys.stderr)


if __name__ == "__main__":
    main()
