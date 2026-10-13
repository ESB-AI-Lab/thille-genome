#!/usr/bin/env python3
"""Call coverage SV hits at candidate genes from the plot_gene_stacks.py cache.

Hit definition (Edwin, 2026-09-29): a gene overlaps a contiguous run of >= --min-run bp (default 1,000)
whose positions are all in the deletion set (low or half) or all in the duplication set (high or outlier).
Low and half positions next to each other form one run (likewise high and outlier). Runs are measured on the
chromosome, not clipped to the gene. Classes: M = int(chromosome median), same bins as the figures.

Usage (from the repo root):
    python3 scripts/coverage_sv_hits.py --bed flowering_beds/igv/flowering_genes_thille.bed \\
        --tracks gwen.clr.thille thille.hifi.thille --out summaries/coverage_sv_hits_thille.tsv
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coverage_class_summary import NAME_COLUMNS, class_index, read_bed, split_name  # noqa: E402
import plot_gene_stacks as pgs  # noqa: E402

SETS = {"deletion": (0, 1), "duplication": (3, 4)}  # low, half / high, outlier
FLANK = 1_000_000  # search this far past the gene for the ends of runs that touch it


def runs_overlapping(mask, offset, start, stop):
    """Runs of True in mask (mask[0] is chromosome position `offset`) that overlap [start, stop]; 1-based inclusive."""
    padded = np.concatenate(([False], mask, [False]))
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    out = []
    for a, b in zip(edges[::2], edges[1::2]):
        rs, re_ = offset + a, offset + b - 1
        if rs <= stop and re_ >= start:
            out.append((rs, re_))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True)
    p.add_argument("--tracks", nargs="+", required=True)
    p.add_argument("--cache-dir", default="cache")
    p.add_argument("--min-run", type=int, default=1_000)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    cols = NAME_COLUMNS + ["locus", "gene_bp"]
    for t in a.tracks:
        for s in SETS:
            cols += [f"{t}_{s}_hit", f"{t}_{s}_longest_run_bp", f"{t}_{s}_gene_bp_in_runs", f"{t}_{s}_runs"]
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for name, contig, start, stop in read_bed(a.bed):
            row = split_name(name) + [f"{contig}:{start}-{stop}", str(stop - start + 1)]
            for t in a.tracks:
                depth, meta = pgs.load_track(a.cache_dir, t, contig)
                lo, hi = max(1, start - FLANK), min(len(depth), stop + FLANK)
                classes = class_index(np.asarray(depth[lo - 1:hi]), int(meta["median"]))
                for s, members in SETS.items():
                    runs = runs_overlapping(np.isin(classes, members), lo, start, stop)
                    for rs, re_ in runs:
                        if (rs == lo and lo > 1) or (re_ == hi and hi < len(depth)):
                            raise RuntimeError(f"{name} {t} {s}: run {rs}-{re_} reaches the {FLANK:,} bp search flank")
                    long_runs = [(rs, re_) for rs, re_ in runs if re_ - rs + 1 >= a.min_run]
                    in_gene = sum(min(re_, stop) - max(rs, start) + 1 for rs, re_ in long_runs)
                    longest = max((re_ - rs + 1 for rs, re_ in long_runs), default=0)
                    row += ["yes" if long_runs else "no", str(longest), str(in_gene),
                            ";".join(f"{rs}-{re_}" for rs, re_ in long_runs)]
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
