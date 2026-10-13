#!/usr/bin/env python3
"""Count coverage classes genome-wide and per gene from the plot_gene_stacks.py cache.

Classes and median are the ones the figures use: M = int(chromosome median), bins
(0, 0.25M], (0.25M, 0.75M], (0.75M, 1.25M], (1.25M, 1.5M], >1.5M, with depth 0 in "low".
The cache (cache/<track>.<contig>.depth.npy + .json) is filled by plot_gene_stacks.ensure_cache.

Usage (from the repo root):
    python3 scripts/coverage_class_summary.py --tracks gwen.clr.hass_uoq thille.hifi.hass_uoq \\
        --contigs CM056809.1 ... CM056820.1 --bed flowering_beds/igv/flowering_genes_hass.bed \\
        --genome-out summaries/coverage_class_counts_hass.tsv --genes-out summaries/flowering_gene_coverage_classes_hass.tsv
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import plot_gene_stacks as pgs  # noqa: E402

LABELS = pgs.CLASS_LABELS


def class_index(depth, median):
    """Same bins as pgs.classify (pd.cut, right-closed, include_lowest), as integer class indices."""
    edges = np.array(pgs.BIN_FACTORS[1:-1]) * median
    return np.searchsorted(edges, depth, side="left")


NAME_COLUMNS = ["gene", "gwen_gene", "model", "match", "swissprot", "id", "scov"]


def split_name(name):
    """'NRX1|ChrU.ver1.g436710|MRB53_001351|6822of6822|NRX1_ARATH_id50.08_scov1.0' -> the 7 NAME_COLUMNS.

    A fifth field without the _id<x>_scov<y> suffix (e.g. no_named_hit) is kept whole in swissprot.
    """
    f = (name.split("|") + [""] * 5)[:5]
    m = re.fullmatch(r"(.+)_id([0-9.]+)_scov([0-9.]+)", f[4])
    return f[:4] + ([m.group(1), m.group(2), m.group(3)] if m else [f[4], "", ""])


def read_bed(path):
    genes = []
    for line in open(path):
        if line.startswith(("track", "#")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        genes.append((f[3], f[0], int(f[1]) + 1, int(f[2])))  # 1-based inclusive
    return genes


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tracks", nargs="+", required=True)
    p.add_argument("--contigs", nargs="+", required=True)
    p.add_argument("--bed", help="genes to summarize (BED, 0-based)")
    p.add_argument("--coverage-dir", default="covfiles")
    p.add_argument("--cache-dir", default="cache")
    p.add_argument("--jobs", type=int, default=12)
    p.add_argument("--genome-out", required=True)
    p.add_argument("--genes-out")
    a = p.parse_args()

    pgs.ensure_cache({(t, c) for t in a.tracks for c in a.contigs}, a.coverage_dir, a.cache_dir, a.jobs)

    genes = read_bed(a.bed) if a.bed else []
    with open(a.genome_out, "w") as out:
        out.write("alignment\tchromosome\tmedian\tlength\t" + "\t".join(LABELS) + "\n")
        for t in a.tracks:
            totals = np.zeros(len(LABELS), dtype=np.int64)
            for c in a.contigs:
                depth, meta = pgs.load_track(a.cache_dir, t, c)
                median = int(meta["median"])
                counts = np.bincount(class_index(np.asarray(depth), median), minlength=len(LABELS))
                totals += counts
                out.write(f"{t}\t{c}\t{median}\t{len(depth)}\t" + "\t".join(map(str, counts)) + "\n")
            out.write(f"{t}\tgenome\t\t{totals.sum()}\t" + "\t".join(map(str, totals)) + "\n")

    if genes and a.genes_out:
        with open(a.genes_out, "w") as out:
            cols = NAME_COLUMNS + ["locus", "gene_bp"]
            for t in a.tracks:
                cols += [f"{t}_{x}" for x in LABELS] + [f"{t}_median_depth", f"{t}_chrom_median"]
            out.write("\t".join(cols) + "\n")
            for name, contig, start, stop in genes:
                row = split_name(name) + [f"{contig}:{start}-{stop}", str(stop - start + 1)]
                for t in a.tracks:
                    depth, meta = pgs.load_track(a.cache_dir, t, contig)
                    median = int(meta["median"])
                    d = np.asarray(depth[start - 1:stop])
                    frac = np.bincount(class_index(d, median), minlength=len(LABELS)) / len(d)
                    row += [f"{x:.3f}" for x in frac] + [f"{np.median(d):.1f}", str(median)]
                out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
