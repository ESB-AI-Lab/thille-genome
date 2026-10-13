#!/usr/bin/env python3
"""Coverage classes over the 100 kb figure window (gene + flanks) for every read set on every genome.

For each A/B panel gene and each reference (Thille, Hass UoQ, West Indian), the window is the one the figures use
(plot_gene_stacks.window: WINDOW_SIZE centred on the gene, start rounded to kb). Per track it reports the fraction of
the window in each class, the window median depth / chromosome median, and the class runs >= 1 kb inside the window
(same classes and label encoding as scripts/coverage_class_runs.py: 0 low .. 4 outlier).

Usage (from the repo root):
    python3 scripts/ab_window_classes.py      # writes summaries/ab_window_classes.tsv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coverage_class_summary import class_index, read_bed, split_name  # noqa: E402
import plot_gene_stacks as pgs  # noqa: E402

GENES = ["SDMYB", "SPA1", "AP2", "NAC25", "MIK2", "SUD1", "Y5713", "OTU7", "APT3", "KCS11", "PAR1", "CCR1"]
REFS = {
    "Thille": ("flowering_beds/igv/flowering_genes_thille.bed", ["gwen.clr.thille", "thille.hifi.thille", "westindian.hifi.thille"]),
    "Hass UoQ": ("flowering_beds/igv/flowering_genes_hass.bed", ["gwen.clr.hass_uoq", "thille.hifi.hass_uoq", "westindian.hifi.hass_uoq"]),
    "West Indian": ("flowering_beds/igv/flowering_genes_wi.bed", ["gwen.clr.wi", "thille.hifi.wi", "westindian.hifi.wi"]),
}
READS = {"gwen": "Gwen (A)", "thille": "Thille (B)", "westindian": "West Indian"}
LABELS = pgs.CLASS_LABELS
MIN_RUN = 1_000


def main():
    rows = []
    for ref, (bed, tracks) in REFS.items():
        seen = set()
        for name, contig, start, stop in read_bed(bed):
            fields = split_name(name)
            gene, gwen_gene = fields[0], fields[1]
            if gene not in GENES or (gene == "MIK2" and not gwen_gene.startswith("Chr12")) or (contig, start, stop) in seen:
                continue
            seen.add((contig, start, stop))
            g = pgs.Gene(name, contig, start, stop, "+", gene)
            lo, hi = pgs.window(g, pgs.WINDOW_SIZE)
            for t in tracks:
                if not pgs.cache_paths("cache", t, contig)[1].is_file():
                    continue
                depth, meta = pgs.load_track("cache", t, contig)
                median = int(meta["median"])
                d = np.asarray(depth[max(lo, 1) - 1:hi])
                cls = class_index(d, median)
                frac = np.bincount(cls, minlength=5) / len(cls)
                change = np.flatnonzero(np.diff(cls)) + 1
                s = np.concatenate(([0], change)); e = np.concatenate((change, [len(cls)]))
                long_runs = [f"{lo + a}-{lo + b - 1}:{cls[a]}" for a, b in zip(s, e) if b - a >= MIN_RUN]
                rows.append({"gene": gene, "gwen_gene": gwen_gene, "reference": ref, "reads": READS[t.split(".")[0]],
                             "track": t, "window": f"{contig}:{lo}-{hi}", "gene_span": f"{start}-{stop}",
                             **{f"frac_{c}": round(float(x), 3) for c, x in zip(LABELS, frac)},
                             "window_depth_ratio": round(float(np.median(d)) / median, 3),
                             "runs_ge_1kb": ";".join(long_runs)})
    pd.DataFrame(rows).to_csv("summaries/ab_window_classes.tsv", sep="\t", index=False)
    print(f"{len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
