#!/usr/bin/env python3
"""Write every contiguous run of one coverage class as a BED line, genome-wide, per alignment track.

Classes are the figures' classes (M = int(chromosome median); plot_gene_stacks.classify), label encoded:
    0 low (d <= 0.25M), 1 half (<= 0.75M), 2 medium (<= 1.25M), 3 high (<= 1.5M), 4 outlier (> 1.5M)
Each output line: contig, start (0-based), end (exclusive), class code, run length (bp), strand "+".
So 1-based positions 500..503 all in class 4 become "contig 499 503 4 4 +".
Depths come from the plot_gene_stacks.py cache (filled from covfiles/ when missing).

Usage (from the repo root):
    python3 scripts/coverage_class_runs.py --track gwen.clr.thille --contigs chr{1..12}
    # writes class_runs/<track>.class_runs.bed.gz
"""

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coverage_class_summary import class_index  # noqa: E402
import plot_gene_stacks as pgs  # noqa: E402


def runs(classes):
    """(start0, end, class) for every maximal run of equal values."""
    change = np.flatnonzero(np.diff(classes)) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [len(classes)]))
    return starts, ends, classes[starts]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--track", required=True)
    p.add_argument("--contigs", nargs="+", required=True)
    p.add_argument("--coverage-dir", default="covfiles")
    p.add_argument("--cache-dir", default="cache")
    p.add_argument("--outdir", default="class_runs")
    p.add_argument("--jobs", type=int, default=12)
    a = p.parse_args()

    pgs.ensure_cache({(a.track, c) for c in a.contigs}, a.coverage_dir, a.cache_dir, a.jobs)
    Path(a.outdir).mkdir(exist_ok=True)
    out_path = Path(a.outdir) / f"{a.track}.class_runs.bed.gz"
    tmp = out_path.with_name(out_path.name + ".tmp")
    total = 0
    with gzip.open(tmp, "wt") as out:
        out.write(f"#track {a.track}; classes 0 low <=0.25M, 1 half <=0.75M, 2 medium <=1.25M, 3 high <=1.5M, 4 outlier >1.5M; "
                  "columns: contig start0 end class run_bp strand\n")
        for contig in a.contigs:
            depth, meta = pgs.load_track(a.cache_dir, a.track, contig)
            median = int(meta["median"])
            classes = class_index(np.asarray(depth), median).astype(np.int8)
            s, e, c = runs(classes)
            if e[-1] != len(depth) or (e - s).sum() != len(depth):
                raise RuntimeError(f"{a.track} {contig}: runs do not tile 1..{len(depth)}")
            out.write(f"#{contig} median {median} length {len(depth)} runs {len(s)}\n")
            out.writelines(f"{contig}\t{x}\t{y}\t{k}\t{y - x}\t+\n" for x, y, k in zip(s.tolist(), e.tolist(), c.tolist()))
            total += len(s)
    tmp.replace(out_path)
    print(f"{out_path}: {total:,} runs", file=sys.stderr)


if __name__ == "__main__":
    main()
