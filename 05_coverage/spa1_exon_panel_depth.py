#!/usr/bin/env python3
"""SPA1 exons 1 and 2 on Gwen contig tig00020836 (9,680,000-9,685,000; Edwin's dissertation figure) in the phenotyped
Illumina panel of the Gwen study (Edwin, 2026-10-12). Exons lifted from Thille g6249: exon 1 9,681,873-9,683,384,
exon 2 9,680,819-9,680,976. Per cultivar: median depth over the window, exon 1 and exon 2, the flank on the same contig
(9,645,000-9,660,000), the ratios to the flank, and the share of window bases with at most 2 reads.
Usage (repo root, on SLURM): python3 scripts/spa1_exon_panel_depth.py  -> summaries/v4/spa1_exon_panel_depth.tsv"""
import statistics
import subprocess
import sys
sys.path.insert(0, "scripts")
from allele_phylo import SAMTOOLS, TYPE_A, TYPE_B, bam  # noqa: E402

REG = {"window": "tig00020836:9680000-9685000", "exon1": "tig00020836:9681873-9683384",
       "exon2": "tig00020836:9680819-9680976", "flank": "tig00020836:9645000-9660000"}


def depths(b, r):
    out = subprocess.run([SAMTOOLS, "depth", "-a", "-r", r, str(b)], capture_output=True, text=True, check=True).stdout
    return [int(l.split("\t")[2]) for l in out.splitlines()]


with open("summaries/v4/spa1_exon_panel_depth.tsv", "w") as fh:
    fh.write("sample\ttype\twindow\texon1\texon2\tflank\twindow_ratio\texon1_ratio\texon2_ratio\twindow_bases_le2\n")
    for typ, ss in (("B", TYPE_B), ("A", TYPE_A)):
        for s in ss:
            d = {k: depths(bam(s), r) for k, r in REG.items()}
            m = {k: statistics.median(v) for k, v in d.items()}
            le2 = sum(x <= 2 for x in d["window"]) / len(d["window"])
            f = m["flank"] or float("nan")
            line = (f"{s}\t{typ}\t{m['window']}\t{m['exon1']}\t{m['exon2']}\t{m['flank']}\t{m['window'] / f:.2f}\t"
                    f"{m['exon1'] / f:.2f}\t{m['exon2'] / f:.2f}\t{le2:.2f}")
            print(line, flush=True); fh.write(line + "\n")
