#!/usr/bin/env python3
"""SPA1 5' region in the phenotyped Illumina panel of the Gwen study (Edwin, 2026-10-12).
The Gwen primary contigs hold two versions of the SPA1 5' region: tig00020836 (~9,678,000-9,690,000; the version on
which type B Illumina reads drop in Edwin's dissertation figure) and tig00000203:320,663-339,083 (the version the Thille
SPA1 5' region aligns to uniquely, MAPQ 60). Per cultivar: median depth over the tig00020836 SPA1 5' stretch, over a
flank on the same contig (9,645,000-9,660,000), and over the tig00000203 version; ratio = stretch / flank.
BAMs: gatkrecal BAMs on the Gwen primary contigs (scripts/allele_phylo.py lists; Hass-2 and Hass-R excluded).
Usage (repo root, on SLURM): python3 scripts/spa1_5p_panel_depth.py  -> summaries/v4/spa1_5p_panel_depth.tsv"""
import statistics
import subprocess
import sys
sys.path.insert(0, "scripts")
from allele_phylo import SAMTOOLS, TYPE_A, TYPE_B, bam  # noqa: E402

REG = {"A_version": "tig00020836:9678000-9690000", "flank": "tig00020836:9645000-9660000",
       "B_version": "tig00000203:320663-339083"}


def med(b, r):
    out = subprocess.run([SAMTOOLS, "depth", "-a", "-r", r, str(b)], capture_output=True, text=True, check=True).stdout
    v = [int(l.split("\t")[2]) for l in out.splitlines()]
    return statistics.median(v) if v else float("nan")


with open("summaries/v4/spa1_5p_panel_depth.tsv", "w") as fh:
    fh.write("sample\tflowering_type\tA_version_median\tflank_median\tB_version_median\tA_version_over_flank\n")
    for typ, ss in (("B", TYPE_B), ("A", TYPE_A)):
        for s in ss:
            m = {k: med(bam(s), r) for k, r in REG.items()}
            ratio = m["A_version"] / m["flank"] if m["flank"] else float("nan")
            line = f"{s}\t{typ}\t{m['A_version']}\t{m['flank']}\t{m['B_version']}\t{ratio:.2f}"
            print(line, flush=True); fh.write(line + "\n")
