#!/usr/bin/env python3
"""Figure 1: cumulative assembly length of the Thille chromosomes (v4) and of the hifiasm Hi-C contigs, drawn with
HapSolo2's assembly_stats.py (the tool that drew the original Figure 1; Edwin, 2026-10-10), run inside
/gpool/bin/singularityimages/hapsolo2.sif. Same inputs, genome size and chromosome count as the original
(annotation/repeatmodeler/thille.hic/asmstats/thille_assemblies_assembly_stats.txt: 12 chromosomes, 885,358,400 bp;
2,021 contigs, 989,332,266 bp; expected genome 896 Mb). HapSolo2 always draws a title, so --no-text (Edwin,
2026-10-09: no title inside figures) suppresses Axes.set_title around its plot call; nothing else in its plot changes.
The legend name "Chromosomes," of the original is written without the stray comma.

Usage (from the repo root):
  singularity exec -B /gpool,/epool /gpool/bin/singularityimages/hapsolo2.sif \\
    python3 scripts/plot_cumulative_length.py --no-text -o Results_v4_notext/figure1/thille_cumulative_length_notext
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, "/opt/HapSolo")   # HapSolo2 inside the container
import assembly_stats as A  # noqa: E402
import matplotlib.axes  # noqa: E402

CHR_FASTA = "genome/thille.hifiv4.asm.hic.scaffold.masked.chr.only.fasta"
CTG_FASTA = "/epool/edwin/pamericana/avoflow/hifiasm/thille/hifi.hic.gcluster/thille.hifiv3.asm.hic.p_ctg.fasta"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chr-fasta", default=CHR_FASTA)
    p.add_argument("--ctg-fasta", default=CTG_FASTA)
    p.add_argument("--genome-size", default="896m", help="expected genome size (default: %(default)s)")
    p.add_argument("-n", "--nchroms", type=int, default=12)
    p.add_argument("--title", default="Persea americana Type B Flower: Thille Genome",
                   help="HapSolo2 appends 'Cumulative Length Graph' (default: %(default)s)")
    p.add_argument("--no-text", action="store_true", help="no title (default off)")
    p.add_argument("-o", "--out", default="Results_v4_notext/figure1/thille_cumulative_length",
                   help="output prefix (default: %(default)s); writes .pdf and .png")
    a = p.parse_args()
    gs = A.parse_genome_size(a.genome_size)
    data = [("Chromosomes", A.compute_assembly(a.chr_fasta, genome_size=gs, nchroms=a.nchroms)),
            ("Contigs", A.compute_assembly(a.ctg_fasta, genome_size=gs, nchroms=a.nchroms))]
    for name, d in data:
        print(f"{name}: {len(d['sizes_sorted'])} sequences, {sum(d['sizes_sorted']):,} bp, first {a.nchroms}: "
              f"{d['nxsome']['total']:,} bp")
    assemblies = [(name, d["sizes_sorted"], d["nxsome"]) for name, d in data]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    set_title = matplotlib.axes.Axes.set_title
    if a.no_text:
        matplotlib.axes.Axes.set_title = lambda self, *args, **kw: None
    try:
        for ext in ("png", "pdf"):
            A.plot_cumulative_length(assemblies, a.title, f"{a.out}.{ext}", genome_size=gs, nchroms=a.nchroms)
    finally:
        matplotlib.axes.Axes.set_title = set_title


if __name__ == "__main__":
    main()
