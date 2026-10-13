"""Table 2 (genome annotation) source numbers: gene-model statistics from the Tiberius GTF and protein FASTA,
total masked from the RepeatMasker .tbl of the 12 chromosomes, repeat classes as
grouped in Figure 2 (thille_avocado_composition.tsv), and BUSCO from the two short_summary JSON files.

Usage: python3 scripts/annotation_table.py [--assembly v3|v4] > summaries/tables/annotation[_v4].tsv
"""
import json, re, statistics, sys
from pathlib import Path

A = Path("/epool/edwin/pamericana/avoflow/annotation")
import argparse
_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--assembly", choices=["v3", "v4"], default="v3", help="Thille assembly version (default: v3, the original Table 2)")
ASM = _ap.parse_args().assembly
if ASM == "v3":
    GTF = A / "tiberius/thille.hic/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
    PROT = A / "tiberius/thille.hic/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.proteins.fasta"
    TBL = A / "repeatmodeler/thille.hic/thille.top12.fasta.tbl"
    COMP = A / "repeatmodeler/thille.hic/thille_avocado_composition.tsv"  # Figure 2 data
    BUSCO = {"Predicted proteins": A / "tiberius/thille.hic/busco/thille_hic_proteins/short_summary.specific.embryophyta_odb12.2.thille_hic_proteins.json",
             "12-chromosome assembly": A / "tiberius/thille.hic/busco/thille_hic_genome/short_summary.specific.embryophyta_odb12.2.thille_hic_genome.json"}
else:  # v4 (complete chr8, 2026-10-04)
    _v = "thille.hifiv4.asm.hic.scaffold.masked.chr.only"
    GTF = A / f"tiberius/thille.hic/{_v}.annotation.tiberius.gtf"
    PROT = A / f"tiberius/thille.hic/{_v}.annotation.tiberius.proteins.fasta"
    TBL = A / "repeatmodeler/thille.hic/repeatmask/thille.hifiv4.asm.hic.scaffold.chr.only.fasta.tbl"
    COMP = A / "repeatmodeler/thille.hic/thille_genome_composition.tsv"  # v4 composition (Figure 2 data)
    _b = A / "tiberius/thille.hic/busco"
    BUSCO = {"Predicted proteins": _b / f"{_v}.annotation.tiberius.proteins_embryophyta_odb12.2/short_summary.specific.embryophyta_odb12.2.{_v}.annotation.tiberius.proteins_embryophyta_odb12.2.json",
             "12-chromosome assembly": _b / f"{_v}_embryophyta_odb12.2/short_summary.specific.embryophyta_odb12.2.{_v}_embryophyta_odb12.2.json"}


def union_bp(intervals):
    total, end = 0, 0
    for s, e in sorted(intervals):
        if e > end:
            total += e - max(s - 1, end)
            end = e
    return total


genes, exons_by_tx, ex_iv = {}, {}, {}
for line in open(GTF):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    s, e = int(f[3]), int(f[4])
    if f[2] == "gene":
        genes[re.search(r'gene_id "([^"]+)"|^(\S+)$', f[8]).group(1) or f[8]] = (f[0], s, e)
    elif f[2] == "exon":
        tx = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons_by_tx[tx] = exons_by_tx.get(tx, 0) + 1
        ex_iv.setdefault(f[0], []).append((s, e))
glen = [e - s + 1 for _c, s, e in genes.values()]
gene_iv = {}
for c, s, e in genes.values():
    gene_iv.setdefault(c, []).append((s, e))
asm = int(re.search(r"total length:\s+(\d+) bp", TBL.read_text()).group(1))
span = sum(union_bp(v) for v in gene_iv.values())
exon_bp = sum(union_bp(v) for v in ex_iv.values())
nex = list(exons_by_tx.values())

plen, seq = [], 0
for line in open(PROT):
    if line.startswith(">"):
        if seq:
            plen.append(seq)
        seq = 0
    else:
        seq += len(line.strip().rstrip("*"))
plen.append(seq)

rows = [("Genes", "Protein-coding genes (one transcript each)", f"{len(genes):,}"),
        ("Genes", "Gene length, mean / median (bp)", f"{statistics.mean(glen):,.0f} / {statistics.median(glen):,.0f}"),
        ("Genes", "Gene length, range", f"{min(glen):,} bp to {max(glen) / 1000:,.1f} kb"),
        ("Genes", "Assembly covered by genes / by exons (%)", f"{100 * span / asm:.2f} / {100 * exon_bp / asm:.2f}"),
        ("Genes", "Exons per gene, mean / median", f"{statistics.mean(nex):.1f} / {statistics.median(nex):.0f}"),
        ("Genes", "Single-exon genes", f"{sum(n == 1 for n in nex):,} ({100 * sum(n == 1 for n in nex) / len(nex):.1f}%)"),
        ("Genes", "Protein length, mean / median (aa)", f"{statistics.mean(plen):.0f} / {statistics.median(plen):.0f}"),
        ("Genes", "Protein length, range (aa)", f"{min(plen):,} to {max(plen):,}")]

tbl = TBL.read_text()
bp, pct = re.search(r"bases masked:\s+(\d+) bp \(\s*([\d.]+) %\)", tbl).groups()
rows.append(("Repeats", "Total masked", f"{int(bp) / 1e6:.1f} Mb ({float(pct):.2f}%)"))
# classes as grouped in Figure 2 (same RepeatMasker run; percent of the assembly excluding N)
for line in list(open(COMP))[1:]:
    cls, sub, cbp, cpct = line.rstrip("\n").split("\t")
    label = f"  {sub}" if sub else cls
    rows.append(("Repeats", label, f"{int(cbp) / 1e6:.1f} Mb ({float(cpct):.2f}%)"))

for label, path in BUSCO.items():
    r = json.load(open(path))["results"]
    rows.append(("BUSCO (embryophyta_odb12.2, n = 1,658)", label,
                 f"C {r['Complete BUSCOs']:,} ({r['Complete percentage']}%) [S {r['Single copy BUSCOs']:,}, "
                 f"D {r['Multi copy BUSCOs']}], F {r['Fragmented BUSCOs']}, M {r['Missing BUSCOs']}"))

print("section\tmetric\tvalue")
for row in rows:
    print("\t".join(row))
