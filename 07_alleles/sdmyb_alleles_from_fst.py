#!/usr/bin/env python3
"""Rebuild the SDMYB allele files (summaries/sdmyb_alleles/{alleles.fa,codons.tsv,variants.tsv}) from the general
Fst-gene allele pipeline (summaries/fst_alleles/SDMYB/), which supersedes scripts/sdmyb_allele_tables.py.

Why: sdmyb_allele_tables.py found Gwen's allele-2 variants only at sites where the minority base reached 30% of all
reads; allele-2 reads are 29% (37/128), so most of its variants were missed (10 SNPs, a 2-bp deletion, 188 aa).
The per-haplotype calls (scripts/gene_alleles_reads.py) give Gwen allele 2 = Puebla hap2 = Hass alt (190 aa).
The old files are kept in summaries/sdmyb_alleles/_pre_gwen_correction/.

Allele names (SDMYB only): recessive = SDMYB^a (Thille reference CDS), dominant-type = SDMYB^b.

Usage (from the repo root): python3 scripts/sdmyb_alleles_from_fst.py
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fst_common import translate  # noqa: E402
from gene_allele_summary import read_fa  # noqa: E402

SRC = Path("summaries/fst_alleles/SDMYB")
OUT = Path("summaries/sdmyb_alleles")


def main():
    fa = read_fa(SRC / "alleles.fa")
    ref = fa["SDMYB^a_Thille_reference_cds"]
    seqs = {"Thille_reference (SDMYB^a)": ref}
    letters = {}
    for line in open(SRC / "alleles.fa"):
        if line.startswith(">") and "allele=" in line:
            name = line[1:].split()[0]
            letters[name] = line.split("allele=")[1].split()[0]
    for name, s in fa.items():
        if name in letters:
            seqs[f"{name} ({letters[name]})"] = s
    with open(OUT / "alleles.fa", "w") as fh:
        for n, s in seqs.items():
            p = translate(s)
            stop = p.find("*")
            fh.write(f">{n.replace(' ', '_')}_cds len={len(s)}\n{s}\n>{n.replace(' ', '_')}_protein "
                     f"len={stop if stop >= 0 else len(p)}\n{p[:stop + 1] if stop >= 0 else p}\n")
        # the dominant-type CDS is frameshifted past the annotated end; its full protein (to its own stop) comes from
        # the Puebla hap2 assembly, read on past exon 3 (scripts/puebla_sdmyb_alleles.py)
        pb = read_fa("summaries/sdmyb_puebla/puebla_sdmyb_alleles.fa")
        fh.write(f">SDMYB^b_full_protein_from_Puebla_hap2_assembly len={len(pb['Puebla_hap2_protein']) - 1}\n"
                 f"{pb['Puebla_hap2_protein']}\n")
    names = list(seqs)
    n_cod = max(len(s) for s in seqs.values()) // 3
    with open(OUT / "codons.tsv", "w") as fh:
        fh.write("codon_index\t" + "\t".join(f"{n}_codon\t{n}_aa" for n in names) + "\n")
        for c in range(n_cod):
            cells = []
            for n in names:
                cod = seqs[n][c * 3:c * 3 + 3]
                cells += [cod, translate(cod) if len(cod) == 3 else ""]
            fh.write(f"{c + 1}\t" + "\t".join(cells) + "\n")
    shutil.copyfile(SRC / "variants_all.tsv", OUT / "variants.tsv")
    for n in names:
        p = translate(seqs[n]); stop = p.find("*")
        print(f"{n:<34} {len(seqs[n])} nt, {stop if stop >= 0 else 'no stop in CDS'} aa")


if __name__ == "__main__":
    main()
