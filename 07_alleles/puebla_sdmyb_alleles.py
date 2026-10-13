#!/usr/bin/env python3
"""SDMYB coding sequence and protein of the two Puebla haplotype assemblies.

The three Thille g30321 CDS exons (summaries/sdmyb_puebla/thille_exons.fa, Thille plus strand) are found in each
Puebla region (summaries/sdmyb_puebla/puebla_hap{1,2}_region.fa) with blastn (best HSP per exon). Exons are joined in
transcript order on the gene strand, and the sequence after exon 3 is read on until the first stop codon, so a
frameshifted allele is translated to its own stop.

Usage (from the repo root): python3 scripts/puebla_sdmyb_alleles.py
Writes summaries/sdmyb_puebla/puebla_sdmyb_alleles.fa (CDS through the stop codon, and protein, per haplotype)
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sdmyb_haplotypes as H  # noqa: E402

S = Path("summaries/sdmyb_puebla")
rc = lambda x: x.translate(H.COMP)[::-1]


def fasta_seq(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def main():
    out = []
    for h in ("hap1", "hap2"):
        reg = fasta_seq(S / f"puebla_{h}_region.fa")
        res = subprocess.run(["singularity", "exec", H.IMG, "blastn", "-query", str(S / "thille_exons.fa"),
                              "-subject", str(S / f"puebla_{h}_region.fa"), "-outfmt", "6 qseqid sstart send",
                              "-evalue", "1e-5", "-max_hsps", "1"], capture_output=True, text=True, check=True).stdout
        best = {}
        for l in res.splitlines():
            q, ss, se = l.split("\t")
            best.setdefault(q, (int(ss), int(se)))
        # ss > se: the Puebla plus strand runs opposite to Thille's plus strand, i.e. already in the CDS direction
        cds = ""
        for q in ("exon1", "exon2", "exon3"):
            ss, se = best[q]
            seg = reg[min(ss, se) - 1:max(ss, se)]
            cds += seg if ss > se else rc(seg)
        ss, se = best["exon3"]
        tail = reg[max(ss, se):max(ss, se) + 300] if ss > se else rc(reg[min(ss, se) - 301:min(ss, se) - 1])
        full, prot = cds + tail, ""
        for i in range(0, len(full) - 2, 3):
            prot += H.CODE.get(full[i:i + 3], "X")
            if prot[-1] == "*":
                break
        nt = full[:len(prot) * 3]
        out.append(f">Puebla_{h}_cds len={len(nt)} exon_nt={len(cds)}\n{nt}\n>Puebla_{h}_protein len={len(prot) - 1}\n{prot}\n")
        print(f"Puebla {h}: exon nt {len(cds)}, CDS to stop {len(nt)}, protein {len(prot) - 1} aa")
    (S / "puebla_sdmyb_alleles.fa").write_text("".join(out))


if __name__ == "__main__":
    main()
