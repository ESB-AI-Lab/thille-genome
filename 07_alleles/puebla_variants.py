#!/usr/bin/env python3
"""Variants of the two Puebla haplotype assemblies against the Thille SDMYB window, in variants.tsv format.

Each Puebla haplotype region (hap1 CM191118.1:40,620,000-40,700,000; hap2 CM191106.1:37,160,000-37,260,000) was
aligned to Thille chr10:10,300,000-10,330,000 with minimap2 2.30 (-c --cs=long -x asm20). SNPs and indels >= 2 bp
come from the cs tags of the primary alignment blocks; between consecutive blocks, unaligned Thille sequence is a
deletion and unaligned Puebla sequence an insertion (both placed at the block junction). An assembly haplotype is
one sequence, so a carried variant is written as 1/1 reads and the other haplotype's columns as 0/0.

Usage (from the repo root): python3 scripts/puebla_variants.py   # writes summaries/sdmyb_alleles/variants_puebla.tsv
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sdmyb_haplotypes as H  # noqa: E402
from sdmyb_allele_tables import region_of  # noqa: E402

S = Path("summaries/sdmyb_puebla")
W0 = 10300000  # window start (1-based position of target offset 0)
COMP = str.maketrans("ACGT", "TGCA")


def variants_from_paf(path):
    blocks, out = [], []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        blocks.append((ts, te, qs, qe))
        t = ts
        for op, val in re.findall(r"([=:*+\-~])([A-Za-z0-9]*)", cs):
            if op == "=":
                t += len(val)
            elif op == ":":
                t += int(val)
            elif op == "*":
                out.append(("SNP", W0 + t, val[0].upper(), val[1].upper(), 1)); t += 1
            elif op == "-":
                if len(val) >= 2:
                    out.append(("DEL", W0 + t - 1, "", "", len(val)))
                t += len(val)
            elif op == "+":
                if len(val) >= 2:
                    out.append(("INS", W0 + t - 1, "", "", len(val)))
    blocks.sort()
    for (t1s, t1e, q1s, q1e), (t2s, t2e, q2s, q2e) in zip(blocks, blocks[1:]):
        if t2s - t1e >= 2:
            out.append(("DEL", W0 + t1e - 1, "", "", t2s - t1e))
        qgap = abs(q2s - q1e) if q2s >= q1e else abs(q1s - q2e)
        if qgap >= 2:
            out.append(("INS", W0 + t1e - 1, "", "", qgap))
    return out, blocks


def main():
    cds_pos = [p for cs, ce in H.CDS for p in range(ce, cs - 1, -1)]
    ref = "".join(H.sh("samtools", "faidx", H.REF, f"chr10:{W0}-{W0 + 30000}").splitlines()[1:]).upper()
    refcds = "".join(ref[p - W0].translate(COMP) for p in cds_pos)
    rows = []
    for h in (1, 2):
        vs, blocks = variants_from_paf(S / f"puebla_hap{h}_vs_thille_window.paf")
        for typ, pos, rb, ab, ln in vs:
            rec = {"read_set": "Puebla", "technology": "assembly (Revio HiFi, hifiasm)", "pos": pos,
                   "region": region_of(pos), "type": typ, "ref": rb, "alt": ab, "length": ln,
                   f"hap{h}_with": 1, f"hap{h}_covering": 1, f"hap{3 - h}_with": 0, f"hap{3 - h}_covering": 0}
            if typ == "SNP" and pos in cds_pos:
                i = cds_pos.index(pos); c = i // 3; codon = refcds[c * 3:c * 3 + 3]
                new = list(codon); new[i % 3] = ab.translate(COMP); new = "".join(new)
                rec.update(codon_index=c + 1, ref_codon=codon, alt_codon=new,
                           aa_change=f"{H.CODE[codon]}{c + 1}{H.CODE.get(new, '?')}")
            rows.append(rec)
        print(f"Puebla hap{h}: {len(blocks)} aligned blocks; {sum(r[0]=='SNP' for r in vs)} SNPs, "
              f"{sum(r[0]=='DEL' for r in vs)} deletions, {sum(r[0]=='INS' for r in vs)} insertions (>= 2 bp)")
    cols = ["read_set", "technology", "pos", "region", "type", "ref", "alt", "length", "hap1_with", "hap1_covering",
            "hap2_with", "hap2_covering", "codon_index", "ref_codon", "alt_codon", "aa_change"]
    with open("summaries/sdmyb_alleles/variants_puebla.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: r["pos"]))


if __name__ == "__main__":
    main()
