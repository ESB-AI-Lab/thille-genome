#!/usr/bin/env python3
"""Write the SDMYB allele data (Thille, West Indian and Gwen reads on Thille) as tables for plotting.

Haplotype assignment per read set:
  Thille, WI : two groups from scripts/sdmyb_haplotypes.py (summaries/sdmyb_haplotypes/<label>.groups.tsv)
  Gwen (CLR) : reads carrying >= 60% of the 10 linked CDS alternative alleles they cover (>= 2 covered) are allele 2,
               <= 20% allele 1 (scripts/sdmyb_cds_variants.py; noise masked as in sdmyb_haplotypes.py --clr)
Variants: heterozygous substitution sites from summaries/sdmyb_haplotypes/<label>.tsv and indel clusters from
<label>.indels.tsv (>= 2 bp). For every variant and haplotype: reads carrying it / reads covering it.
Coordinates are Thille chr10 (1-based); SDMYB = g30321 (minus strand), CDS exons from the Tiberius GTF.

Usage (from the repo root): python3 scripts/sdmyb_allele_tables.py
Writes summaries/sdmyb_alleles/{variants.tsv, codons.tsv, alleles.fa}
"""
import collections
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sdmyb_haplotypes as H  # noqa: E402

M = "/epool/edwin/pamericana/avoflow/mappingback"
# THILLE_ASSEMBLY=v4 reads the alignments to the v4 assembly (chr10 is identical in v3 and v4); default v3
SETS = {"v3": {"Thille": (f"{M}/minimap2/thille.hifi.thille.chr10.bam", "Thille_B", "HiFi"),
               "WI": (f"{M}/minimap2/westindian.hifi.thille.chr10.bam", "WI", "HiFi"),
               "Gwen": (f"{M}/gwen.clr.thille.chr10.bam", "Gwen_A_CLR", "CLR")},
        "v4": {"Thille": (f"{M}/v4/thille.hifi.thille_v4.chr10.bam", "Thille_B", "HiFi"),
               "WI": (f"{M}/v4/westindian.hifi.thille_v4.chr10.bam", "WI", "HiFi"),
               "Gwen": (f"{M}/v4/gwen.clr.thille_v4.chr10.bam", "Gwen_A_CLR", "CLR")}}[os.environ.get("THILLE_ASSEMBLY", "v3")]
GWEN_ALT = {10311829: "T", 10311908: "T", 10311912: "T", 10311923: "C", 10311924: "A", 10319296: "T",
            10319339: "A", 10319351: "G", 10319606: "A", 10319646: "C"}
GWEN_CDS_DEL = (10311778, 10311779)  # 2 bp deletion in exon 3 on Gwen allele 2 (CIGAR D after 10,311,777)
LO, HI = 10300000, 10330000
OUT = Path("summaries/sdmyb_alleles")


def region_of(pos):
    names = {0: "exon1", 1: "exon2", 2: "exon3"}
    for i, (s, e) in enumerate(H.CDS):
        if s <= pos <= e:
            return names[i] + "_CDS"
    if H.CDS[2][1] < pos < H.CDS[1][0]:
        return "intron2"
    if H.CDS[1][1] < pos < H.CDS[0][0]:
        return "intron1"
    return "3prime_flank" if pos < H.CDS[2][0] else "5prime_flank"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ref = "".join(H.sh("samtools", "faidx", H.REF, f"chr10:{LO}-{HI}").splitlines()[1:]).upper()
    refb = lambda p: ref[p - LO]
    cds_pos = [p for cs, ce in H.CDS for p in range(ce, cs - 1, -1)]
    refcds = "".join(refb(p).translate(H.COMP) for p in cds_pos)
    tr = lambda nt: "".join(H.CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))

    rows = []
    for name, (bam, label, tech) in SETS.items():
        reads = H.read_events(bam, f"chr10:{LO}-{HI}")
        if name == "Gwen":
            grp = {}
            for n, (s, e, b, _) in reads.items():
                cov = [p for p in GWEN_ALT if p in b and b[p] != "-"]
                if len(cov) >= 2:
                    f = sum(b[p] == GWEN_ALT[p] for p in cov) / len(cov)
                    grp[n] = 2 if f >= 0.6 else 1 if f <= 0.2 else None
            grp = {k: v for k, v in grp.items() if v}
        else:
            grp = {l.split()[0]: int(l.split()[1]) + 1 for l in open(f"summaries/sdmyb_haplotypes/{label}.groups.tsv")}
        sites = [(int(l.split("\t")[2]), l.split("\t")[3], l.split("\t")[4].strip())
                 for l in open(f"summaries/sdmyb_haplotypes/{label}.tsv") if l.startswith("site\tB")]
        for pos, a1, a2 in sites:
            if not LO <= pos <= HI:
                continue
            alt = a2 if a1 == refb(pos) else a1
            rec = {"read_set": name, "technology": tech, "pos": pos, "region": region_of(pos), "type": "SNP",
                   "ref": refb(pos), "alt": alt, "length": 1}
            for h in (1, 2):
                cov = [r for n, r in reads.items() if grp.get(n) == h and r[0] <= pos <= r[1] and r[2].get(pos, "-") != "-"]
                rec[f"hap{h}_with"] = sum(r[2].get(pos) == alt for r in cov)
                rec[f"hap{h}_covering"] = len(cov)
            if pos in cds_pos:
                i = cds_pos.index(pos); c = i // 3; codon = refcds[c * 3:c * 3 + 3]
                new = list(codon); new[i % 3] = alt.translate(H.COMP); new = "".join(new)
                rec.update(codon_index=c + 1, ref_codon=codon, alt_codon=new,
                           aa_change=f"{H.CODE[codon]}{c + 1}{H.CODE.get(new, '?')}")
            rows.append(rec)
        for r in csv.DictReader((l for l in open(f"summaries/sdmyb_haplotypes/{label}.indels.tsv") if not l.startswith("#")),
                                delimiter="\t"):
            pos = int(r["pos"])
            if not LO <= pos <= HI:
                continue
            w1, c1 = map(int, r["group0_with/spanning"].split("/"))
            w2, c2 = map(int, r["group1_with/spanning"].split("/"))
            rows.append({"read_set": name, "technology": tech, "pos": pos, "region": region_of(pos), "type": r["type"],
                         "ref": "", "alt": "", "length": int(r["len"]),
                         "hap1_with": w1, "hap1_covering": c1, "hap2_with": w2, "hap2_covering": c2})
    cols = ["read_set", "technology", "pos", "region", "type", "ref", "alt", "length", "hap1_with", "hap1_covering",
            "hap2_with", "hap2_covering", "codon_index", "ref_codon", "alt_codon", "aa_change"]
    with open(OUT / "variants.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["read_set"], r["pos"])))

    # alleles: Thille 1/2 and WI 1/2 have the reference CDS (0 differences, scripts/sdmyb_haplotypes.py);
    # Gwen allele 2 = 10 SNPs + the exon-3 deletion, translated through to the first stop
    def build(alt=None, dele=()):
        s = {LO + i: b for i, b in enumerate(ref)}
        s.update(alt or {})
        for d in dele:
            s.pop(d, None)
        nt = "".join(s[p] for cs, ce in H.CDS[:-1] for p in range(ce, cs - 1, -1) if p in s)
        nt += "".join(s[p] for p in range(H.CDS[-1][1], H.CDS[-1][0] - 600, -1) if p in s)
        nt = nt.translate(H.COMP)
        prot = ""
        for i in range(0, len(nt) - 2, 3):
            prot += H.CODE.get(nt[i:i + 3], "X")
            if prot[-1] == "*":
                return nt[:i + 3], prot
        return nt, prot
    alleles = {"Thille_ref": build()}
    for n in ("Thille_allele1", "Thille_allele2", "WI_allele1", "WI_allele2", "Gwen_allele1"):
        alleles[n] = alleles["Thille_ref"]
    alleles["Gwen_allele2"] = build(GWEN_ALT, GWEN_CDS_DEL)
    with open(OUT / "alleles.fa", "w") as fh:
        for n, (nt, prot) in alleles.items():
            fh.write(f">{n}_cds len={len(nt)}\n{nt}\n>{n}_protein len={len(prot) - 1}\n{prot}\n")
    with open(OUT / "codons.tsv", "w") as fh:
        names = list(alleles)
        fh.write("codon_index\t" + "\t".join(f"{n}_codon\t{n}_aa" for n in names) + "\n")
        n_codons = max(len(nt) // 3 for nt, _ in alleles.values())
        for c in range(n_codons):
            cells = []
            for n in names:
                nt, prot = alleles[n]
                cells += [nt[c * 3:c * 3 + 3] or "", prot[c] if c < len(prot) else ""]
            fh.write(f"{c + 1}\t" + "\t".join(cells) + "\n")
    print(f"{len(rows)} variant rows; alleles: " + ", ".join(f"{n} {len(p) - 1} aa" for n, (_, p) in alleles.items()))


if __name__ == "__main__":
    main()
