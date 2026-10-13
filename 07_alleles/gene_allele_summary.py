#!/usr/bin/env python3
"""Merge the per-haplotype calls of each Fst flowering gene and assign superscript allele letters.

Inputs per gene (summaries/fst_alleles/<label>/): <set>.variants.tsv and <set>.fa from scripts/gene_alleles_reads.py
(Thille, WI, Gwen; h1 = group0, h2 = group1; read sets that could not be phased: h1 = consensus, h2 = h1 with
the second allele at exon SNPs, from the hap1_carried / hap2_carried calls) and <asm>.variants.tsv / <asm>.cds.fa from
scripts/gene_alleles_assembly.py (Puebla_hap1, Puebla_hap2, Hass_pri, Hass_alt).

Allele letters (Edwin, 2026-10-01): GENE^a is the Thille reference CDS (the Tiberius model on the Thille assembly);
every other distinct CDS gets the next letter (b, c, ...) in row order Thille h1, Thille h2, WI h1, WI h2, Gwen h1,
Gwen h2, Puebla hap1, Puebla hap2, Hass pri, Hass alt. Two haplotypes share a letter only when their CDS (exons) are identical, so
any coding difference (nonsynonymous or synonymous SNP, insertion, deletion) is a new allele. A CDS with unresolved
bases (N: no reads of that haplotype over part of the exons) gets "?". SDMYB keeps the dominant/recessive names in the
text; its letters are listed for completeness. Genes without a Thille CDS model get no letters.

Allele frequency (AF%, Edwin 2026-10-01): the share of a variety's reads carrying a row's exon allele. Sites: CDS
positions where the variety's two haplotypes differ (a variant carried by one haplotype and not the other, from
variants_all.tsv). At each site only the top two alleles count (the carried variant and the reference allele, which the
phasing groups separate); other bases and gaps (read errors) are ignored. Counts are pooled over the sites, using all
of the variety's reads covering each site (phased or not); indel sites use the indel cluster's reads with / spanning.
When two phased haplotypes do not differ in the exons, the allele is at 100%. A read set that cannot be phased and
whose second allele is < 25% gets the majority allele on both rows (no other cultivar's allele is ever used); its AF is
the real majority-allele frequency at CDS positions with a >= 10% second allele (100% when there is none); when their exon sequences differ only by
changes not listed as variants (1-bp indels), AF is n.d. Assemblies (two haplotypes): 50% or 100%.
Intronic and flanking variation is not used.

Usage (from the repo root): python3 scripts/gene_allele_summary.py
Writes summaries/fst_alleles/<label>/variants_all.tsv, <label>/alleles.fa and summaries/fst_alleles/allele_summary.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fst_common import OUT, CODE, COMP, load_models, ref_window, translate  # noqa: E402
import sdmyb_haplotypes as H  # noqa: E402
from gene_alleles_reads import SETS  # noqa: E402


def carried(r, h):
    if r.get(f"hap{h}_carried", "") != "":  # read sets: the call made by gene_alleles_reads.py
        return str(r[f"hap{h}_carried"]) in ("1", "1.0")
    c = int(r[f"hap{h}_covering"] or 0)
    return c > 0 and int(r[f"hap{h}_with"] or 0) / c >= 0.5


def read_af(g, rs, rows, d):
    """AF of h1 and h2 exon alleles for a read set (see module docstring); None when no reads cover the sites."""
    if not g.exons:
        return None
    vr = [r for r in rows if r["read_set"] == rs and r["region"][0] == "E" and carried(r, 1) != carried(r, 2)]
    if not vr:
        return {1: 1.0, 2: 1.0}
    lo, hi = min(min(e) for e in g.exons), max(max(e) for e in g.exons)
    reads = H.read_events(SETS[rs][0].format(chrom=g.chrom), f"{g.chrom}:{lo}-{hi}")
    ind = {}
    f = d / f"{rs}.indels.tsv"
    for r in csv.DictReader((l for l in open(f) if not l.startswith("#")), delimiter="\t"):
        ind[(int(r["pos"]), r["type"], int(r["len"]))] = (int(r["reads_with"]), int(r["reads_spanning"]))
    n = {1: 0, 2: 0}
    for r in vr:
        h_alt = 1 if carried(r, 1) else 2
        pos = int(r["pos"])
        if r["type"] == "SNP":
            bases = [rd[2].get(pos) for rd in reads.values() if rd[0] <= pos <= rd[1]]
            alt, ref = sum(b == r["alt"] for b in bases), sum(b == r["ref"] for b in bases)
        else:
            w, s = ind.get((pos, r["type"], int(r["length"])), (0, 0))
            alt, ref = w, s - w
        n[h_alt] += alt; n[3 - h_alt] += ref
    tot = n[1] + n[2]
    return {h: n[h] / tot for h in (1, 2)} if tot else None

def majority_af(g, rs, minor=0.10):
    """Unphased read set whose two rows carry the same exon allele (second allele < 25%): the real frequency of that
    majority allele, pooled over CDS positions where a second allele is in >= `minor` of the reads (top two alleles
    only, gaps excluded); 1.0 when no such position exists. Edwin, 2026-10-01: e.g. 90/10 -> both rows get the 90%
    allele, reported as 90%, not 100%."""
    lo, hi = min(min(e) for e in g.exons), max(max(e) for e in g.exons)
    reads = H.read_events(SETS[rs][0].format(chrom=g.chrom), f"{g.chrom}:{lo}-{hi}")
    top = sec = 0
    for pos in g.cds_pos:
        c = {}
        for rd in reads.values():
            if rd[0] <= pos <= rd[1]:
                b = rd[2].get(pos)
                if b and b != "-":
                    c[b] = c.get(b, 0) + 1
        n = sorted(c.values(), reverse=True)
        if len(n) > 1 and n[1] / (n[0] + n[1]) >= minor:
            top += n[0]; sec += n[1]
    return top / (top + sec) if top + sec else 1.0


ROWS = [("Thille", 1), ("Thille", 2), ("WI", 1), ("WI", 2), ("Gwen", 1), ("Gwen", 2), ("Puebla", 1), ("Puebla", 2),
        ("Hass", 1), ("Hass", 2)]
TECH = {"Puebla": "assembly (Revio HiFi, hifiasm)", "Hass": "assembly (UoQ pseudohaplotypes GCA_031208685.1 / GCA_031208655.1)"}
# Hass h1 = principal pseudohaplotype, h2 = alternate. GCA_029852735.1 (the Hass coverage reference) is a mosaic of the
# two at these genes, so it is not used as a row.
ASM = (("Puebla_hap1", "Puebla", 1), ("Puebla_hap2", "Puebla", 2), ("Hass_pri", "Hass", 1), ("Hass_alt", "Hass", 2))


def read_fa(path):
    seqs, k = {}, None
    for l in open(path):
        l = l.strip()
        if l.startswith(">"):
            k = l[1:].split()[0]; seqs[k] = ""
        elif k:
            seqs[k] += l
    return seqs


def aa_diffs(ref_p, p):
    d = [f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(ref_p, p)) if a != b]
    if len(p) != len(ref_p):
        d.append(f"length {len(ref_p)}->{len(p)}")
    return d


def main():
    models = load_models()
    summary = []
    for lab, g in models.items():
        d = OUT / lab
        ref = ref_window(g)
        refb = lambda pos: ref[pos - g.win_start]
        refcds = "".join(g.orient("".join(refb(p) for p in range(s, e + 1))) for s, e in g.exons) if g.exons else ""
        # --- variants: read sets as they are; assemblies converted to with/covering = 1/1 on their haplotype
        rows = []
        for rs in ("Thille", "WI", "Gwen"):
            f = d / f"{rs}.variants.tsv"
            rows += list(csv.DictReader(open(f), delimiter="\t")) if f.exists() else []
        asm = {}
        for name, rs, h in ASM:
            f = d / f"{name}.variants.tsv"
            if not f.exists() or f.read_text().startswith("status"):
                continue
            for r in csv.DictReader(open(f), delimiter="\t"):
                key = (rs, r["type"], int(r["pos"]), int(r["length"]), r["alt"])
                e = asm.setdefault(key, {"read_set": rs, "technology": TECH[rs], "pos": int(r["pos"]),
                                         "region": g.region_of(int(r["pos"])), "type": r["type"], "ref": r["ref"],
                                         "alt": r["alt"], "length": int(r["length"]), "hap1_with": 0, "hap1_covering": 0,
                                         "hap2_with": 0, "hap2_covering": 0})
                e[f"hap{h}_with"] = 1; e[f"hap{h}_covering"] = 1
        rows += list(asm.values())
        for r in rows:  # codon effect of single SNPs inside the CDS (relative to the reference codon)
            pos = int(r["pos"])
            if r["type"] == "SNP" and pos in g.cds_pos:
                i = g.cds_pos.index(pos); c = i // 3; codon = refcds[c * 3:c * 3 + 3]
                new = list(codon); new[i % 3] = g.base_on_gene_strand(r["alt"]); new = "".join(new)
                r.update(codon_index=c + 1, ref_codon=codon, alt_codon=new,
                         aa_change=f"{CODE.get(codon, '?')}{c + 1}{CODE.get(new, '?')}")
        cols = ["read_set", "technology", "pos", "region", "type", "ref", "alt", "length", "hap1_with", "hap1_covering",
                "hap2_with", "hap2_covering", "hap1_carried", "hap2_carried", "unphased", "codon_index", "ref_codon", "alt_codon", "aa_change"]
        with open(d / "variants_all.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="", extrasaction="ignore")
            w.writeheader(); w.writerows(sorted(rows, key=lambda r: (r["read_set"], int(r["pos"]))))
        # --- CDS per haplotype and allele letters
        cds = {}
        for rs in ("Thille", "WI", "Gwen"):
            f = d / f"{rs}.fa"
            if not f.exists() or not g.exons:
                continue
            fa = read_fa(f)
            g0, g1 = fa.get(f"{rs}_group0"), fa.get(f"{rs}_group1")
            cds[(rs, 1)] = g0
            cds[(rs, 2)] = g1 if g1 is not None else g0
            if g1 is None and g0 is not None and len(g0) == len(refcds):
                # not phased: h2 = h1 with the second allele at exon SNPs where the two calls differ
                h2 = list(g0)
                for r in rows:
                    if r["read_set"] == rs and r["type"] == "SNP" and int(r["pos"]) in g.cds_pos and carried(r, 1) != carried(r, 2):
                        i = g.cds_pos.index(int(r["pos"]))
                        h2[i] = g.base_on_gene_strand(r["alt"] if carried(r, 2) else r["ref"])
                cds[(rs, 2)] = "".join(h2)
        for name, rs, h in ASM:
            f = d / f"{name}.cds.fa"
            if f.exists():
                cds[(rs, h)] = read_fa(f)[f"{name}_cds"]
        letters, seen = {}, {}
        if g.exons:
            seen[refcds] = "a"
            nxt = iter("bcdefghijklmnopqrstuvwxyz")
            for row in ROWS:
                s = cds.get(row)
                if s is None:
                    letters[row] = "n/a"
                elif "N" in s:
                    letters[row] = "?"
                else:
                    if s not in seen:
                        seen[s] = next(nxt)
                    letters[row] = seen[s]
        af = {}
        for rs in ("Thille", "WI", "Gwen"):
            x = read_af(g, rs, rows, d) if (d / f"{rs}.variants.tsv").exists() else None
            same = letters.get((rs, 1)) == letters.get((rs, 2)) and letters.get((rs, 1)) not in ("?", "n/a")
            unph = any(r["read_set"] == rs and str(r.get("unphased", "")) in ("1", "1.0") for r in rows)
            if same and unph and g.exons:  # one diploid pool, same allele on both rows: report its real frequency
                m = majority_af(g, rs)
                x = {1: m, 2: m}
            elif x and same:  # two phased haplotypes with the same exon allele
                x = {1: 1.0, 2: 1.0}
            elif x == {1: 1.0, 2: 1.0}:  # exon alleles differ only by changes not in variants_all.tsv (e.g. 1-bp indels)
                x = "n.d."
            for h in (1, 2):
                af[(rs, h)] = x if x == "n.d." else (x[h] if x else None)
        for rs in ("Puebla", "Hass"):
            for h in (1, 2):
                l = letters.get((rs, h))
                af[(rs, h)] = (0.5 + 0.5 * (l == letters.get((rs, 3 - h)))) if l and l not in ("?", "n/a") else None
        ref_p = translate(refcds) if refcds else ""
        with open(d / "alleles.fa", "w") as fh:
            if refcds:
                fh.write(f">{g.gene}^a_Thille_reference_cds len={len(refcds)}\n{refcds}\n>{g.gene}^a_protein\n{ref_p}\n")
            for row in ROWS:
                s = cds.get(row)
                if s is not None:
                    fh.write(f">{row[0]}_h{row[1]} allele={g.gene}^{letters.get(row, '')} len={len(s)}\n{s}\n")
        for row in ROWS:
            s = cds.get(row)
            p = translate(s) if s and "N" not in s else ""
            stop = p.find("*")
            summary.append({"label": lab, "gene": g.gene, "thille_model": g.thille_model, "haplotype": f"{row[0]} h{row[1]}",
                            "allele": f"{g.gene}^{letters[row]}" if row in letters else "no CDS model",
                            "allele_freq_pct": ("n.d." if af.get(row) == "n.d." else f"{100 * af[row]:.0f}")
                            if af.get(row) is not None else "",
                            "cds_len": len(s) if s else "", "frame_ok": (len(s) % 3 == 0) if s else "",
                            "first_stop_aa": (stop + 1 if stop >= 0 else "none") if p else "",
                            "aa_diffs_vs_a": ",".join(aa_diffs(ref_p, p)[:15]) if p else "",
                            "n_aa_diffs": len(aa_diffs(ref_p, p)) if p else "",
                            "nt_diffs_vs_a": sum(x != y for x, y in zip(refcds, s)) + abs(len(refcds) - len(s)) if s and refcds else ""})
        print(f"{lab:<14} " + " ".join(f"{r[0][:3]}{r[1]}:{letters.get(r, '-')}" for r in ROWS))
    cols = list(summary[0])
    with open(OUT / "allele_summary.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(summary)


if __name__ == "__main__":
    main()
