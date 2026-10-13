#!/usr/bin/env python3
"""Flowering-type allele classes under the SDMYB model, scored without any reference default.

Model (Groh et al. 2026 for SDMYB): B-type individuals (Thille, Lowland/WI) are homozygous for one allele, and each
A-type individual (Gwen, Hass MH2022, Puebla) carries that allele on one haplotype and a different allele on the other.

Sites: every variant in summaries/fst_alleles/<label>/variants_all.tsv (key: position, type, length, SNP alternative).
State of a haplotype at a site: carried (the haplotype's call), not carried (no call but the haplotype has data at the
position: >= 3 reads for read haplotypes, an aligned block for assemblies), or missing. A site is scored only when all
10 haplotypes have a state; it follows the pattern when the four B haplotypes share one state and each A individual
has exactly one haplotype in the other state. Missing data are never replaced by the Thille reference (Edwin,
2026-10-02); sites with any missing haplotype are counted as not scored.

Alleles: the CDS letters of summaries/fst_alleles/allele_summary.tsv. A gene matches the model at the allele level when
the four B haplotypes carry one allele, each A individual carries that allele on exactly one haplotype, and the A
individuals' other allele is the same allele in all three; tested on CDS (letters) and on proteins (translated CDS,
so synonymous differences are ignored). Haplotypes with '?' or 'n/a' make the gene not scorable.

Usage (from the repo root): python3 scripts/allele_class_pattern.py
Writes summaries/fst_alleles/allele_class_pattern.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fst_common import OUT, load_models, translate  # noqa: E402
from gene_allele_summary import read_fa  # noqa: E402

B = [("Thille", 1), ("Thille", 2), ("WI", 1), ("WI", 2)]
A = {"Gwen": [("Gwen", 1), ("Gwen", 2)], "Hass": [("Hass", 1), ("Hass", 2)], "Puebla": [("Puebla", 1), ("Puebla", 2)]}
ROWS = B + [r for v in A.values() for r in v]
ASM = {("Puebla", 1): "Puebla_hap1", ("Puebla", 2): "Puebla_hap2", ("Hass", 1): "Hass_pri", ("Hass", 2): "Hass_alt"}


def coverage(g, d):
    cov = defaultdict(list)
    for rs in ("Thille", "WI", "Gwen"):
        f = d / f"{rs}.covered.tsv"
        if f.exists():
            for r in csv.DictReader(open(f), delimiter="\t"):
                cov[(rs, int(r["hap"]))].append((int(r["start"]), int(r["end"])))
    for row, name in ASM.items():
        f = d / f"{name}.paf"
        if f.exists():
            for line in open(f):
                x = line.split("\t")
                if "tp:A:P" in x[12:]:
                    cov[row].append((g.win_start + int(x[7]), g.win_start + int(x[8]) - 1))
    return cov


def flag(r, h):
    v = r.get(f"hap{h}_carried", "")
    if v != "":
        return str(v) in ("1", "1.0")
    c = int(r[f"hap{h}_covering"] or 0)
    return c > 0 and int(r[f"hap{h}_with"] or 0) / c >= 0.5


def main():
    models = load_models()
    letters = defaultdict(dict)
    for r in csv.DictReader(open(OUT / "allele_summary.tsv"), delimiter="\t"):
        cv, h = r["haplotype"].split(" h")
        letters[r["label"]][(cv, int(h))] = r["allele"].split("^")[1] if "^" in r["allele"] else ""
    out = []
    for lab, g in models.items():
        d = OUT / lab
        f = d / "variants_all.tsv"
        if not f.exists():
            continue
        cov = coverage(g, d)
        covered = lambda row, pos: any(s <= pos <= e for s, e in cov.get(row, []))
        carriers = defaultdict(set)
        region = {}
        for r in csv.DictReader(open(f), delimiter="\t"):
            key = (int(r["pos"]), r["type"], int(r["length"]), r["alt"] if r["type"] == "SNP" else "")
            region[key] = r["region"]
            for h in (1, 2):
                if flag(r, h) and (r["read_set"], h) in ROWS:
                    carriers[key].add((r["read_set"], h))
        n_sites = n_scored = n_pat = n_pat_cds = 0
        for key, car in carriers.items():
            if not car:
                continue
            n_sites += 1
            state = {}
            for row in ROWS:
                state[row] = True if row in car else (False if covered(row, key[0]) else None)
            if any(v is None for v in state.values()):
                continue
            n_scored += 1
            bs = {state[r] for r in B}
            if len(bs) == 1:
                s = bs.pop()
                if all(sum(state[r] != s for r in rows) == 1 for rows in A.values()):
                    n_pat += 1
                    n_pat_cds += region[key][0] == "E"
        L = letters.get(lab, {})
        res = {"cds": "not scorable", "protein": "not scorable"}
        if g.exons and L and all(L.get(r) not in (None, "", "?", "n/a") for r in ROWS):
            fa = read_fa(d / "alleles.fa")
            seq = {}
            for name, s in fa.items():
                if "_h" in name and not name.startswith(g.gene + "^"):
                    cv, h = name.rsplit("_h", 1)
                    seq[(cv, int(h))] = s
            prot = {r: translate(seq[r]).split("*")[0] for r in ROWS if r in seq}
            for kind, val in (("cds", L), ("protein", prot)):
                if len(val) < len(ROWS):
                    continue
                bset = {val[r] for r in B}
                if len(bset) != 1:
                    res[kind] = "no (B haplotypes differ)"
                    continue
                b = bset.pop()
                others = []
                ok = True
                for rows in A.values():
                    diff = [val[r] for r in rows if val[r] != b]
                    if len(diff) != 1:
                        ok = False
                        break
                    others.append(diff[0])
                res[kind] = ("yes" if len(set(others)) == 1 else "partial (A alleles differ)") if ok else "no"
        out.append({"label": lab, "gene": g.gene, "variant_sites": n_sites, "sites_scored": n_scored,
                    "sites_missing_data": n_sites - n_scored, "pattern_sites": n_pat, "pattern_sites_cds": n_pat_cds,
                    "allele_match_cds": res["cds"], "allele_match_protein": res["protein"]})
    cols = list(out[0])
    with open(OUT / "allele_class_pattern.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(out)
    for r in sorted(out, key=lambda r: -r["pattern_sites"]):
        print(f"{r['label']:<18} sites {r['variant_sites']:>5} scored {r['sites_scored']:>5} pattern {r['pattern_sites']:>4} "
              f"(CDS {r['pattern_sites_cds']:>2})  CDS {r['allele_match_cds']:<26} protein {r['allele_match_protein']}")


if __name__ == "__main__":
    main()
