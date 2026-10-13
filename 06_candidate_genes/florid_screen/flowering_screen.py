#!/usr/bin/env python3
"""Screen genes for flowering-associated evidence (tiers T1 FLOR-ID, T2 GO, T3 text keywords).

Input: long evidence table from gather_gene_evidence.py plus UniProt annotation tables.
Rules (Edwin, 2026-09-24):
  T1  any hit whose Arabidopsis locus (Araport AGI) is one of the FLOR-ID flowering-time genes
  T2  flowering GO biological-process terms from ANY hit's UniProt entry, or from Blast2GO GO names
  (light signalling, i.e. photomorphogenesis/COP1/blue-red-far-red light, counts since 2026-09-24)
  T3  flowering keywords in protein names / UniProt function text / Blast2GO description; T3 hits are
      a pre-filter for LLM-assisted review (verdicts are added separately, never by this script)
  TE-derived hits (transposon polyproteins etc.) are excluded as evidence, so TE-flagged genes are
  screened on their non-TE hits. BLAST hits need e-value <= --max-evalue when an e-value is given.

Usage:
  scripts/flowering_screen.py --evidence results/final_candidates/evidence_long.tsv \
      --annot reference/uniprot/uniprotkb_2026_03_swissprot_annotations.tsv \
      --annot reference/uniprot_supp_20260924/uniprotkb_2026_03_swissprot_annotations.tsv \
      --florid reference/florid/florid_flowering_genes.tsv --out results/final_candidates/screen.tsv
"""
import argparse
import collections
import csv
import os
import re
import sys

csv.field_size_limit(sys.maxsize)
TE = re.compile(r"retrovirus-related pol|retrotransposon|transposon|transposase|reverse transcriptase|gag-pol|copia|"
                r"line-1|gypsy|polyprotein", re.I)
GO = re.compile(r"flower|photoperiod|circadian|rhythm|vernaliz|floral|meristem identity|short-day|long-day|gibberellin|"
                r"photomorphogenesis|blue light|red light|far-red|red, far-red", re.I)  # light signalling added 2026-09-24 (Edwin)
KW = re.compile(r"flower|floral|photoperiod|circadian|clock|florigen|vernaliz|heading date|day length|"
                r"CONSTANS|FLOWERING LOCUS|TERMINAL FLOWER|LEAFY|APETALA|AGAMOUS|SEPALLATA|PISTILLATA|MADS|"
                r"SUPPRESSOR OF OVEREXPRESSION OF CO|SHORT VEGETATIVE PHASE|SQUAMOSA|GIGANTEA|EARLY FLOWERING|"
                r"LATE ELONGATED HYPOCOTYL|CIRCADIAN CLOCK ASSOCIATED|TIMING OF CAB|PSEUDO-RESPONSE REGULATOR|"
                r"PHYTOCHROME|CRYPTOCHROME|SUPPRESSOR OF PHYA|ELONGATED HYPOCOTYL|CYCLING DOF|FLAVIN-BINDING, KELCH|"
                r"ZEITLUPE|TIME FOR COFFEE|gibberellin|DELLA|pollen|anther|stamen|carpel|pistil|stigma|tapetum|"
                r"dichogamy|benzenoid|floral scent|EMISSION OF BENZENOIDS|MYB21|MYB24|MYB57|MYB305|MYB340|"
                r"photomorphogen|COP1|CONSTITUTIVE PHOTOMORPHOGENIC|blue light|red light|far-red|phototropin|photoreceptor", re.I)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--annot", action="append", required=True)
    ap.add_argument("--florid", required=True)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"ERROR: {a.out} exists; not overwriting")
    ann = {}
    for p in a.annot:
        for r in csv.DictReader(open(p), delimiter="\t"):
            ann.setdefault(r["requested_accession"], r)
    flor = {r["AGI"]: r for r in csv.DictReader(open(a.florid), delimiter="\t")}

    per = collections.defaultdict(lambda: {"set": "", "T1": [], "T2": [], "T3": [], "te_hits": 0, "hits": 0, "flag": []})
    for r in csv.DictReader(open(a.evidence), delimiter="\t"):
        g = per[r["key"]]; g["set"] = r["gene_set"]
        if r["src"] == "E_paper_flowering_flag":
            g["flag"].append(r["note"]); continue
        if r["src"] == "F_blast2go":
            txt = r["name"] + " " + r["note"]
            if TE.search(r["name"]):
                continue
            go = [t.strip() for t in r["note"].split("|")[0].replace("GO:", "", 1).split(";") if t.strip().startswith("P:") and GO.search(t)]
            if go:
                g["T2"].append(f"B2G: {'; '.join(go)}")
            if KW.search(r["name"]):
                g["T3"].append(f"B2G desc: {r['name']}")
            continue
        if r["evalue"]:
            try:
                if float(r["evalue"]) > a.max_evalue:
                    continue
            except ValueError:
                pass
        x = ann.get(r["acc"])
        if not x:
            continue
        g["hits"] += 1
        pname, func = x.get("Protein names", ""), x.get("Function [CC]", "")
        if TE.search(pname):
            g["te_hits"] += 1; continue
        tag = f"{x['Entry Name']} ({r['src'][0]})"
        agis = [s.strip() for s in x.get("Araport", "").split(";") if s.strip()]
        for agi in agis:
            if agi in flor:
                g["T1"].append(f"{tag} {agi} FLOR-ID:{flor[agi]['short_names']} [{flor[agi]['keyword'][:40]}]")
        go = [t.split(" [")[0].strip() for t in x.get("Gene Ontology (biological process)", "").split(";") if GO.search(t)]
        if go:
            g["T2"].append(f"{tag} GO: {'; '.join(sorted(set(go)))}")
        m = KW.search(pname) or KW.search(func)
        if m:
            snip = pname if KW.search(pname) else func[max(0, m.start() - 60): m.end() + 60]
            g["T3"].append(f"{tag}: {snip}")
    cols = ["key", "gene_set", "tiers", "n_hits_used", "n_te_hits_excluded", "paper_flowering_flag", "T1_florid", "T2_go", "T3_keyword"]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
        for k, g in per.items():
            tiers = ",".join(t for t in ("T1", "T2", "T3") if g[t])
            w.writerow([k, g["set"], tiers, g["hits"], g["te_hits"], ";".join(g["flag"]),
                        " || ".join(sorted(set(g["T1"]))), " || ".join(sorted(set(g["T2"]))), " || ".join(sorted(set(g["T3"])))])
    c = collections.Counter()
    for k, g in per.items():
        for t in ("T1", "T2", "T3"):
            if g[t]:
                c[(g["set"], t)] += 1
        if g["T1"] or g["T2"] or g["T3"]:
            c[(g["set"], "any")] += 1
    print("genes with evidence:", len(per), "|", dict(sorted(c.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
