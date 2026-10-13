#!/usr/bin/env python3
"""Collect every SwissProt / Blast2GO evidence record for a set of Gwen genes (long format).

Genes: the source-flagged union table (ver1 and/or contig genes) plus optional annex genes
(paper coverage-deltas genes outside the union). Evidence sources, each keyed through the verified
gene-ID map (never bare numbers):
  A  2021 A/B contig blastx vs SwissProt (all HSP subjects)          contig gene, 1-based locus
  B  genome-wide contig blastx vs SwissProt (pamer…blastx.outfmt6)  contig gene, 0-based locus
  C  genome-wide contig blastp id mappings (…idmappings.txt)         contig gene, 0-based locus
  D  paper merged scaffold SwissProt blastp (best hit per gene)      scaffold gene_id (GeneID)
  E  paper flowering-flag column (Edwin's manual 'Flowering?' flags) scaffold gene_id
  F  Blast2GO scaffold file (description, GO names, InterPro)        ver1 ID

Usage:
  scripts/gather_gene_evidence.py --union results/final_candidates/avb_top1_union_genes_flags.tsv \
      --out results/final_candidates/evidence_long.tsv --missing-acc results/final_candidates/accessions_missing.txt
"""
import argparse
import collections
import csv
import os
import sys

csv.field_size_limit(sys.maxsize)
P = "publication_data/"
SRC = {
    "A": "gwen_braker_rnaseq_hapsolo_annotation.genes.avo_typeA_v_typeB.top1pct.20k.windowed.weir.spblastx.hits",
    "B": "pamer.contigs.c21.consensus.consensus_pilon_pilon.gene.annotation.genes_swissprot_uniprot.blastx.outfmt6",
    "C": "pamer.contigs.c21.consensus.consensus_pilon_pilon.gene.annotation.genes.swissprot_uniprot.blastp.idmappings.txt",
    "D": P + "gwen.scaffold.gene.annotation.genes.swissprot_uniprot.blastp.sorted.swissprotids.mergedresults.xlsx - gwen scaffold gene annotation g.csv",
    "E": P + "gwen.scaffold.gene.annotation.genes.swissprot_uniprot.blastp.sorted.swissprotids.floweringtype.mergedresults.xlsx - gwen scaffold gene annotation g.csv",
    "F": P + "gwen.scaffolds.20kNs.gene.annotation.genes.b2g.csv",
    "COV": P + "gwen.scaffold.TypeAvsB.top1pct.coverage.deltas.genes.SP.mapped.results.xlsx - gwen.scaffold.TypeAvsB.top1pct.csv",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--union", required=True)
    ap.add_argument("--map-dir", default="reference/gene_id_map")
    ap.add_argument("--annot", default="reference/uniprot/uniprotkb_2026_03_swissprot_annotations.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--missing-acc", required=True, help="write accessions lacking UniProt annotation here")
    a = ap.parse_args()
    for p in (a.out, a.missing_acc):
        if os.path.exists(p):
            sys.exit(f"ERROR: {p} exists; not overwriting")

    ver1 = {r["ver1_id"]: r for r in csv.DictReader(open(os.path.join(a.map_dir, "ver1_genes.tsv")), delimiter="\t")}
    gid2v = {r["scaf_gene_id"]: v for v, r in ver1.items()}
    cg = {r["contig_gene_id"]: r for r in csv.DictReader(open(os.path.join(a.map_dir, "contig_genes.tsv")), delimiter="\t")}
    loc1 = {r["contig_locus_1b"]: j for j, r in cg.items()}
    loc0 = {r["contig_locus_0b"]: j for j, r in cg.items()}
    links_v = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(a.map_dir, "ver1_jg_links.tsv")), delimiter="\t"):
        if float(r["frac_of_jg"]) >= 0.5 or float(r["frac_of_ver1"]) >= 0.5:
            links_v[r["ver1_id"]].append(r["contig_gene_id"])

    genes = []  # (key, set, ver1, [contig genes])
    for r in csv.DictReader(open(a.union), delimiter="\t"):
        js = [j for j in (r["contig_genes_in_top1"].split(";") + r["contig_genes_linked_all"].split(";")) if j]
        genes.append((r["key"], "union", r["ver1_id"], sorted(set(js))))
    in_union = {g[2] for g in genes if g[2]}
    annex = 0
    for r in csv.DictReader(open(SRC["COV"], encoding="utf-8-sig")):
        v = gid2v.get("file_1_file_1_" + r["GeneID"].strip())
        if v and v not in in_union:
            genes.append((v, "annex_coverage", v, sorted(set(links_v.get(v, []))))); annex += 1
    want_j = {j for g in genes for j in g[3]}
    want_v = {g[2] for g in genes if g[2]}
    want_gid = {ver1[v]["scaf_gene_id"].replace("file_1_file_1_", ""): v for v in want_v if ver1[v]["scaf_gene_id"].startswith("file_1_file_1_")}

    ev = collections.defaultdict(list)  # key per gene -> list of dicts
    jkey = collections.defaultdict(list)
    for key, st, v, js in genes:
        for j in js:
            jkey[j].append(key)
    vkey = collections.defaultdict(list)
    for key, st, v, js in genes:
        if v:
            vkey[v].append(key)

    def add(keys, **d):
        for k in keys:
            ev[k].append(d)

    for line in open(SRC["A"]):
        c = line.rstrip("\n").split("\t")
        j = loc1.get(c[0])
        if j in want_j and c[1].startswith("sp|"):
            add(jkey[j], src="A_2021_contig_blastx", gene=j, acc=c[1].split("|")[1], name=c[1].split("|")[2], evalue=c[10], bits=c[11], pident=c[2])
    for line in open(SRC["B"]):
        c = line.rstrip("\n").split("\t")
        j = loc0.get(c[0])
        if j in want_j and c[1].startswith("sp|"):
            add(jkey[j], src="B_genomewide_contig_blastx", gene=j, acc=c[1].split("|")[1], name=c[1].split("|")[2], evalue=c[10], bits=c[11], pident=c[2])
    for line in open(SRC["C"]):
        l0, acc = line.rstrip("\n").split("\t")
        j = loc0.get(l0)
        if j in want_j:
            add(jkey[j], src="C_genomewide_contig_blastp_idmap", gene=j, acc=acc, name="", evalue="", bits="", pident="")
    for r in csv.DictReader(open(SRC["D"], encoding="utf-8-sig")):
        v = want_gid.get(r["GeneID"].strip())
        if v:
            add(vkey[v], src="D_paper_scaffold_blastp", gene=v, acc=r["Entry"], name=r["Entry name"], evalue="", bits="", pident="")
    rows = list(csv.reader(open(SRC["E"], encoding="utf-8-sig")))
    for r in rows[1:]:
        v = want_gid.get(r[-1].strip())
        if v and r[0].strip():
            add(vkey[v], src="E_paper_flowering_flag", gene=v, acc=r[2], name=r[3], evalue="", bits="", pident="", note=f"Flowering?={r[0].strip()}")
    for r in csv.DictReader(open(SRC["F"], encoding="utf-8-sig")):
        if r["GeneID"] in want_v:
            add(vkey[r["GeneID"]], src="F_blast2go", gene=r["GeneID"], acc="", name=r["Description"], evalue=r["e-Value"], bits="", pident=r["sim mean"],
                note=f"GO: {r['GO Names']} | IPR: {r['InterPro IDs']}")

    annot = {r["requested_accession"] for r in csv.DictReader(open(a.annot), delimiter="\t")}
    missing = sorted({d["acc"] for L in ev.values() for d in L if d["acc"] and d["acc"] not in annot})
    cols = ["key", "gene_set", "src", "gene", "acc", "name", "evalue", "bits", "pident", "note"]
    gset = {g[0]: g[1] for g in genes}
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for key, st, v, js in genes:
            for d in ev.get(key, []):
                w.writerow({"key": key, "gene_set": st, "note": "", **d})
    with open(a.missing_acc, "w") as fh:
        for acc in missing:
            fh.write(f"q\tsp|{acc}|X\n")  # pseudo BLAST line for fetch_uniprot.py
    by_src = collections.Counter(d["src"] for L in ev.values() for d in L)
    print(f"genes {len(genes)} (union {len(genes) - annex}, annex {annex}) | genes with any evidence {sum(1 for g in genes if ev.get(g[0]))} | "
          f"records {dict(by_src)} | accessions without UniProt annotation: {len(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()
