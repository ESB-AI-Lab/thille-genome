#!/usr/bin/env python3
"""Thille gene model, CDS exons and plotting window for each top-1% Fst flowering gene.

Genes: --genes (one "NAME|gwen_gene_id" per line). The Thille placement and Tiberius model come from
flowering_beds/igv/flowering_genes_thille.bed (name field 3 = Tiberius gene id, or "unannotated"). CDS exons are the
CDS features of transcript <gene>_t1 in the Tiberius GTF, listed in transcript order (minus-strand genes from the
highest coordinate down). The window is the gene model span (or the placed locus when there is no model) +- --flank bp.
Gwen genes that place on the same Thille model are merged into one row (gwen_ids lists them all).

--extra NAME:MODEL adds a gene that is not in the flowering BED, placed by its Tiberius model alone (e.g. the
stay-green homologs found by blastp of UniProt SGR proteins against the Thille proteome); existing rows are kept.

Usage (from the repo root):
    python3 scripts/gene_models.py --genes summaries/fst_alleles/fst_genes.txt [--extra SGR1:g24740 --extra SGRL:g32057]
Writes summaries/fst_alleles/gene_models.tsv
"""
import argparse
import collections
from pathlib import Path

BED = "flowering_beds/igv/flowering_genes_thille.bed"
GTF = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
OUT = Path("summaries/fst_alleles")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True)
    p.add_argument("--flank", type=int, default=10000)
    p.add_argument("--extra", action="append", default=[], help="NAME:TIBERIUS_GENE_ID (repeatable)")
    a = p.parse_args()
    wanted = [l.strip() for l in open(a.genes) if l.strip()]
    placed = {}
    for l in open(BED):
        if l.startswith("track"):
            continue
        c, s, e, name, _, strand = l.rstrip("\n").split("\t")[:6]
        f = name.split("|")
        placed.setdefault(f"{f[0]}|{f[1]}", (c, int(s) + 1, int(e), f[2], strand))
    missing = [g for g in wanted if g not in placed]
    if missing:
        raise SystemExit(f"not in {BED}: {missing}")
    for x in a.extra:
        name, model = x.split(":")
        placed[f"{name}|{model}"] = (None, None, None, model, None)
        wanted.append(f"{name}|{model}")
    models = {placed[g][3] for g in wanted} - {"unannotated"}
    gtf_loc = {}
    cds = collections.defaultdict(list)
    for l in open(GTF):
        x = l.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] != "CDS":
            continue
        tid = x[8].split('transcript_id "')[1].split('"')[0]
        gid = tid.rsplit("_t", 1)[0]
        if gid in models and tid == f"{gid}_t1":
            cds[gid].append((int(x[3]), int(x[4])))
            gtf_loc[gid] = (x[0], x[6])
    rows = collections.OrderedDict()
    for g in wanted:
        c, s, e, model, strand = placed[g]
        if c is None:  # --extra: chromosome and strand from the GTF
            c, strand = gtf_loc[model]
        name, gwen = g.split("|")
        key = (name, model if model != "unannotated" else f"{c}:{s}-{e}")
        if key in rows:
            rows[key]["gwen_ids"].append(gwen)
            continue
        ex = sorted(cds.get(model, []), reverse=(strand == "-"))
        if model != "unannotated":
            if not ex:
                raise SystemExit(f"{g}: no CDS for {model}_t1 in {GTF}")
            s, e = min(x[0] for x in ex), max(x[1] for x in ex)
        rows[key] = {"gene": name, "gwen_ids": [gwen], "thille_model": model, "chrom": c, "strand": strand,
                     "gene_start": s, "gene_end": e, "win_start": max(1, s - a.flank), "win_end": e + a.flank,
                     "cds_exons": ",".join(f"{x[0]}-{x[1]}" for x in ex) or "none",
                     "cds_bp": sum(x[1] - x[0] + 1 for x in ex)}
    # unique label per figure: gene name, plus the model when one name has several models
    n_listed = len(rows) - len(a.extra)
    by_name = collections.Counter(r["gene"] for r in list(rows.values())[:n_listed])
    a_extra_models = {x.split(":")[1] for x in a.extra}
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["label", "gene", "gwen_ids", "thille_model", "chrom", "strand", "gene_start", "gene_end", "win_start",
            "win_end", "cds_exons", "cds_bp"]
    with open(OUT / "gene_models.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows.values():
            # --extra genes never rename the listed ones: an extra whose name is taken gets the model suffix
            taken = by_name[r["gene"]] > 1 or (r["gwen_ids"][0] in a_extra_models and by_name[r["gene"]] == 1)
            r["label"] = r["gene"] if not taken else f"{r['gene']}_{r['thille_model']}"
            r["gwen_ids"] = ";".join(r["gwen_ids"])
            fh.write("\t".join(str(r[k]) for k in cols) + "\n")
            n_ex = 0 if r["cds_exons"] == "none" else r["cds_exons"].count(",") + 1
            print(f"{r['label']:<22} {r['chrom']}:{r['gene_start']}-{r['gene_end']} {r['strand']} model {r['thille_model']:<12} "
                  f"exons {n_ex:>2} CDS {r['cds_bp']:>5} bp {'(no Thille model: locus only)' if not n_ex else ''}")


if __name__ == "__main__":
    main()
