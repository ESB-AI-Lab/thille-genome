#!/usr/bin/env python3
"""TE overlap of every tandem-array member (Edwin, 2026-10-02: copy sizes, which members are genes rather than TEs,
and array region size).

For each member gene model in summaries/tandem_array_members.tsv: gene span and CDS length (from that table), and the
fraction of the gene span and of the CDS exons covered by the union of RepeatMasker TE intervals
(annotation/...TE.gff). A member is called TE-like when >= 50% of its CDS is TE-covered.
Usage (from the repo root): python3 scripts/tandem_member_te.py
Writes summaries/tandem_array_member_te.tsv
"""
import csv
from collections import defaultdict

TE = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.TE.gff"
GTF = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
mem = list(csv.DictReader(open("summaries/tandem_array_members.tsv"), delimiter="\t"))
chroms = {m["contig"] for m in mem}
te = defaultdict(list)
for l in open(TE):
    if l.startswith("#"):
        continue
    x = l.split("\t")
    if len(x) > 4 and x[0] in chroms:
        te[x[0]].append((int(x[3]), int(x[4])))
for c in te:
    te[c].sort()
cds = defaultdict(list)
tx = {m["transcript"] for m in mem}
for l in open(GTF):
    x = l.rstrip("\n").split("\t")
    if len(x) > 8 and x[2] == "CDS":
        t = x[8].split('transcript_id "')[1].split('"')[0]
        if t in tx:
            cds[t].append((int(x[3]), int(x[4])))


def covered(c, ivs):
    """bp of ivs covered by TE union on contig c"""
    tot = 0
    for s, e in ivs:
        cur = s
        for a, b in te[c]:
            if b < cur:
                continue
            if a > e:
                break
            lo, hi = max(a, cur), min(b, e)
            if lo <= hi:
                tot += hi - lo + 1
                cur = hi + 1
    return tot


out = []
for m in mem:
    c, s, e = m["contig"], int(m["start"]), int(m["end"])
    span = e - s + 1
    cd = cds[m["transcript"]]
    cdbp = sum(b - a + 1 for a, b in cd)
    te_gene = covered(c, [(s, e)]) / span
    te_cds = covered(c, cd) / cdbp if cdbp else float("nan")
    out.append({"array": m["array"], "model": m["model"], "start": s, "end": e, "strand": m["strand"],
                "gene_bp": span, "cds_bp": cdbp, "exons": m["exons"], "exon_identity_best_partner": m["exon_identity"],
                "te_cover_gene": f"{te_gene:.2f}", "te_cover_cds": f"{te_cds:.2f}",
                "call": "TE-like (CDS >= 50% TE)" if te_cds >= 0.5 else "gene"})
with open("summaries/tandem_array_member_te.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
