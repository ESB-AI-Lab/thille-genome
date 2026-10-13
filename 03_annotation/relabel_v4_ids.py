#!/usr/bin/env python3
"""Write a v4 copy of a text file with every Thille Tiberius gene ID replaced by its v4 ID (Edwin, 2026-10-04).
IDs are matched as standalone tokens 'g<digits>' (also inside 'g123_t1', 'NAME_g123', 'thille_g123'); Gwen IDs such as
Chr3.ver1.g154150, jg35492 and m1-g59046 are not touched (preceded by '.', 'j' or '-'). Uses
summaries/gene_id_map_v3_v4.tsv (status 'identical' only); an ID without a one-to-one v4 match stops the script.
Gwen ver1 IDs written without their prefix (e.g. 'g414860') are left as they are; v3 chr8 genes re-predicted in v4 are
mapped explicitly (EXTRA), each to the v4 gene over its lifted span.
Usage: python3 scripts/relabel_v4_ids.py <in> <out>   (out must not exist)"""
import re, sys, os
src, dst = sys.argv[1], sys.argv[2]
SKIP = {"g414860", "g440620", "g440630", "g94530", "g154550", "g477040"}   # Gwen ver1 IDs written without 'ChrN.ver1.'
EXTRA = {"g33824": "g23529"}   # FTIP3 annex gene (not promoted, not in the manuscript): v3 chr8 model re-predicted in v4
if os.path.exists(dst): sys.exit(f"exists: {dst}")
m = {}
for l in open("summaries/gene_id_map_v3_v4.tsv"):
    f = l.rstrip("\n").split("\t")
    if f[0] != "v3_gene" and f[6] == "identical": m[f[0]] = f[1]
PAT = re.compile(r"(?<![A-Za-z0-9.\-])g(\d+)(?![0-9])")
missing = set(); n = 0
def sub(mt):
    global n
    k = "g" + mt.group(1)
    if k in SKIP: return k
    if k in EXTRA: n += 1; return EXTRA[k]
    if k not in m: missing.add(k); return k
    n += k != m[k]; return m[k]
out = PAT.sub(sub, open(src, newline="").read())   # keep the original line endings
if missing: sys.exit(f"{src}: no identical v4 model for {sorted(missing)[:20]}")
open(dst, "w", newline="").write(out); print(f"{dst}: {n} IDs changed")
