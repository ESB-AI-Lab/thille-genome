#!/usr/bin/env python3
"""Build cooler chromsizes files for the Thille Hi-C map from the YaHS final AGP and the
scaffold->chromosome mapping.

Scaffold length = largest object_end in the AGP. Outputs (in <outdir>):
  chr12.sizes  mapped scaffolds only, renamed to chromosome names, ordered by chromosome number
  all.sizes    the same chromosomes first, then every other scaffold by length (desc), name kept
  rename.tsv   scaffold -> name used in the maps (identity for unmapped scaffolds)

  <name>.sizes for each --extra-scope NAME:scafA,scafB  (the mapped chromosomes, then the listed scaffolds)

Usage: make_chromsizes.py <scaffolds_final.agp> <scaffold_to_chrom.tsv> <outdir> [--extra-scope chr12scf13:scaffold_13]
"""
import csv
import os
import re
import sys

args = sys.argv[1:]
extra = []
while "--extra-scope" in args:
    i = args.index("--extra-scope")
    if i + 1 >= len(args) or ":" not in args[i + 1]:
        sys.exit(__doc__)
    name, scafs = args[i + 1].split(":", 1)
    extra.append((name, [x for x in scafs.split(",") if x]))
    del args[i:i + 2]
if len(args) != 3:
    sys.exit(__doc__)
agp, mapping, outdir = args

lengths = {}
with open(agp) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        lengths[f[0]] = max(lengths.get(f[0], 0), int(f[2]))

rename = {}
with open(mapping) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        scaf, chrom = r["scaffold"], r["chr_name"]
        if scaf not in lengths:
            sys.exit(f"ERROR: {scaf} from {mapping} not in {agp}")
        if int(r["scaffold_len"]) != lengths[scaf]:
            sys.exit(f"ERROR: {scaf} length {r['scaffold_len']} in mapping != {lengths[scaf]} in AGP")
        if chrom in rename.values():
            sys.exit(f"ERROR: {chrom} assigned to more than one scaffold")
        rename[scaf] = chrom


def chrom_key(name):
    m = re.match(r"^chr(\d+)$", name)
    if not m:
        sys.exit(f"ERROR: cannot order chromosome name '{name}'")
    return int(m.group(1))


chroms = sorted(rename.items(), key=lambda kv: chrom_key(kv[1]))
others = sorted((s for s in lengths if s not in rename), key=lambda s: (-lengths[s], s))

os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "chr12.sizes"), "w") as fh:
    for scaf, chrom in chroms:
        fh.write(f"{chrom}\t{lengths[scaf]}\n")
with open(os.path.join(outdir, "all.sizes"), "w") as fh:
    for scaf, chrom in chroms:
        fh.write(f"{chrom}\t{lengths[scaf]}\n")
    for scaf in others:
        fh.write(f"{scaf}\t{lengths[scaf]}\n")
for name, scafs in extra:
    bad = [x for x in scafs if x not in lengths or x in rename]
    if bad:
        sys.exit(f"ERROR: --extra-scope {name}: not an unmapped scaffold in the AGP: {bad}")
    with open(os.path.join(outdir, f"{name}.sizes"), "w") as fh:
        for scaf, chrom in chroms:
            fh.write(f"{chrom}\t{lengths[scaf]}\n")
        for scaf in scafs:
            fh.write(f"{scaf}\t{lengths[scaf]}\n")
with open(os.path.join(outdir, "rename.tsv"), "w") as fh:
    for scaf, chrom in chroms:
        fh.write(f"{scaf}\t{chrom}\n")
    for scaf in others:
        fh.write(f"{scaf}\t{scaf}\n")
print(f"{len(lengths)} scaffolds ({sum(lengths.values())} bp); {len(chroms)} mapped "
      f"({sum(lengths[s] for s, _ in chroms)} bp); {len(others)} unmapped", file=sys.stderr)
