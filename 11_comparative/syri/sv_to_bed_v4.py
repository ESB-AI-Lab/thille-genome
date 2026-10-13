#!/usr/bin/env python3
# v4 copy (2026-10-05): v4 SyRI, Gwen SV, annotation and centromere inputs; outputs under compare/v4/.
# 2026-10-07: Puebla added (its SyRI against Thille existed but the v3 genome list was kept when porting), and the
# Hass (reference) vs Puebla SyRI written as its own pair of BEDs, both sides in original GenBank coordinates.
"""Write the structural variants between Thille and each other assembly as BED files, one per assembly, each in that
assembly's own coordinates.

SyRI (Hass, West Indian, Puebla): every structural annotation in <name>syri.out except SNPs, syntenic regions and the
per-alignment rows (*AL). SyRI coordinates are 1-based inclusive; BED is 0-based half-open. The query FASTA given to
SyRI had chromosomes renamed chr1..chr12 and some reverse-complemented (see sbatch_syri.sh), so query coordinates
are mapped back to the original GenBank ids and forward-strand positions. Insertions have no Thille span (and
deletions no query span); those rows are written as the 1-bp flanking position SyRI reports, with size_bp giving the
event size on the other genome. Rows where the assembly has no interval (NOTAL on the other genome) are omitted from
that assembly's file.

Gwen: rows of gwen_sv_candidates.py (Thille coordinates only; Gwen is contig-level).

Columns (BED6+): chrom start end name score strand type size_bp partner_genome partner_chrom partner_start
partner_end detail. name = <genome>_<SyRI id or kind>; strand '-' marks inverted events; size_bp = larger of the two
spans. All-size files are *.all.bed; >= --min-size files are *.ge100kb.bed.

Usage: sv_to_bed.py <outdir> [--min-size 100000]
"""
import argparse
import os

SKIP = {"SNP", "SYN"}
INVERTED = {"INV", "INVTR", "INVDP"}
# per-alignment rows inside SYN/INV/TRANS/DUP blocks; NOTAL also ends in "AL" but is a real annotation
ALIGNMENT_ROWS = {"SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"}
QUERIES = {
    "hass": ("hass_UoQ_chr_only.fasta.fai",
             ["CM056809.1", "CM056810.1", "CM056811.1", "CM056812.1", "CM056813.1", "CM056814.1",
              "CM056815.1", "CM056816.1", "CM056817.1", "CM056818.1", "CM056819.1", "CM056820.1"]),
    "westindian": ("westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna.fai",
                   ["CP152278.1", "CP152279.1", "CP152280.1", "CP152281.1", "CP152282.1", "CP152283.1",
                    "CP152284.1", "CP152285.1", "CP152286.1", "CP152287.1", "CP152288.1", "CP152289.1"]),
    "puebla": ("compare/v4/puebla/puebla.chr.orig.fasta.fai",
               ["CM191109.1", "CM191110.1", "CM191111.1", "CM191112.1", "CM191113.1", "CM191114.1",
                "CM191115.1", "CM191116.1", "CM191117.1", "CM191118.1", "CM191119.1", "CM191120.1"]),
}
# SyRI runs whose reference is not Thille: name -> (reference genome, query genome); both FASTAs were renamed and
# strand-flipped the same way (hass.chr.samestrand.fasta and puebla.chr.samestrand.fasta)
OTHER_PAIRS = {"hass_puebla": ("hass", "puebla")}
RC = {"chr1", "chr4", "chr5", "chr6", "chr10", "chr11"}   # reverse-complemented in sbatch_syri.sh


def to_orig(chrom, s, e, ids, lens):
    """1-based inclusive renamed-query interval -> (original id, 0-based start, end)."""
    oid = ids[int(chrom[3:]) - 1]
    if chrom in RC:
        L = lens[oid]
        s, e = L - e + 1, L - s + 1
    return oid, s - 1, e


def syri_rows(name, fai, ids, ref_fai=None, ref_ids=None):
    lens = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
    ref_lens = {l.split()[0]: int(l.split()[1]) for l in open(ref_fai)} if ref_fai else None
    for line in open(f"compare/v4/syri/{name}/{name}syri.out"):
        f = line.rstrip("\n").split("\t")
        typ = f[10]
        if typ in SKIP or typ in ALIGNMENT_ROWS:
            continue
        rc, rs, re_, qc, qs, qe = f[0], f[1], f[2], f[5], f[6], f[7]
        if rs == "-":
            ref = None
        elif ref_ids:
            ref = to_orig(rc, min(int(rs), int(re_)), max(int(rs), int(re_)), ref_ids, ref_lens)
        else:
            ref = (rc, int(rs) - 1, int(re_))
        qry = None if qs == "-" or qc == "-" else to_orig(qc, min(int(qs), int(qe)), max(int(qs), int(qe)), ids, lens)
        if ref and ref[2] <= ref[1]:
            ref = (ref[0], ref[2] - 1, ref[1] + 1)
        size = max(ref[2] - ref[1] if ref else 0, qry[2] - qry[1] if qry else 0)
        strand = "-" if typ in INVERTED else "."
        detail = f[11] if f[11] != "-" else ""
        yield typ, f[8], ref, qry, size, strand, detail


def fmt(chrom, s, e, name, strand, typ, size, pg, partner, detail):
    pc, ps, pe = partner if partner else (".", ".", ".")
    return f"{chrom}\t{s}\t{e}\t{name}\t0\t{strand}\t{typ}\t{size}\t{pg}\t{pc}\t{ps}\t{pe}\t{detail or '.'}\n"


HEADER = ("#chrom\tstart\tend\tname\tscore\tstrand\ttype\tsize_bp\tpartner_genome\tpartner_chrom\tpartner_start\t"
          "partner_end\tdetail\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--min-size", type=int, default=100000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    thille = {"all": [], "big": []}
    for name, (fai, ids) in QUERIES.items():
        own = {"all": [], "big": []}
        n = 0
        for typ, sid, ref, qry, size, strand, detail in syri_rows(name, fai, ids):
            n += 1
            keys = ["all"] + (["big"] if size >= a.min_size else [])
            for k in keys:
                if ref:
                    thille[k].append((ref[0], ref[1], fmt(*ref, f"{name}_{sid}", strand, typ, size, name, qry, detail)))
                if qry:
                    own[k].append((qry[0], qry[1], fmt(*qry, f"thille_{sid}", strand, typ, size, "thille", ref, detail)))
        for k, suf in (("all", "all"), ("big", "ge100kb")):
            path = f"{a.outdir}/{name}.sv.{suf}.bed"
            with open(path, "w") as fh:
                fh.write(HEADER)
                fh.writelines(r[2] for r in sorted(own[k], key=lambda r: (r[0], r[1])))
            print(f"{path}: {len(own[k])} rows")
        print(f"  ({n} SyRI structural rows read for {name})")
    for name, (rg, qg) in OTHER_PAIRS.items():
        (rfai, rids), (qfai, qids) = QUERIES[rg], QUERIES[qg]
        side = {rg: {"all": [], "big": []}, qg: {"all": [], "big": []}}
        n = 0
        for typ, sid, ref, qry, size, strand, detail in syri_rows(name, qfai, qids, rfai, rids):
            n += 1
            for k in ["all"] + (["big"] if size >= a.min_size else []):
                if ref:
                    side[rg][k].append((ref[0], ref[1], fmt(*ref, f"{qg}_{sid}", strand, typ, size, qg, qry, detail)))
                if qry:
                    side[qg][k].append((qry[0], qry[1], fmt(*qry, f"{rg}_{sid}", strand, typ, size, rg, ref, detail)))
        for g in (rg, qg):
            for k, suf in (("all", "all"), ("big", "ge100kb")):
                path = f"{a.outdir}/{name}.{g}_coords.sv.{suf}.bed"
                with open(path, "w") as fh:
                    fh.write(HEADER)
                    fh.writelines(r[2] for r in sorted(side[g][k], key=lambda r: (r[0], r[1])))
                print(f"{path}: {len(side[g][k])} rows")
        print(f"  ({n} SyRI structural rows read for {name})")
    # Gwen candidates, Thille coordinates only
    gwen = []
    for line in open("compare/v4/gwen_sv/gwen_vs_thille.sv100kb.tsv"):
        if line.startswith("chrom\t"):
            continue
        f = line.rstrip("\n").split("\t")
        typ = {"uncovered": "GWEN_UNCOVERED", "multi": "GWEN_MULTICONTIG"}[f[4]]
        detail = f"read_class={f[10]};depth_ratio={f[9]};contigs={f[6] or '.'};at_chrom_end={f[11]}"
        row = (f[0], int(f[1]), fmt(f[0], f[1], f[2], f"gwen_{f[4]}_{f[0]}_{f[1]}", ".", typ, f[3], "gwen", None, detail))
        gwen.append(row)
        thille["all"].append(row); thille["big"].append(row)
    path = f"{a.outdir}/gwen.sv.ge100kb.thille_coords.bed"
    with open(path, "w") as fh:
        fh.write(HEADER); fh.writelines(r[2] for r in sorted(gwen, key=lambda r: (r[0], r[1])))
    print(f"{path}: {len(gwen)} rows")
    order = {f"chr{i}": i for i in range(1, 13)}
    for k, suf in (("all", "all"), ("big", "ge100kb")):
        path = f"{a.outdir}/thille.sv.{suf}.bed"
        with open(path, "w") as fh:
            fh.write(HEADER)
            fh.writelines(r[2] for r in sorted(thille[k], key=lambda r: (order.get(r[0], 99), r[1])))
        print(f"{path}: {len(thille[k])} rows")


if __name__ == "__main__":
    main()
