#!/usr/bin/env python3
"""Plot windows for every SyRI v4 inversion between Thille and Hass UoQ, the Unnamed Lowland genome and Puebla.

Each INV row of compare/v4/syri/<pair>/<pair>syri.out gives the inversion on Thille (reference) and on the compared
genome (query, renamed chr1..chr12 on Thille's strand). Query coordinates are converted back to the original
accession and strand, so they match the BAMs and annotations of that genome. Per genome, an inversion of at most
WIN bp gets one WIN-bp window centered on it; a larger one gets two WIN-bp windows, one centered on each breakpoint.
Windows that would run past a chromosome end are shifted inward so every window stays exactly WIN bp.
Coordinates are 0-based half-open.

Usage (from assemblies/thille.hic): compare/v4/inversions/make_windows.py [--win 100000] [--min-size 100]
Writes compare/v4/inversions/inversions.tsv and windows.tsv.
"""
import argparse
import json

THILLE_FAI = "thille.hifiv4.asm.hic.scaffold.chr.only.fasta.fai"
RC_DEFAULT = {"chr1", "chr4", "chr5", "chr6", "chr10", "chr11"}  # reverse-complemented in compare/sbatch_syri_v4.sh
PAIRS = {
    "hass": dict(fai="hass_UoQ_chr_only.fasta.fai",
                 ids=[f"CM0568{9 + i:02d}.1" for i in range(12)]),
    "westindian": dict(fai="westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna.fai",
                       ids=[f"CP1522{78 + i}.1" for i in range(12)]),
    "puebla": dict(fai="/gpool/edwin/avocados/thille_paper/genome/puebla/GCA_060308345.1/ncbi_dataset/data/"
                       "GCA_060308345.1/GCA_060308345.1_ASM6030834v1_genomic.fna.fai",
                   ids=[f"CM1911{9 + i:02d}.1" for i in range(12)]),
}
OUT = "compare/v4/inversions"


def lengths(fai):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}


def puebla_rc():
    rc = set()
    for line in list(open("compare/v4/puebla/puebla_strand.tsv"))[1:]:
        f = line.rstrip("\n").split("\t")
        if f[5] == "yes":
            rc.add(f[1])
    return rc


def windows(s, e, L, win):
    """WIN-bp windows for interval [s, e) on a chromosome of length L."""
    def place(center):
        a = max(0, min(center - win // 2, L - win))
        return a, a + win
    if e - s <= win:
        return [("whole", *place((s + e) // 2))]
    return [("left", *place(s)), ("right", *place(e))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--win", type=int, default=100000)
    ap.add_argument("--min-size", type=int, default=100)
    a = ap.parse_args()
    tl = lengths(THILLE_FAI)
    invs, wins = [], []
    for pair, cfg in PAIRS.items():
        ql = lengths(cfg["fai"])
        rc = puebla_rc() if pair == "puebla" else RC_DEFAULT
        ids = cfg["ids"]
        assert all(i in ql for i in ids), (pair, [i for i in ids if i not in ql])
        for line in open(f"compare/v4/syri/{pair}/{pair}syri.out"):
            f = line.rstrip("\n").split("\t")
            if f[10] != "INV":
                continue
            ts, te = int(f[1]) - 1, int(f[2])
            qchr, qs, qe = f[5], int(f[6]), int(f[7])
            oid = ids[int(qchr[3:]) - 1]
            L = ql[oid]
            if qchr in rc:
                qs, qe = L - qe + 1, L - qs + 1
            qs -= 1
            if max(te - ts, qe - qs) < a.min_size:
                continue
            iid = f"{pair}_{f[8]}"
            invs.append([iid, pair, f[0], ts, te, te - ts, oid, qs, qe, qe - qs, "yes" if qchr in rc else "no"])
            for genome, chrom, s, e, LL in (("thille", f[0], ts, te, tl[f[0]]), (pair, oid, qs, qe, L)):
                for side, ws, we in windows(s, e, LL, a.win):
                    assert we - ws == a.win
                    wins.append([iid, pair, genome, side, chrom, ws, we, s, e])
    with open(f"{OUT}/inversions.tsv", "w") as fh:
        fh.write("inv_id\tpair\tthille_chrom\tthille_start\tthille_end\tthille_len\tother_chrom\tother_start\t"
                 "other_end\tother_len\tother_chrom_revcomp_for_syri\n")
        for r in invs:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(f"{OUT}/windows.tsv", "w") as fh:
        fh.write("inv_id\tpair\tgenome\tside\tchrom\twin_start\twin_end\tinv_start\tinv_end\n")
        for r in wins:
            fh.write("\t".join(map(str, r)) + "\n")
    by = {}
    for r in invs:
        by[r[1]] = by.get(r[1], 0) + 1
    print(f"inversions: {len(invs)} {json.dumps(by)}; windows: {len(wins)} "
          f"(Thille {sum(1 for w in wins if w[2] == 'thille')}, other {sum(1 for w in wins if w[2] != 'thille')})")


if __name__ == "__main__":
    main()
