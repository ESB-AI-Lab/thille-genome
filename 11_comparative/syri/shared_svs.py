#!/usr/bin/env python3
"""SVs shared by Puebla and the Unnamed Lowland genome relative to Thille (v4 SyRI runs; Thille is the reference in both).

Events of the same class are matched on Thille coordinates by reciprocal overlap >= --ro (insertions, which have no
Thille span, by breakpoint distance <= --ins-dist). Classes: inversion (INV), translocation (TRANS, INVTR),
duplication (DUP, INVDP), deletion (DEL), insertion (INS), copy gain (CPG), copy loss (CPL), Thille-only sequence
(NOTAL with a Thille interval). Also reports whether Hass UoQ carries the same event (same rule).
Usage (from assemblies/thille.hic): compare/v4/shared_svs.py [--min-size 50] [--ro 0.5] [--ins-dist 100]
Writes compare/v4/shared_svs.tsv and prints counts by class and size.
"""
import argparse
import bisect
import collections

CLS = {"INV": "inversion", "TRANS": "translocation", "INVTR": "translocation", "DUP": "duplication",
       "INVDP": "duplication", "DEL": "deletion", "INS": "insertion", "CPG": "copy gain", "CPL": "copy loss",
       "NOTAL": "Thille-only sequence"}


def events(pair, min_size):
    out = collections.defaultdict(list)
    for l in open(f"compare/v4/syri/{pair}/{pair}syri.out"):
        f = l.rstrip("\n").split("\t")
        t = f[10]
        if t not in CLS or f[1] == "-":
            continue
        s, e = int(f[1]) - 1, int(f[2])
        q = int(f[7]) - int(f[6]) + 1 if f[6] not in ("-", "") else 0
        size = max(e - s, q)
        if t == "NOTAL" and f[5] != "-":
            continue  # NOTAL rows with a query interval are query-only sequence
        if size < min_size:
            continue
        out[(CLS[t], f[0])].append((s, e, size, f[8]))
    for v in out.values():
        v.sort()
    return out


def match(a, pool, cls, ro, ins_dist):
    s, e = a[0], a[1]
    starts = [x[0] for x in pool]
    i = bisect.bisect_left(starts, s - max(ins_dist, 10_000_000))
    best = None
    for x in pool[i:]:
        if x[0] > e + ins_dist:
            break
        if cls == "insertion":
            ok = abs(x[0] - s) <= ins_dist
        else:
            ov = min(e, x[1]) - max(s, x[0])
            ok = ov > 0 and ov >= ro * (e - s) and ov >= ro * (x[1] - x[0])
        if ok:
            best = x
            break
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--ro", type=float, default=0.5)
    ap.add_argument("--ins-dist", type=int, default=100)
    a = ap.parse_args()
    ev = {p: events(p, a.min_size) for p in ("westindian", "puebla", "hass")}
    counts = collections.Counter()
    with open("compare/v4/shared_svs.tsv", "w") as fh:
        fh.write("class\tthille_chrom\tthille_start\tthille_end\tlowland_id\tlowland_size\tpuebla_id\tpuebla_size\t"
                 "also_in_hass\n")
        for (cls, chrom), lst in sorted(ev["westindian"].items()):
            for x in lst:
                bin_ = "<1 kb" if x[2] < 1000 else ("1-100 kb" if x[2] < 100000 else ">=100 kb")
                counts[(cls, bin_, "lowland")] += 1
                y = match(x, ev["puebla"].get((cls, chrom), []), cls, a.ro, a.ins_dist)
                if y:
                    h = match(x, ev["hass"].get((cls, chrom), []), cls, a.ro, a.ins_dist)
                    counts[(cls, bin_, "shared")] += 1
                    counts[(cls, bin_, "shared+hass")] += bool(h)
                    fh.write(f"{cls}\t{chrom}\t{x[0]}\t{x[1]}\t{x[3]}\t{x[2]}\t{y[3]}\t{y[2]}\t{'yes' if h else 'no'}\n")
        for (cls, chrom), lst in ev["puebla"].items():
            for x in lst:
                bin_ = "<1 kb" if x[2] < 1000 else ("1-100 kb" if x[2] < 100000 else ">=100 kb")
                counts[(cls, bin_, "puebla")] += 1
    print("class\tsize\tLowland_vs_Thille\tPuebla_vs_Thille\tshared_Lowland_and_Puebla\tof_shared_also_in_Hass")
    for cls in dict.fromkeys(CLS.values()):
        for b in ("<1 kb", "1-100 kb", ">=100 kb"):
            if counts[(cls, b, "lowland")] or counts[(cls, b, "puebla")]:
                print(f"{cls}\t{b}\t{counts[(cls, b, 'lowland')]}\t{counts[(cls, b, 'puebla')]}\t"
                      f"{counts[(cls, b, 'shared')]}\t{counts[(cls, b, 'shared+hass')]}")


if __name__ == "__main__":
    main()
