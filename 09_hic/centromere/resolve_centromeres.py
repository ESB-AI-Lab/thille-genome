#!/usr/bin/env python3
"""Resolve each chromosome's centromere by building candidate regions from all evidence and eliminating them.

Candidates (per chromosome, whole length):
  - every refined Hi-C trans-contact peak (coarse-to-fine search from call_centromeres.py; all local maxima above
    median + 0.5 * (max - median), not only the top one)
  - every satellite cluster: TRF arrays with period >= --min-period, merged when closer than --merge-gap, >= --min-cluster
  Candidates that overlap (or lie within --join-mb, default 0) are merged into one candidate region. A Hi-C peak counts
  for a region when it lies within --hic-near-mb of it.

Centromeric satellite family (genome-wide, data-driven): TRF periods >= --min-period are grouped into families (periods
within --family-tol of the family's seed, seeds taken in order of total bp); the centromeric family is the one with the
most array bp inside the gene-poorest --gene-poor-frac of --window-res windows pooled over all chromosomes.

Elimination, in order (each decision recorded per candidate):
  1 telomere:        within --end-mb of an end that carries >= --min-telomere-bp of TTTAGGG/CCCTAAA array
  2 no_cen_satellite: < --min-cen-bp of the centromeric family (skipped for a chromosome if no candidate has it)
  3 gene_rich:       gene density above the chromosome median (skipped if it would remove every candidate)
  4 resolved:        --rank hic: highest Hi-C peak height, then most centromeric-satellite bp;
                     --rank satellite: most centromeric-satellite bp, then Hi-C height.
                     Survivors whose first rank key is within 10% of the winner's are flagged ambiguous.
Centromere region = winning candidate span. Point: --rank hic -> its Hi-C peak (else satellite midpoint);
--rank satellite -> its Hi-C peak if that lies within the largest centromeric cluster +/- --hic-near-mb, else the
midpoint of that cluster.

Usage:
  resolve_centromeres.py <chr12.mcool> <genes.gtf> <trf_dir> <out_prefix>
Writes <out_prefix>.bed, <out_prefix>.peak.bed, <out_prefix>.tsv, <out_prefix>_candidates.tsv, <out_prefix>.telomeres.tsv
and <out_prefix>.family.tsv.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

import cooler
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from call_centromeres import HicProfile, search, telomere_bp  # noqa: E402


def merge(intervals, gap=0):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered_bp(intervals, lo, hi):
    return sum(max(0, min(e, hi) - max(s, lo)) for s, e in merge(intervals))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mcool")
    ap.add_argument("gtf")
    ap.add_argument("trf_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("--coarse-res", type=int, default=1000000)
    ap.add_argument("--fine-res", default="250000,100000")
    ap.add_argument("--sigma-mb", type=float, default=1.0)
    ap.add_argument("--refine-mb", type=float, default=1.0)
    ap.add_argument("--min-period", type=int, default=100, help="minimum TRF period for satellites (default 100)")
    ap.add_argument("--merge-gap", type=int, default=50000, help="merge satellite arrays closer than this (default 50000)")
    ap.add_argument("--min-cluster", type=int, default=200000, help="minimum satellite cluster span (default 200000)")
    ap.add_argument("--join-mb", type=float, default=0.0, help="merge candidates within this distance (default 0 = overlap only)")
    ap.add_argument("--hic-near-mb", type=float, default=1.0, help="Hi-C peak association distance (default 1 Mb)")
    ap.add_argument("--rank", choices=("satellite", "hic"), required=True, help="step-4 primary evidence")
    ap.add_argument("--family-tol", type=float, default=0.05, help="period family tolerance (default 0.05)")
    ap.add_argument("--window-res", type=int, default=500000, help="gene-density window (default 500000)")
    ap.add_argument("--gene-poor-frac", type=float, default=0.10, help="gene-poorest window fraction (default 0.10)")
    ap.add_argument("--end-mb", type=float, default=2.0)
    ap.add_argument("--min-telomere-bp", type=int, default=1000)
    ap.add_argument("--telomere-window", type=int, default=100000)
    ap.add_argument("--min-cen-bp", type=int, default=100000, help="centromeric-family bp needed to keep a candidate (default 100000)")
    a = ap.parse_args()

    probe = cooler.Cooler(f"{a.mcool}::resolutions/{a.coarse_res}")
    chroms = [c for c in probe.chromnames if c.startswith("chr")]
    L = {c: int(probe.chromsizes[c]) for c in chroms}
    coarse = HicProfile(a.mcool, a.coarse_res, a.sigma_mb, chroms)
    fines = [HicProfile(a.mcool, int(r), a.sigma_mb, chroms) for r in a.fine_res.split(",")]

    # ---- sequence evidence
    arrays = {}
    for c in chroms:
        path = os.path.join(a.trf_dir, f"{c}.dat")
        if not os.path.exists(path + ".done"):
            raise SystemExit(f"ERROR: TRF output {path} is missing or unfinished")
        rows = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                f = line.split()
                if int(f[2]) >= a.min_period:
                    rows.append((int(f[0]) - 1, int(f[1]), int(f[2])))
        arrays[c] = rows
    telo = {c: telomere_bp(os.path.join(a.trf_dir, f"{c}.dat"), L[c], a.telomere_window) for c in chroms}
    gene_starts = defaultdict(list)
    with open(a.gtf) as fh:
        for line in fh:
            f = line.split("\t", 5)
            if len(f) > 4 and f[2] == "gene" and f[0] in L:
                gene_starts[f[0]].append(int(f[3]) - 1)
    for c in chroms:
        gene_starts[c] = np.sort(np.array(gene_starts[c]))

    # ---- centromeric satellite family
    by_period = defaultdict(int)
    for c in chroms:
        for s, e, p in arrays[c]:
            by_period[p] += e - s
    families, left = [], dict(by_period)
    while left:
        seed = max(left, key=left.get)
        members = {p for p in left if abs(p - seed) <= a.family_tol * seed}
        families.append((seed, members))
        for p in members:
            del left[p]
    windows = []
    for c in chroms:
        nw = int(np.ceil(L[c] / a.window_res))
        counts = np.bincount(gene_starts[c] // a.window_res, minlength=nw)[:nw]
        windows += [(counts[k], c, k * a.window_res, min(L[c], (k + 1) * a.window_res)) for k in range(nw)]
    windows.sort(key=lambda w: w[0])
    poor = windows[:max(1, int(round(a.gene_poor_frac * len(windows))))]
    fam_bp = []
    for seed, members in families:
        tot = 0
        for _, c, s, e in poor:
            tot += covered_bp([(x, y) for x, y, p in arrays[c] if p in members], s, e)
        fam_bp.append((tot, seed, members))
    fam_bp.sort(reverse=True)
    cen_seed, cen_members = fam_bp[0][1], fam_bp[0][2]
    with open(a.out_prefix + ".family.tsv", "w") as out:
        out.write("family_seed_period\tperiod_min\tperiod_max\ttotal_bp\tbp_in_gene_poor_windows\tcentromeric\n")
        for tot, seed, members in fam_bp[:15]:
            out.write(f"{seed}\t{min(members)}\t{max(members)}\t{sum(by_period[p] for p in members)}\t{tot}\t"
                      f"{'yes' if seed == cen_seed else 'no'}\n")

    # ---- candidates and elimination
    cand_rows, res_rows, tel_rows = [], [], []
    for c in chroms:
        peaks, P, med = search(c, coarse, fines, a.refine_mb)
        sat = merge([(s, e) for s, e, _ in arrays[c]], a.merge_gap)
        sat = [(s, e) for s, e in sat if e - s >= a.min_cluster]
        raw = [(p, p, "hic") for p, _ in peaks] + [(s, e, "sat") for s, e in sat]
        regions = []
        for s, e, src in sorted(raw):
            if regions and s <= regions[-1]["end"] + a.join_mb * 1e6:
                regions[-1]["end"] = max(regions[-1]["end"], e)
                regions[-1]["src"].add(src)
            else:
                regions.append({"start": s, "end": e, "src": {src}})
        gw = L[c] / 1e6
        chrom_density = len(gene_starts[c]) / gw
        cen_arrays = [(s, e) for s, e, p in arrays[c] if p in cen_members]
        for r in regions:
            s, e = int(r["start"]), int(r["end"])
            near = a.hic_near_mb * 1e6
            hp = [(p, h) for p, h in peaks if s - near <= p <= e + near]
            r["hic_pos"], r["hic_h"] = max(hp, key=lambda x: x[1]) if hp else (None, 0.0)
            r["cen_bp"] = covered_bp(cen_arrays, s, e + 1)
            r["sat_bp"] = covered_bp([(x, y) for x, y, _ in arrays[c]], s, e + 1)
            mid, half = (s + e) / 2, max((e - s) / 2, 5e5)
            lo, hi = max(0, mid - half), min(L[c], mid + half)
            n_genes = np.searchsorted(gene_starts[c], hi) - np.searchsorted(gene_starts[c], lo)
            r["gene_density"] = n_genes / ((hi - lo) / 1e6)
            r["dist_left"], r["dist_right"] = s, L[c] - e
            per = defaultdict(int)
            for x, y, p in arrays[c]:
                ov = min(y, e + 1) - max(x, s)
                if ov > 0:
                    per[p] += ov
            r["dominant_period"] = max(per, key=per.get) if per else "NA"
            r["status"] = "candidate"
        # 1 telomere
        for r in regions:
            if (r["dist_left"] <= a.end_mb * 1e6 and telo[c][0] >= a.min_telomere_bp) or \
               (r["dist_right"] <= a.end_mb * 1e6 and telo[c][1] >= a.min_telomere_bp):
                r["status"] = "eliminated:1_telomere"
        alive = [r for r in regions if r["status"] == "candidate"]
        flags = []
        # 2 centromeric satellite
        if any(r["cen_bp"] >= a.min_cen_bp for r in alive):
            for r in alive:
                if r["cen_bp"] < a.min_cen_bp:
                    r["status"] = "eliminated:2_no_centromeric_satellite"
        else:
            flags.append("no_centromeric_satellite")
        alive = [r for r in regions if r["status"] == "candidate"]
        # 3 gene-rich
        rich = [r for r in alive if r["gene_density"] > chrom_density]
        if rich and len(rich) < len(alive):
            for r in rich:
                r["status"] = "eliminated:3_gene_rich"
        elif rich:
            flags.append("all_candidates_gene_rich")
        alive = [r for r in regions if r["status"] == "candidate"]
        # 4 resolve
        if not alive:
            flags.append("no_candidate")
            win = None
        else:
            key = (lambda r: (r["hic_h"], r["cen_bp"])) if a.rank == "hic" else (lambda r: (r["cen_bp"], r["hic_h"]))
            alive.sort(key=key, reverse=True)
            win = alive[0]
            for r in alive[1:]:
                r["status"] = "eliminated:4_lower_rank"
                if key(win)[0] > 0 and key(r)[0] >= 0.9 * key(win)[0]:
                    flags.append("ambiguous")
            win["status"] = "centromere"
        for i, r in enumerate(regions):
            cand_rows.append((c, i + 1, int(r["start"]), int(r["end"]), "+".join(sorted(r["src"])),
                              "NA" if r["hic_pos"] is None else int(r["hic_pos"]), round(r["hic_h"] / med, 2),
                              r["cen_bp"], r["sat_bp"], r["dominant_period"], round(r["gene_density"], 1),
                              round(chrom_density, 1), int(r["dist_left"]), int(r["dist_right"]), r["status"]))
        for side, bp in (("left", telo[c][0]), ("right", telo[c][1])):
            hic = any(r["status"] == "eliminated:1_telomere" and "hic" in r["src"] and
                      (r["dist_left"] if side == "left" else r["dist_right"]) <= a.end_mb * 1e6 for r in regions)
            tel_rows.append((c, side, bp, "yes" if bp >= a.min_telomere_bp else "no", "yes" if hic else "no"))
        if win is None:
            res_rows.append(dict(chrom=c, length=L[c], flag=",".join(flags)))
            continue
        big = max(((x, y) for x, y in merge(cen_arrays, a.merge_gap) if x <= win["end"] and y >= win["start"]),
                  key=lambda t: t[1] - t[0], default=None)
        near = a.hic_near_mb * 1e6
        hic_in_sat = (win["hic_pos"] is not None and big is not None and big[0] - near <= win["hic_pos"] <= big[1] + near)
        if win["hic_pos"] is not None and (a.rank == "hic" or hic_in_sat or big is None):
            point, basis = int(win["hic_pos"]), "hic_peak"
        elif big is not None:
            point, basis = int((big[0] + big[1]) / 2), "satellite_midpoint"
        else:
            point, basis = int((win["start"] + win["end"]) / 2), "region_midpoint"
        if win["hic_pos"] is None:
            flags.append("no_hic_peak")
        elif basis != "hic_peak":
            flags.append("hic_peak_outside_satellite")
        res_rows.append(dict(
            chrom=c, length=L[c], region_start=int(win["start"]), region_end=int(win["end"]) + 1,
            region_mb=f"{win['start'] / 1e6:.2f}-{(win['end'] + 1) / 1e6:.2f}", point_bp=point,
            point_mb=round(point / 1e6, 2), point_basis=basis, hic_score=round(win["hic_h"] / med, 2),
            centromeric_satellite_bp=win["cen_bp"], dominant_period=win["dominant_period"],
            gene_density_per_mb=round(win["gene_density"], 1), chrom_gene_density_per_mb=round(chrom_density, 1),
            n_candidates=len(regions), telomere_left_bp=telo[c][0], telomere_right_bp=telo[c][1],
            flag=",".join(sorted(set(flags))) or "ok"))

    cols = list(max(res_rows, key=len).keys())
    with open(a.out_prefix + ".tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in res_rows:
            out.write("\t".join(str(r.get(k, "NA")) for k in cols) + "\n")
    with open(a.out_prefix + "_candidates.tsv", "w") as out:
        out.write("chrom\tcandidate\tstart\tend\tevidence\thic_peak_bp\thic_score\tcentromeric_satellite_bp\t"
                  "satellite_bp\tdominant_period\tgene_density_per_mb\tchrom_gene_density_per_mb\tdist_left\t"
                  "dist_right\tstatus\n")
        for r in cand_rows:
            out.write("\t".join(map(str, r)) + "\n")
    with open(a.out_prefix + ".bed", "w") as ob, open(a.out_prefix + ".peak.bed", "w") as op:
        for r in res_rows:
            if "point_bp" not in r:
                continue
            num = re.sub(r"^chr", "", r["chrom"])
            score = min(1000, int(round(max(r["hic_score"], 0) * 100)))
            ob.write(f"{r['chrom']}\t{r['region_start']}\t{r['region_end']}\tcen{num}\t{score}\t.\n")
            op.write(f"{r['chrom']}\t{r['point_bp']}\t{r['point_bp'] + 1}\tcen{num}\t{score}\t.\n")
    with open(a.out_prefix + ".telomeres.tsv", "w") as out:
        out.write("chrom\tend\ttelomere_repeat_bp\ttelomere_repeat_present\thic_telomere_peak\n")
        for r in tel_rows:
            out.write("\t".join(map(str, r)) + "\n")
    print(f"centromeric family seed period {cen_seed} ({min(cen_members)}-{max(cen_members)} bp)")
    for r in res_rows:
        print(r["chrom"], r.get("region_mb"), r.get("point_mb"), r.get("point_basis"), r["flag"])


if __name__ == "__main__":
    main()
