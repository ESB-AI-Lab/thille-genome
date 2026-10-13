#!/usr/bin/env python3
"""Call each chromosome's centromere (and telomere contact peaks) from Hi-C, searching the whole chromosome
coarse to fine, with supporting sequence tracks.

Hi-C signal. In many plant nuclei centromeres cluster and telomeres cluster, so a chromosome's centromere and its
telomeres each form contact hotspots with the other chromosomes. For chromosome i, every inter-chromosome block with
chromosome j is ICE-balanced, divided by its mean non-zero value and Gaussian-smoothed (sigma = --sigma-mb), and the
contact profile is
    P_i(x) = sum_{j != i} max_y OE_ij(x, y)
i.e. each position's strongest hotspot with every other chromosome, wherever it lies along that chromosome.

Search (whole chromosome, no restriction to the interior):
  1. coarse: P_i at --coarse-res over the full chromosome; candidate areas are local maxima with
     P >= median + 0.5 * (max - median);
  2. fine: each candidate is re-evaluated at each --fine-res in turn over the neighbouring coarse bins (+/- --refine-mb),
     keeping the refined peak.
Classification: a refined peak within --end-mb of an end is a telomere peak when that end carries at least
--min-telomere-bp of TTTAGGG/CCCTAAA array (TRF period 7) within --telomere-window of the end; the centromere is the
highest remaining peak (it may lie near an end that has no telomere array, e.g. acrocentric). Interval = contiguous
positions around the centromere peak (finest resolution, +/- --interval-mb) with P >= median + 0.5 * (peak - median);
score = coarse peak / coarse median. The call is repeated with each --stability-res as the coarse resolution and the
largest shift of the final peak is reported.

Supporting tracks (--track-res bins, --smooth-bins moving average, --end-mb excluded at both ends): gene density
minimum, repeat density maximum, satellite (TRF period >= --min-period) density maximum, and the dominant satellite
period inside the interval. A track supports the call if it lies within the interval +/- --support-mb; fewer than 2 of 3
supporting -> flag low_support. Flags never move a call.

Usage:
  call_centromeres.py <chr12.mcool> <genes.gtf> <repeats.gff> <trf_dir> <out_prefix>
Writes <out_prefix>.bed (BED6 intervals), <out_prefix>.peak.bed (1-bp centromere peaks), <out_prefix>.tsv and
<out_prefix>.telomeres.tsv.
"""
import argparse
import os
import re
from collections import defaultdict

import cooler
import numpy as np
from scipy.ndimage import gaussian_filter

TELO_MOTIFS = set()
for m in ("TTTAGGG", "CCCTAAA"):
    for k in range(7):
        TELO_MOTIFS.add(m[k:] + m[:k])


class HicProfile:
    """Per-chromosome trans-contact profile at one resolution, computed for any row window on demand."""

    def __init__(self, mcool, res, sigma_mb, chroms):
        self.clr = cooler.Cooler(f"{mcool}::resolutions/{res}")
        self.res = res
        self.sigma = sigma_mb * 1e6 / res
        self.chroms = chroms
        self.nbins = {c: len(self.clr.bins().fetch(c)) for c in chroms}
        self._blocks = {}

    def _block(self, ci, cj):
        key = (ci, cj)
        if key not in self._blocks:
            B = np.nan_to_num(self.clr.matrix(balance=True).fetch(ci, cj)).astype(np.float32)
            nz = B[B > 0]
            if nz.size == 0:
                raise SystemExit(f"ERROR: no trans contacts between {ci} and {cj} at {self.res} bp")
            self._blocks[key] = gaussian_filter(B / nz.mean(), sigma=self.sigma)
        return self._blocks[key]

    def profile(self, ci, lo=0, hi=None):
        hi = self.nbins[ci] if hi is None else hi
        return sum(self._block(ci, cj)[lo:hi].max(axis=1) for cj in self.chroms if cj != ci)


def local_maxima(p):
    idx = []
    for k in range(len(p)):
        left = p[k - 1] if k > 0 else -np.inf
        right = p[k + 1] if k < len(p) - 1 else -np.inf
        if p[k] >= left and p[k] > right:
            idx.append(k)
    return idx


def search(chrom, coarse, fines, refine_mb):
    """Coarse candidates, each refined through the fine resolutions. Returns [(peak_bp, coarse_height)], coarse profile."""
    P = coarse.profile(chrom)
    med = float(np.median(P))
    thr = med + 0.5 * (P.max() - med)
    peaks = []
    for k in local_maxima(P):
        if P[k] < thr:
            continue
        centre = (k + 0.5) * coarse.res
        for f in fines:
            n = f.nbins[chrom]
            lo = max(0, int((centre - refine_mb * 1e6) // f.res))
            hi = min(n, int((centre + refine_mb * 1e6) // f.res) + 1)
            sub = f.profile(chrom, lo, hi)
            centre = (lo + int(np.argmax(sub)) + 0.5) * f.res
        peaks.append((centre, float(P[k])))
    return peaks, P, med


def telomere_bp(dat, length, window):
    left = right = 0
    with open(dat) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split()
            if int(f[2]) != 7 or f[13] not in TELO_MOTIFS:
                continue
            s, e = int(f[0]), int(f[1])
            if s <= window:
                left += e - s + 1
            if e >= length - window:
                right += e - s + 1
    return left, right


def merged_bp_per_bin(intervals, nbins, res):
    out = np.zeros(nbins)
    if not intervals:
        return out
    intervals.sort()
    cs, ce = intervals[0]
    merged = []
    for s, e in intervals[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    for s, e in merged:
        b = s // res
        while s < e and b < nbins:
            be = min(e, (b + 1) * res)
            out[b] += be - s
            s, b = be, b + 1
    return out


def smooth(v, k):
    return np.convolve(v, np.ones(k) / k, mode="same")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mcool")
    ap.add_argument("gtf")
    ap.add_argument("gff")
    ap.add_argument("trf_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("--coarse-res", type=int, default=1000000, help="coarse resolution (default 1000000)")
    ap.add_argument("--fine-res", default="250000,100000", help="refinement resolutions, coarse to fine (default 250000,100000)")
    ap.add_argument("--stability-res", default="500000,250000", help="alternative coarse resolutions for the stability check")
    ap.add_argument("--sigma-mb", type=float, default=1.0, help="Gaussian smoothing of contact blocks in Mb (default 1)")
    ap.add_argument("--refine-mb", type=float, default=1.0, help="refinement window around each candidate (default 1 Mb)")
    ap.add_argument("--interval-mb", type=float, default=5.0, help="window for the half-maximum interval (default 5 Mb)")
    ap.add_argument("--end-mb", type=float, default=2.0, help="distance from an end for telomere peaks/track exclusion (default 2)")
    ap.add_argument("--min-telomere-bp", type=int, default=1000, help="TTTAGGG array bp needed to call a telomere (default 1000)")
    ap.add_argument("--telomere-window", type=int, default=100000, help="bp from each end searched for telomere repeats (default 100000)")
    ap.add_argument("--track-res", type=int, default=500000)
    ap.add_argument("--smooth-bins", type=int, default=9)
    ap.add_argument("--min-period", type=int, default=100)
    ap.add_argument("--support-mb", type=float, default=5.0)
    a = ap.parse_args()

    probe = cooler.Cooler(f"{a.mcool}::resolutions/{a.coarse_res}")
    chroms = [c for c in probe.chromnames if c.startswith("chr")]
    L = {c: int(probe.chromsizes[c]) for c in chroms}
    fine_res = [int(x) for x in a.fine_res.split(",")]
    fines = [HicProfile(a.mcool, r, a.sigma_mb, chroms) for r in fine_res]
    finest = fines[-1]
    coarse_sets = {a.coarse_res: HicProfile(a.mcool, a.coarse_res, a.sigma_mb, chroms)}
    for r in (int(x) for x in a.stability_res.split(",") if x):
        coarse_sets[r] = HicProfile(a.mcool, r, a.sigma_mb, chroms)

    telo = {}
    for c in chroms:
        path = os.path.join(a.trf_dir, f"{c}.dat")
        if not os.path.exists(path + ".done"):
            raise SystemExit(f"ERROR: TRF output {path} is missing or unfinished")
        telo[c] = telomere_bp(path, L[c], a.telomere_window)

    def classify(c, peaks):
        end_bp = a.end_mb * 1e6
        tl, tr = telo[c]
        telomere_peaks, rest = [], []
        for pos, h in peaks:
            at_left, at_right = pos <= end_bp, pos >= L[c] - end_bp
            if (at_left and tl >= a.min_telomere_bp) or (at_right and tr >= a.min_telomere_bp):
                telomere_peaks.append((pos, h, "left" if at_left else "right"))
            else:
                rest.append((pos, h))
        cen = max(rest, key=lambda x: x[1]) if rest else None
        return cen, telomere_peaks, rest

    calls, tel_rows = {}, []
    for c in chroms:
        peaks, P, med = search(c, coarse_sets[a.coarse_res], fines, a.refine_mb)
        cen, tpeaks, rest = classify(c, peaks)
        alt = []
        for r, prof in coarse_sets.items():
            if r == a.coarse_res:
                continue
            ap_, _, _ = search(c, prof, fines, a.refine_mb)
            ac, _, _ = classify(c, ap_)
            alt.append(ac[0] if ac else None)
        for side, bp in (("left", telo[c][0]), ("right", telo[c][1])):
            has_peak = any(t[2] == side for t in tpeaks)
            near = any((p <= a.end_mb * 1e6) if side == "left" else (p >= L[c] - a.end_mb * 1e6) for p, _ in peaks)
            tel_rows.append((c, side, bp, "yes" if bp >= a.min_telomere_bp else "no", "yes" if near else "no",
                             "yes" if has_peak else "no"))
        calls[c] = (cen, rest, P, med, alt)

    # ---- supporting tracks
    tres = a.track_res
    nb = {c: int(np.ceil(L[c] / tres)) for c in chroms}
    genes = {c: np.zeros(nb[c]) for c in chroms}
    with open(a.gtf) as fh:
        for line in fh:
            f = line.split("\t", 5)
            if len(f) > 4 and f[2] == "gene" and f[0] in genes:
                genes[f[0]][min((int(f[3]) - 1) // tres, nb[f[0]] - 1)] += 1
    rep_iv = defaultdict(list)
    with open(a.gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if len(f) > 4 and f[0] in nb:
                rep_iv[f[0]].append((int(f[3]) - 1, int(f[4])))
    sat_iv, sat_rows = defaultdict(list), defaultdict(list)
    for c in chroms:
        with open(os.path.join(a.trf_dir, f"{c}.dat")) as fh:
            for line in fh:
                if line.startswith("@"):
                    continue
                f = line.split()
                s, e, period = int(f[0]) - 1, int(f[1]), int(f[2])
                if period >= a.min_period:
                    sat_iv[c].append((s, e))
                    sat_rows[c].append((s, e, period))
    end_b = int(round(a.end_mb * 1e6 / tres))

    rows = []
    for c in chroms:
        cen, rest, P, med, alt = calls[c]
        if cen is None:
            rows.append(dict(chrom=c, length=L[c], peak_mb="NA", flag="no_centromere_peak"))
            continue
        pk_bp, pk_h = cen
        n = finest.nbins[c]
        k = min(n - 1, int(pk_bp // finest.res))
        lo = max(0, k - int(a.interval_mb * 1e6 // finest.res))
        hi = min(n, k + int(a.interval_mb * 1e6 // finest.res) + 1)
        sub = finest.profile(c, lo, hi)
        fmed = float(np.median(finest.profile(c)))
        thr = fmed + 0.5 * (sub[k - lo] - fmed)
        s_, e_ = k - lo, k - lo
        while s_ > 0 and sub[s_ - 1] >= thr:
            s_ -= 1
        while e_ < len(sub) - 1 and sub[e_ + 1] >= thr:
            e_ += 1
        iv_s, iv_e = (lo + s_) * finest.res, min((lo + e_ + 1) * finest.res, L[c])
        shift = max((abs(x - pk_bp) for x in alt if x is not None), default=0.0)
        second = sorted((h for p, h in rest if abs(p - pk_bp) > a.refine_mb * 1e6), reverse=True)

        def pos(track, fn):
            v = smooth(track, a.smooth_bins)
            inner = v[end_b:len(v) - end_b]
            return (fn(inner) + end_b + 0.5) * tres / 1e6

        gmin = pos(genes[c], np.argmin)
        rmax = pos(merged_bp_per_bin(rep_iv[c], nb[c], tres), np.argmax)
        smax = pos(merged_bp_per_bin(sat_iv[c], nb[c], tres), np.argmax)
        per_bp = defaultdict(int)
        for s, e, p in sat_rows[c]:
            ov = min(e, iv_e) - max(s, iv_s)
            if ov > 0:
                per_bp[p] += ov
        dom = max(per_bp, key=per_bp.get) if per_bp else "NA"
        w0, w1 = iv_s / 1e6 - a.support_mb, iv_e / 1e6 + a.support_mb
        support = [w0 <= x <= w1 for x in (gmin, rmax, smax)]
        flags = []
        if sum(support) < 2:
            flags.append("low_support")
        if pk_bp <= a.end_mb * 1e6 or pk_bp >= L[c] - a.end_mb * 1e6:
            flags.append("near_chromosome_end")
        if second and second[0] >= 0.9 * pk_h:
            flags.append("ambiguous_second_peak")
        if shift > 1e6:
            flags.append("unstable_across_resolutions")
        rows.append(dict(
            chrom=c, length=L[c], peak_bp=int(pk_bp), peak_mb=round(pk_bp / 1e6, 2), interval_start=iv_s,
            interval_end=iv_e, interval_mb=f"{iv_s / 1e6:.2f}-{iv_e / 1e6:.2f}", score=round(pk_h / med, 2),
            n_candidate_peaks=len(rest), second_peak_ratio=round(second[0] / pk_h, 2) if second else "NA",
            max_shift_other_res_mb=round(shift / 1e6, 2),
            telomere_left_bp=telo[c][0], telomere_right_bp=telo[c][1],
            gene_min_mb=round(gmin, 2), repeat_max_mb=round(rmax, 2), satellite_max_mb=round(smax, 2),
            dominant_satellite_period=dom, gene_supports=support[0], repeat_supports=support[1],
            satellite_supports=support[2], n_supporting=sum(support), flag=",".join(flags) or "ok"))

    cols = list(max(rows, key=len).keys())
    with open(a.out_prefix + ".tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(k, "NA")) for k in cols) + "\n")
    with open(a.out_prefix + ".bed", "w") as ob, open(a.out_prefix + ".peak.bed", "w") as op:
        for r in rows:
            if r.get("peak_bp") is None:
                continue
            num = re.sub(r"^chr", "", r["chrom"])
            score = min(1000, int(round(r["score"] * 100)))
            ob.write(f"{r['chrom']}\t{r['interval_start']}\t{r['interval_end']}\tcen{num}\t{score}\t.\n")
            op.write(f"{r['chrom']}\t{r['peak_bp']}\t{r['peak_bp'] + 1}\tcen{num}\t{score}\t.\n")
    with open(a.out_prefix + ".telomeres.tsv", "w") as out:
        out.write("chrom\tend\ttelomere_repeat_bp\ttelomere_repeat_present\thic_peak_near_end\thic_telomere_peak\n")
        for r in tel_rows:
            out.write("\t".join(map(str, r)) + "\n")
    print("done:", ", ".join(f"{r['chrom']} {r.get('peak_mb')} Mb ({r['flag']})" for r in rows))


if __name__ == "__main__":
    main()
