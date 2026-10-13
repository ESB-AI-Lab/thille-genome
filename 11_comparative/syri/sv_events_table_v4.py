#!/usr/bin/env python3
# v4 copy (2026-10-05): v4 SyRI, Gwen SV, annotation and centromere inputs; outputs under compare/v4/.
# 2026-10-07: Puebla rows and Puebla coordinate columns added.
"""One spreadsheet-friendly table of every SV >= 100 kb between Thille and Hass UoQ, West Indian T2T, Puebla and Gwen.

Rows come from the SV BEDs written by sv_to_bed.py (compare/sv_bed/):
  - thille.sv.ge100kb.bed: every event with a Thille interval (Hass, West Indian, Puebla and Gwen rows)
  - hass/westindian/puebla.sv.ge100kb.bed: events with no Thille interval (sequence present only in that
    genome, SyRI NOTAL on the query side), added so the table holds all events
Columns give the coordinates in each genome side by side (empty where the genome has no interval), plus, for rows
with a Thille interval: Tiberius genes overlapped, the fraction covered by RepeatMasker TEs (union of intervals) and
the distance to the Thille centromere call (0 = overlaps it). Coordinates are BED-style (0-based start, end).

Usage (from assemblies/thille.hic): sv_events_table.py <out.tsv>
"""
import bisect
import collections
import sys

SV = "compare/v4/sv_bed"
GTF = "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf"
TE = "thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.TE.gff"
CEN = "centromere/v4/thille_centromeres.satellitefirst.bed"
COLS = ["event_id", "compared_genome", "type", "size_bp", "thille_chrom", "thille_start", "thille_end",
        "hass_chrom", "hass_start", "hass_end", "westindian_chrom", "westindian_start", "westindian_end",
        "puebla_chrom", "puebla_start", "puebla_end",
        "n_thille_genes", "thille_genes", "thille_te_fraction", "dist_to_thille_centromere_bp", "detail"]


def read_bed(path):
    for line in open(path):
        if not line.startswith("#"):
            yield line.rstrip("\n").split("\t")


def merged_intervals(path, chrom_col, s_col, e_col, feature=None):
    iv = collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) <= e_col or (feature and f[2] != feature):
            continue
        iv[f[chrom_col]].append((int(f[s_col]) - 1, int(f[e_col])))
    out = {}
    for c, v in iv.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = m
    return out


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    genes = collections.defaultdict(list)
    for line in open(GTF):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "gene":
            genes[f[0]].append((int(f[3]) - 1, int(f[4]), f[8].split('"')[1]))
    te = merged_intervals(TE, 0, 3, 4)
    te_starts = {c: [x[0] for x in v] for c, v in te.items()}
    cen = {f[0]: (int(f[1]), int(f[2])) for f in (l.split() for l in open(CEN))}

    def te_frac(c, s, e):
        iv = te.get(c, [])
        i = max(0, bisect.bisect_right(te_starts.get(c, []), s) - 1)
        t = 0
        while i < len(iv) and iv[i][0] < e:
            t += max(0, min(e, iv[i][1]) - max(s, iv[i][0]))
            i += 1
        return t / (e - s) if e > s else 0.0

    rows = []
    for f in read_bed(f"{SV}/thille.sv.ge100kb.bed"):
        c, s, e = f[0], int(f[1]), int(f[2])
        src, typ, size = f[8], f[6], f[7]
        partner = {"hass": ("", "", ""), "westindian": ("", "", ""), "puebla": ("", "", "")}
        if src in partner and f[9] != ".":
            partner[src] = (f[9], f[10], f[11])
        g = [n for gs, ge, n in genes[c] if gs < e and ge > s]
        cs, ce = cen[c]
        dist = 0 if s < ce and e > cs else min(abs(s - ce), abs(cs - e))
        rows.append([f[3], src, typ, size, c, s, e, *partner["hass"], *partner["westindian"], *partner["puebla"], len(g),
                     ",".join(g), f"{te_frac(c, s, e):.3f}", dist, f[12]])
    for src in ("hass", "westindian", "puebla"):
        for f in read_bed(f"{SV}/{src}.sv.ge100kb.bed"):
            if f[9] != ".":
                continue                      # has a Thille interval, already in the table
            own = (f[0], f[1], f[2])
            cols = [own if g == src else ("", "", "") for g in ("hass", "westindian", "puebla")]
            rows.append([f"{src}_{f[3].split('_', 1)[1]}", src, f[6], f[7], "", "", "", *(x for c in cols for x in c),
                         "", "", "", "", f[12]])
    order = {f"chr{i}": i for i in range(1, 13)}
    rows.sort(key=lambda r: (r[4] == "", order.get(r[4], 99), r[5] if r[5] != "" else 0, r[1]))
    with open(sys.argv[1], "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    n = collections.Counter((r[1], r[4] != "") for r in rows)
    print(f"{len(rows)} rows -> {sys.argv[1]}")
    for (src, has_t), k in sorted(n.items()):
        print(f"  {src:10s} {'with Thille interval' if has_t else 'query-only'}: {k}")


if __name__ == "__main__":
    main()
