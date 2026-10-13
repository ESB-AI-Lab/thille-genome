#!/usr/bin/env python3
"""Whole-genome read-mapping statistics from a minimap2 BAM (Edwin, 2026-10-02), reading every record once.

reads            primary records (mapped or unmapped; secondary/supplementary excluded) = reads in the FASTQ
mapping rate     primary mapped / reads
MAPQ             mean and % at 60 over primary mapped reads
error rates (each also as QV = -10 log10 rate):
  nm_all         NM summed over all non-secondary alignments (primary + supplementary) / their aligned bases
                 (M/=/X/I), which is what samtools stats reports as "error rate"
  nm_primary     NM / aligned bases, primary alignments only
  subst_primary  (NM - inserted - deleted bases) / aligned bases, primary only (substitutions only)
  gapcomp        minimap2 gap-compressed divergence pooled over primary alignments: sum(de * columns) / sum(columns)
median read identity = 1 - median per-read de (primary)
Usage: samtools view -@ 2 <bam> | python3 scripts/mapping_stats_full.py <label> > <label>.tsv
"""
import sys, re, math, statistics
lab = sys.argv[1]
OPS = re.compile(r"(\d+)([MIDNSHP=X])")
reads = mapped = supp = mq_sum = mq60 = 0
al_all = nm_all = al_p = nm_p = ins_p = del_p = 0
de_w = de_cols = 0.0; des = []
for line in sys.stdin:
    f = line.rstrip("\n").split("\t"); flag = int(f[1])
    if flag & 256: continue
    is_supp = bool(flag & 2048)
    if not is_supp:
        reads += 1
        if flag & 4: continue
        mapped += 1; mq = int(f[4]); mq_sum += mq; mq60 += mq >= 60
    else:
        supp += 1
    ops = OPS.findall(f[5]); al = ins = dl = 0; cols = 0
    for n, o in ops:
        n = int(n)
        if o in "M=X": al += n; cols += n
        elif o == "I": al += n; ins += n; cols += 1
        elif o == "D": dl += n; cols += 1
    tags = "\t".join(f[11:])   # all optional fields (NM, de, ... are separate columns)
    nm = int(re.search(r"NM:i:(\d+)", tags).group(1)); al_all += al; nm_all += nm
    if not is_supp:
        al_p += al; nm_p += nm; ins_p += ins; del_p += dl
        m = re.search(r"de:f:([0-9.eE-]+)", tags)
        if m:
            d = float(m.group(1)); des.append(d); de_w += d * cols; de_cols += cols
qv = lambda e: -10 * math.log10(e) if e > 0 else float("inf")
r = {"label": lab, "reads": reads, "primary_mapped": mapped, "mapping_rate_pct": 100 * mapped / reads,
     "mean_mapq": mq_sum / mapped, "mapq60_pct": 100 * mq60 / mapped, "supplementary": supp,
     "supp_per_mapped_read": supp / mapped,
     "err_nm_all": nm_all / al_all, "err_nm_primary": nm_p / al_p, "err_subst_primary": (nm_p - ins_p - del_p) / al_p,
     "err_gapcomp_primary": de_w / de_cols, "median_read_identity_pct": 100 * (1 - statistics.median(des))}
for k in ("err_nm_all", "err_nm_primary", "err_subst_primary", "err_gapcomp_primary"):
    r["qv_" + k[4:]] = qv(r[k])
print("\t".join(r)); print("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()))
