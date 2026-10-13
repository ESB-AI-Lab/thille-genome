#!/usr/bin/env python3
"""Sum the per-chromosome samtools stats parts (run_stats_parallel.sh) into Table Z.

mapping rate = reads mapped / (raw total sequences over chromosomes + unplaced unmapped reads)
mean MAPQ    = mean of the samtools stats MAPQ histogram (mapped primary reads; no secondary/supplementary/dup/QC-fail)
error rate   = mismatches (NM: mismatches + inserted + deleted bases) / bases mapped (cigar)
QV           = -10 log10(error rate)
Usage (from the repo root): python3 summaries/mapping_stats/aggregate_stats.py > summaries/mapping_stats/table_z.tsv
"""
import collections, glob, math, re
KEYS = {"raw total sequences", "reads mapped", "reads MQ0", "supplementary alignments", "bases mapped (cigar)", "mismatches"}
names = ["gwen.clr.thille", "thille.hifi.thille", "gwen.clr.hass_uoq", "thille.hifi.hass_uoq", "gwen.clr.wi", "thille.hifi.wi", "westindian.hifi.thille", "westindian.hifi.hass_uoq", "westindian.hifi.wi"]
print("alignment\tchromosome_parts\ttotal_reads\treads_mapped\tmapping_rate_pct\tmean_mapq\tmapq60_pct\tmapq0_pct\tbases_mapped_cigar\tNM\terror_rate\tQV")
for n in names:
    sn, mq = collections.Counter(), collections.Counter()
    parts = sorted(glob.glob(f"summaries/mapping_stats/parts/{n}.*.stats.txt"))
    for p in parts:
        for line in open(p):
            f = line.rstrip("\n").split("\t")
            if f[0] == "SN" and f[1].rstrip(":") in KEYS:
                sn[f[1].rstrip(":")] += int(float(f[2]))
            elif f[0] == "MAPQ":
                mq[int(f[1])] += int(f[2])
    unplaced = int(open(f"summaries/mapping_stats/parts/{n}.unmapped_unplaced.count").read())
    total = sn["raw total sequences"] + unplaced
    nmq = sum(mq.values())
    err = sn["mismatches"] / sn["bases mapped (cigar)"]
    print(f"{n}\t{len(parts)}\t{total}\t{sn['reads mapped']}\t{100*sn['reads mapped']/total:.2f}\t"
          f"{sum(k*v for k,v in mq.items())/nmq:.2f}\t{100*sum(v for k,v in mq.items() if k>=60)/nmq:.2f}\t{100*mq[0]/nmq:.2f}\t"
          f"{sn['bases mapped (cigar)']}\t{sn['mismatches']}\t{err:.6f}\t{-10*math.log10(err):.2f}")
