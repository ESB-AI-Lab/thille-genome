#!/usr/bin/env python3
"""Raw-read pileups and breakpoint support for the inversion windows (compare/v4/inversions/windows.tsv).

For every window and every read set aligned to that window's genome:
  mpileup  : samtools mpileup -aa -B -Q 0 -q 0 -s -f <ref> (-aa: every position, also where no read aligns) (no mapping- or base-quality filter, so poorly mapped reads
             stay visible; -s adds the MAPQ of each read) -> mpileup/<inv>_<genome>_<side>.<readset>.mpileup.gz
  bins     : per BIN bp: mean depth of all reads, mean depth of reads with MAPQ >= MQ, and the fraction of read bases
             that differ from the reference (substitutions, insertions and deletions per covering read base)
             -> bins/<same name>.bins.tsv
  support  : for each inversion breakpoint inside the window (with the other breakpoint of the same inversion as
             partner), counted over primary alignments with MAPQ >= MQ:
               spanning   reads aligned continuously over the breakpoint with >= SPAN bp on both sides
               clipped    reads whose alignment ends within +/- CLIPWIN bp of the breakpoint with a clip >= CLIP bp
               split      clipped reads with a supplementary alignment (SA tag) within +/- PARTNER bp of the partner
                          breakpoint, on the same chromosome and the opposite strand (the inversion signature)
             -> support/task<N>.tsv; clip and split positions for plotting -> clips/<same name>.clips.tsv
Chromosome median depth (all reads) per read set and chromosome, from 30 random 10 kb windows (seed 1), for the
depth / median scale -> chrom_medians/<readset>__<genome>.tsv (computed once, reused).

Usage (inside the asmbenchmarkers container, from assemblies/thille.hic):
  inv_pileup.py --task N --ntasks M [--select genome:readset ...] [--tag _x] [--bin 500] [--mq 20] [--span 5000] [--clip 500] [--clipwin 1000] [--partner 10000]
"""
import argparse
import collections
import gzip
import os
import random
import re
import subprocess
import sys

import pysam

D = "compare/v4/inversions"
SAM = "/gpool/bin/samtools-1.24/samtools"
MB = "/epool/edwin/pamericana/avoflow/mappingback"
REF = {"thille": "thille.hifiv4.asm.hic.scaffold.chr.only.fasta",
       "hass": "hass_UoQ_chr_only.fasta",
       "westindian": "westindian/ncbi/ncbi_dataset/data/GCA_051132755.1/GCA_051132755.1_XTBG_Pame_1.0_genomic.fna",
       "puebla": "/gpool/edwin/avocados/thille_paper/genome/puebla/GCA_060308345.1/ncbi_dataset/data/GCA_060308345.1/"
                 "GCA_060308345.1_ASM6030834v1_genomic.fna"}
BAMS = {"thille": [("thille_hifi", f"{MB}/v4/thille.hifi.thille_v4.bam"),
                   ("lowland_hifi", f"{MB}/v4/westindian.hifi.thille_v4.bam"),
                   ("gwen_clr", f"{MB}/v4/gwen.clr.thille_v4.bam"),
                   ("puebla_hifi", f"{MB}/puebla/minimap2/puebla.hifi.thille_v4.bam")],
        "hass": [("thille_hifi", f"{MB}/minimap2/thille.hifi.hass_uoq.bam"),
                 ("lowland_hifi", f"{MB}/minimap2/westindian.hifi.hass_uoq.bam"),
                 ("gwen_clr", f"{MB}/gwen.clr.hass_uoq.bam")],
        "westindian": [("thille_hifi", f"{MB}/minimap2/thille.hifi.wi.bam"),
                       ("lowland_hifi", f"{MB}/minimap2/westindian.hifi.wi.bam"),
                       ("gwen_clr", f"{MB}/gwen.clr.wi.bam")],
        "puebla": [("thille_hifi", f"{MB}/puebla/minimap2/thille.hifi.puebla.bam"),
                   ("lowland_hifi", f"{MB}/puebla/minimap2/westindian.hifi.puebla.bam"),
                   ("gwen_clr", f"{MB}/puebla/gwen.clr.puebla.bam"),
                   ("puebla_hifi", f"{MB}/puebla/minimap2/puebla.hifi.puebla.bam")]}
INDEL = re.compile(r"[+-](\d+)")


def strip_indels(bases):
    """Remove +N<seq>/-N<seq> runs (their letters are not substitutions); return (rest, number of indels)."""
    out, i, k = [], 0, 0
    while i < len(bases):
        m = INDEL.match(bases, i)
        if m:
            i = m.end() + int(m.group(1))
            k += 1
        else:
            out.append(bases[i])
            i += 1
    return "".join(out), k


def chrom_median(rs, bam, genome, chrom):
    path = f"{D}/chrom_medians/{rs}__{genome}.tsv"
    med = {}
    if os.path.exists(path):
        med = dict((l.split("\t")[0], float(l.split("\t")[1])) for l in open(path))
    if chrom in med:
        return med[chrom]
    L = dict(zip(pysam.AlignmentFile(bam).references, pysam.AlignmentFile(bam).lengths))[chrom]
    rng = random.Random(1)
    vals = []
    for _ in range(30):
        s = rng.randrange(0, L - 10000)
        out = subprocess.run([SAM, "depth", "-a", "-r", f"{chrom}:{s + 1}-{s + 10000}", bam],
                             capture_output=True, text=True, check=True).stdout
        vals += [int(l.split("\t")[2]) for l in out.splitlines()]
    vals.sort()
    m = float(vals[len(vals) // 2])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + f".tmp{os.getpid()}", "w") as fh:  # per-process temp name; parallel tasks may race
        for k, v in {**med, chrom: m}.items():
            fh.write(f"{k}\t{v}\n")
    os.replace(path + f".tmp{os.getpid()}", path)
    return m


def pileup(bam, ref, chrom, s, e, out_gz, bins_tsv, binsize, mq):
    cmd = [SAM, "mpileup", "-aa", "-B", "-Q", "0", "-q", "0", "-s", "-f", ref, "-r", f"{chrom}:{s + 1}-{e}", bam]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    nb = (e - s) // binsize
    dep, depq, mism, cov = [0] * nb, [0] * nb, [0] * nb, [0] * nb
    n = 0
    with gzip.open(out_gz + ".tmp", "wt") as gz:
        for line in p.stdout:
            gz.write(line)
            f = line.rstrip("\n").split("\t")
            b = (int(f[1]) - 1 - s) // binsize
            if not 0 <= b < nb:
                continue
            n += 1
            d = int(f[3])
            dep[b] += d
            depq[b] += sum(1 for c in f[6] if ord(c) - 33 >= mq) if d else 0
            if d:
                bases, indel = strip_indels(re.sub(r"\^.", "", f[4]))
                sub = sum(1 for c in bases if c in "ACGTNacgtn")
                dels = bases.count("*") + bases.count("#")
                mism[b] += sub + dels + indel
                cov[b] += d
    if p.wait() != 0:
        sys.exit(f"samtools mpileup failed for {chrom}:{s + 1}-{e} {bam}")
    os.replace(out_gz + ".tmp", out_gz)
    with open(bins_tsv, "w") as fh:
        fh.write("bin_start\tbin_end\tmean_depth\tmean_depth_mq\tmismatch_frac\n")
        for i in range(nb):
            fh.write(f"{s + i * binsize}\t{s + (i + 1) * binsize}\t{dep[i] / binsize:.2f}\t{depq[i] / binsize:.2f}\t"
                     f"{(mism[i] / cov[i]) if cov[i] else 0:.4f}\n")
    return n


def support(bam, chrom, bp, partner, a):
    """Spanning, clipped and split-read counts at breakpoint bp (0-based position between bases)."""
    span = clip = split = 0
    marks = []
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(chrom, max(0, bp - a.partner), bp + a.partner):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < a.mq:
                continue
            if r.reference_start <= bp - a.span and r.reference_end >= bp + a.span:
                span += 1
                continue
            ct = r.cigartuples
            ends = []
            if ct[0][0] in (4, 5) and ct[0][1] >= a.clip:
                ends.append(r.reference_start)
            if ct[-1][0] in (4, 5) and ct[-1][1] >= a.clip:
                ends.append(r.reference_end)
            for x in ends:
                if abs(x - bp) > a.clipwin:
                    continue
                clip += 1
                is_split = False
                if r.has_tag("SA"):
                    for sa in r.get_tag("SA").rstrip(";").split(";"):
                        c, pos, strand = sa.split(",")[:3]
                        if (c == chrom and abs(int(pos) - 1 - partner) <= a.partner
                                and (strand == "-") != r.is_reverse):
                            is_split = True
                split += is_split
                marks.append((x, "split" if is_split else "clip"))
    return span, clip, split, marks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--ntasks", type=int, required=True)
    ap.add_argument("--bin", type=int, default=500)
    ap.add_argument("--mq", type=int, default=20)
    ap.add_argument("--span", type=int, default=5000)
    ap.add_argument("--clip", type=int, default=500)
    ap.add_argument("--clipwin", type=int, default=1000)
    ap.add_argument("--partner", type=int, default=10000)
    ap.add_argument("--select", nargs="*", help="only these genome:readset pairs ('*' allowed for readset), "
                    "e.g. puebla:* thille:puebla_hifi; default all")
    ap.add_argument("--tag", default="", help="suffix of the support file (support/task<N><tag>.tsv)")
    a = ap.parse_args()
    for sub in ("mpileup", "bins", "support", "clips", "chrom_medians"):
        os.makedirs(f"{D}/{sub}", exist_ok=True)
    rows = [l.rstrip("\n").split("\t") for l in list(open(f"{D}/windows.tsv"))[1:]]
    mine = [r for i, r in enumerate(rows) if i % a.ntasks == a.task]
    out = []
    for inv, pair, genome, side, chrom, ws, we, s, e in mine:
        ws, we, s, e = int(ws), int(we), int(s), int(e)
        name = f"{inv}_{genome}_{side}"
        for rs, bam in BAMS[genome]:
            if a.select and not any(x in (f"{genome}:{rs}", f"{genome}:*") for x in a.select):
                continue
            n = pileup(bam, REF[genome], chrom, ws, we, f"{D}/mpileup/{name}.{rs}.mpileup.gz",
                       f"{D}/bins/{name}.{rs}.bins.tsv", a.bin, a.mq)
            if n != we - ws:
                sys.exit(f"{name} {rs}: {n} pileup positions, expected {we - ws}")
            med = chrom_median(rs, bam, genome, chrom)
            marks = []
            for bp, partner, which in ((s, e, "start"), (e, s, "end")):
                if not ws <= bp < we:
                    continue
                sp, cl, sl, mk = support(bam, chrom, bp, partner, a)
                marks += mk
                out.append([inv, pair, genome, side, chrom, which, bp, rs, med, sp, cl, sl])
            with open(f"{D}/clips/{name}.{rs}.clips.tsv", "w") as fh:
                fh.write("pos\tkind\n")
                for x, k in sorted(marks):
                    fh.write(f"{x}\t{k}\n")
        print(f"done {name}", flush=True)
    with open(f"{D}/support/task{a.task}{a.tag}.tsv", "w") as fh:
        fh.write("inv_id\tpair\tgenome\tside\tchrom\tbreakpoint\tpos\treadset\tchrom_median_depth\tspanning\tclipped\t"
                 "split_to_partner\n")
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
