"""Shared helpers for the Fst flowering-gene allele analysis (scripts/gene_alleles*.py, plot_gene_matrix.py).

Coordinates are Thille chromosome positions (1-based). Exons are CDS exons in transcript order, from
summaries/fst_alleles/gene_models.tsv (scripts/gene_models.py).
"""
import csv
import re
import subprocess
from pathlib import Path

IMG = "containers/aligners_v3.sif"
REF = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta"
MODELS = Path("summaries/fst_alleles/gene_models.tsv")
OUT = Path("summaries/fst_alleles")
COMP = str.maketrans("ACGTN", "TGCAN")
CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def sh(*cmd):
    return subprocess.run(["singularity", "exec", "--no-home", "--bind", "/epool/edwin/pamericana/avoflow/mappingback:/epool/edwin/pamericana/avoflow/mappingback:ro", IMG, *cmd], check=True, text=True,
                          capture_output=True).stdout


class GeneModel:
    def __init__(self, row):
        self.__dict__.update(row)
        for k in ("gene_start", "gene_end", "win_start", "win_end", "cds_bp"):
            setattr(self, k, int(row[k]))
        self.exons = [] if row["cds_exons"] == "none" else [tuple(map(int, x.split("-"))) for x in row["cds_exons"].split(",")]
        # CDS positions in transcript order
        self.cds_pos = [p for s, e in self.exons for p in (range(e, s - 1, -1) if self.strand == "-" else range(s, e + 1))]

    def orient(self, s):
        return s.translate(COMP)[::-1] if self.strand == "-" else s

    def base_on_gene_strand(self, b):
        return b.translate(COMP) if self.strand == "-" else b

    def region_of(self, pos):
        """5prime_flank, E1, I1, E2, ..., En, 3prime_flank (transcript numbering); 'gene' when there is no CDS model."""
        if not self.exons:
            if self.gene_start <= pos <= self.gene_end:
                return "gene"
        else:
            for i, (s, e) in enumerate(self.exons):
                if s <= pos <= e:
                    return f"E{i + 1}"
            for i in range(len(self.exons) - 1):
                a, b = self.exons[i], self.exons[i + 1]
                lo, hi = (b[1], a[0]) if self.strand == "-" else (a[1], b[0])
                if lo < pos < hi:
                    return f"I{i + 1}"
        upstream = pos > self.gene_end if self.strand == "-" else pos < self.gene_start
        return "5prime_flank" if upstream else "3prime_flank"

    def regions(self):
        if not self.exons:
            return ["5prime_flank", "gene", "3prime_flank"]
        r = ["5prime_flank"]
        for i in range(len(self.exons)):
            r.append(f"E{i + 1}")
            if i < len(self.exons) - 1:
                r.append(f"I{i + 1}")
        return r + ["3prime_flank"]

    def region_span(self, reg):
        """(low, high) Thille coordinates of a region inside the window."""
        if reg == "gene":
            return self.gene_start, self.gene_end
        if reg[0] == "E":
            return self.exons[int(reg[1:]) - 1]
        if reg[0] == "I":
            i = int(reg[1:]) - 1
            a, b = self.exons[i], self.exons[i + 1]
            return (b[1] + 1, a[0] - 1) if self.strand == "-" else (a[1] + 1, b[0] - 1)
        five = reg == "5prime_flank"
        if (self.strand == "-") == five:
            return self.gene_end + 1, self.win_end
        return self.win_start, self.gene_start - 1


def load_models():
    return {r["label"]: GeneModel(r) for r in csv.DictReader(open(MODELS), delimiter="\t")}


def ref_window(g):
    return "".join(sh("samtools", "faidx", REF, f"{g.chrom}:{g.win_start}-{g.win_end}").splitlines()[1:]).upper()


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def cs_target_bases(paf, w0):
    """From a minimap2 PAF with --cs=long (target = Thille window starting at w0): Thille position -> aligned base
    (target-plus-strand orientation, '-' for a deletion) and position -> inserted sequence after it. Primary only."""
    base, ins, blocks = {}, {}, []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
        t = w0 + int(f[7])
        blocks.append((int(f[7]), int(f[8]), int(f[2]), int(f[3])))
        for op, val in re.findall(r"([=:*+\-~])([A-Za-z0-9]*)", cs):
            if op == "=":
                for ch in val:
                    base[t] = ch.upper(); t += 1
            elif op == "*":
                base[t] = val[1].upper(); t += 1
            elif op == "-":
                for _ in val:
                    base[t] = "-"; t += 1
            elif op == "+":
                ins[t - 1] = ins.get(t - 1, "") + val.upper()
    return base, ins, sorted(blocks)
