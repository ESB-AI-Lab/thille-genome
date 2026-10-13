#!/usr/bin/env python3
"""
genome_composition_donut.py

Publication-quality nested donut plot of genome composition, parsed from a
RepeatMasker .tbl summary. Nothing is hard coded; every value comes from the
input files.

  Inner ring : major repeat classes plus non-repetitive sequence (sums to 100%)
  Outer ring : subfamily breakdown for classes RepeatMasker subdivides
  Gene track : optional (--gff). Gene bodies overlap repeats (intronic TEs),
               so gene space is drawn as a separate outer track rather than a
               slice of the same partition.

Percentages use assembly length excluding N/X runs, matching RepeatMasker.

Copied into thille_paper/scripts from /epool/edwin/pamericana/avoflow/annotation/repeatmodeler/thille.hic/ on
2026-10-09 to add --no-text (Edwin: no explanatory text inside figures; the inner-ring percentages and the center
summary go to the legend md, the color key keeps its percentages). Figure 2 (v4) is reproduced by
  python3 scripts/genome_composition_donut.py <thille.hic>/repeatmask/thille.hifiv4.asm.hic.scaffold.chr.only.fasta.tbl \
    --gff <thille.hic>/thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf -o <prefix>
(same TSV as thille_genome_composition.tsv; the PNG differs only in text anti-aliasing).

Usage
  python genome_composition_donut.py thille.tbl -o avocado_composition
  python genome_composition_donut.py thille.tbl --gff genes.gff3 --panel b -o fig2b
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Patch, Wedge  # noqa: E402

# ---------------------------------------------------------------------------
# Class grouping (structure only; values always come from the .tbl)
# label, color (Okabe-Ito, colorblind safe), source node names, subdivide, "other" label
# ---------------------------------------------------------------------------
GROUPS = [
    ("LTR retrotransposons", "#0072B2", ["LTR elements"], True, "Other LTR"),
    ("Non-LTR retrotransposons", "#56B4E9", ["LINEs", "SINEs", "Penelope"], True, "Other non-LTR"),
    ("DNA transposons", "#009E73", ["DNA transposons"], True, "Other DNA"),
    ("Rolling-circle (Helitron)", "#E69F00", ["Rolling-circles"], False, "Other RC"),
    ("Unclassified repeats", "#CC79A7", ["Unclassified"], False, "Other"),
    ("Tandem and low complexity", "#D55E00",
     ["Satellites", "Simple repeats", "Low complexity", "Small RNA"], True, "Other tandem"),
]
NONREP_LABEL, NONREP_COLOR = "Non-repetitive", "#D9D9D9"
GENE_COLOR = "#117733"   # Paul Tol dark green, colorblind safe

SHORT_NAMES = {
    "ty1copia": "Copia", "gypsydirs1": "Gypsy", "belpao": "BEL/Pao",
    "l1cin4": "L1", "rtebovb": "RTE", "l2cr1rex": "L2/CR1/Rex",
    "hoboactivator": "hAT", "touristharbinger": "PIF/Harbinger",
    "tc1is630pogo": "Tc1/Mariner", "enspm": "CACTA", "mudris905": "MuDR",
}

ROW_RE = re.compile(
    r"^(?P<indent>[ \t]*)(?P<name>\S.*?):?\s+(?P<count>\d+)\s+(?P<bp>\d+)\s+bp\s+"
    r"(?P<pct>[\d.]+)\s*%\s*$"
)


@dataclass
class Node:
    name: str
    count: int
    bp: int
    level: int
    children: list["Node"] = field(default_factory=list)


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def indent_level(ws: str) -> int:
    # RepeatMasker indentation is inconsistent (e.g. 4 vs 5 spaces for LINE
    # subfamilies), so bucket widths into hierarchy levels.
    n = len(ws.expandtabs(8))
    return 0 if n == 0 else 1 if n <= 3 else 2 if n <= 6 else 3


def parse_tbl(path: Path) -> tuple[dict, list[Node]]:
    text = path.read_text()

    def grab(pattern, cast=str, required=False):
        m = re.search(pattern, text)
        if not m:
            if required:
                sys.exit(f"ERROR: could not parse /{pattern}/ from {path}")
            return None
        return tuple(cast(g) for g in m.groups()) if len(m.groups()) > 1 else cast(m.group(1))

    total_bp, total_noN = grab(r"total length:\s*(\d+)\s*bp\s*\((\d+)\s*bp excl", int, True)
    meta = {
        "file": grab(r"file name:\s*(\S+)") or path.name,
        "sequences": grab(r"sequences:\s*(\d+)", int),
        "total_bp": total_bp,
        "total_bp_noN": total_noN,
        "masked_bp": grab(r"bases masked:\s*(\d+)\s*bp", int),
        "gc": grab(r"GC level:\s*([\d.]+)", float),
        "rm_version": grab(r"RepeatMasker version\s*([\w.]+)"),
    }

    roots: list[Node] = []
    stack: list[Node] = []
    for line in text.splitlines():
        m = ROW_RE.match(line)
        if not m:
            continue
        node = Node(m["name"].strip().rstrip(":").strip(), int(m["count"]),
                    int(m["bp"]), indent_level(m["indent"]))
        while stack and stack[-1].level >= node.level:
            stack.pop()
        (stack[-1].children if stack else roots).append(node)
        stack.append(node)
    if not roots:
        sys.exit(f"ERROR: no repeat rows parsed from {path}")
    return meta, roots


def index_nodes(roots: list[Node]) -> tuple[dict[str, Node], dict[int, Node | None]]:
    idx: dict[str, Node] = {}
    parent: dict[int, Node | None] = {}

    def walk(nodes, par):
        for n in nodes:
            idx.setdefault(norm(n.name), n)
            parent[id(n)] = par
            walk(n.children, n)

    walk(roots, None)
    return idx, parent


def ancestors(n: Node, parent: dict[int, Node | None]):
    a = parent.get(id(n))
    while a is not None:
        yield a
        a = parent.get(id(a))


def build_partition(meta: dict, roots: list[Node], min_sub_frac: float) -> list[dict]:
    denom = meta["total_bp_noN"]
    idx, parent = index_nodes(roots)

    # resolve source nodes; a node may sit under another group's source
    # (e.g. Rolling-circles nested under DNA transposons in some RM versions)
    spec = []
    selected: set[int] = set()
    for label, color, sources, subdivide, other_label in GROUPS:
        nodes = [idx[norm(s)] for s in sources if norm(s) in idx]
        spec.append((label, color, nodes, subdivide, other_label))
        selected.update(id(n) for n in nodes)

    # effective bp: subtract any selected descendant so nothing is counted twice
    eff = {}
    for _, _, nodes, _, _ in spec:
        for n in nodes:
            eff[id(n)] = n.bp
    for _, _, nodes, _, _ in spec:
        for n in nodes:
            for a in ancestors(n, parent):
                if id(a) in selected:
                    eff[id(a)] -= n.bp
                    print(f"NOTE: '{n.name}' is nested under '{a.name}'; "
                          f"counted once, in its own class", file=sys.stderr)
                    break

    # top-level rows not covered by any group (future RM classes) go to "Other repeats"
    def covered(n: Node) -> bool:
        if id(n) in selected:
            return True
        return any(covered(k) for k in n.children)

    unassigned = [r for r in roots if r.bp > 0 and not covered(r)]
    for r in unassigned:
        print(f"NOTE: unmapped class '{r.name}' ({r.bp:,} bp) plotted as 'Other repeats'",
              file=sys.stderr)

    groups = []
    for label, color, nodes, subdivide, other_label in spec:
        bp = sum(eff[id(n)] for n in nodes)
        if bp <= 0:
            continue
        subs: list[tuple[str, int]] = []
        if subdivide:
            for n in nodes:
                kids = [k for k in n.children if k.bp > 0 and id(k) not in selected]
                if kids:
                    subs += [(SHORT_NAMES.get(norm(k.name), k.name), k.bp) for k in kids]
                    rest = eff[id(n)] - sum(k.bp for k in kids)
                    if rest > 0:
                        subs.append((other_label, rest))
                elif eff[id(n)] > 0:
                    subs.append((SHORT_NAMES.get(norm(n.name), n.name), eff[id(n)]))
            # fold tiny subfamilies into "Other" so the outer ring stays legible
            kept, folded = [], []
            for name, sbp in subs:
                if sbp / denom < min_sub_frac or name == other_label:
                    folded.append((name, sbp))
                else:
                    kept.append((name, sbp))
            kept.sort(key=lambda x: -x[1])
            if len(folded) == 1:
                kept.append(folded[0])
            elif folded:
                kept.append((other_label, sum(b for _, b in folded)))
            subs = kept if len(kept) > 1 else []
        groups.append({"label": label, "color": color, "bp": bp, "subs": subs})
    if unassigned:
        groups.append({"label": "Other repeats", "color": "#999999",
                       "bp": sum(r.bp for r in unassigned), "subs": []})

    class_sum = sum(g["bp"] for g in groups)
    masked = meta["masked_bp"] if meta["masked_bp"] is not None else class_sum
    diff = class_sum - masked
    if diff > 0:
        # RM class totals are per annotation; 'bases masked' is the union.
        # Overlapping calls (e.g. simple repeats inside TEs) are counted in both classes.
        print(f"NOTE: class totals exceed 'bases masked' by {diff:,} bp "
              f"({100 * diff / denom:.3f}% of assembly) from overlapping annotations. "
              f"Non-repetitive is taken from 'bases masked'; arcs are scaled to close the ring.",
              file=sys.stderr)
    elif diff < 0:
        print(f"WARNING: class totals fall short of 'bases masked' by {-diff:,} bp; "
              f"check the .tbl for rows the parser did not read", file=sys.stderr)
    groups.append({"label": NONREP_LABEL, "color": NONREP_COLOR,
                   "bp": denom - masked, "subs": []})
    ring_total = sum(g["bp"] for g in groups)
    for g in groups:
        g["pct"] = 100 * g["bp"] / denom
        g["arc"] = 360 * g["bp"] / ring_total
        g["subs"] = [(n, b, 100 * b / denom, 360 * b / ring_total) for n, b in g["subs"]]
    meta["overlap_bp"] = max(diff, 0)
    return groups


def gene_span(path: Path, feature: str) -> tuple[int, int]:
    """Union of feature intervals (bp) from GFF3/GTF, so overlapping models count once."""
    opener = gzip.open if path.suffix == ".gz" else open
    ivs = defaultdict(list)
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[2] == feature:
                ivs[f[0]].append((int(f[3]) - 1, int(f[4])))
    n = sum(len(v) for v in ivs.values())
    if n == 0:
        sys.exit(f"ERROR: no '{feature}' features in {path}")
    total = 0
    for v in ivs.values():
        v.sort()
        cs, ce = v[0]
        for a, b in v[1:]:
            if a > ce:
                total += ce - cs
                cs, ce = a, b
            else:
                ce = max(ce, b)
        total += ce - cs
    return total, n


# ---------------------------------------------------------------------------
# Exact base-level partition mode (RepeatMasker .out + GTF/GFF + FASTA)
# Priority, highest first: N > CDS > UTR > repeats (higher SW score wins) > intron > intergenic
# TEs inside introns and UTR-free gene space are therefore counted as TE, and
# intronic sequence without a TE is counted as intron. No base is counted twice.
# ---------------------------------------------------------------------------
GENIC_LABEL, GENIC_COLOR = "Genic, non-repeat", GENE_COLOR
INTERGENIC_LABEL = "Intergenic, non-repeat"
OTHER_REP = ("Other repeats", "#999999", "Other")
PART_ORDER = [(g[0], g[1], g[4]) for g in GROUPS] + [OTHER_REP] + \
    [(GENIC_LABEL, GENIC_COLOR, "Other genic")]


def open_any(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def classify(class_family: str) -> tuple[str, str]:
    cls, _, fam = class_family.partition("/")
    c, f = cls.rstrip("?").lower(), fam.rstrip("?").lower()
    if c == "ltr":
        sub = "Gypsy" if "gypsy" in f else "Copia" if "copia" in f else "Other LTR"
        return "LTR retrotransposons", sub
    if c in ("line", "sine", "ple", "retroposon"):
        sub = ("L1" if f.startswith("l1") else "RTE" if f.startswith("rte")
               else "SINE" if c == "sine" else "Penelope" if c == "ple" else "Other non-LTR")
        return "Non-LTR retrotransposons", sub
    if c == "dna":
        for key, name in (("hat", "hAT"), ("pif", "PIF/Harbinger"), ("harbinger", "PIF/Harbinger"),
                          ("cmc", "CACTA"), ("enspm", "CACTA"), ("mule", "MuDR"),
                          ("mudr", "MuDR"), ("tcmar", "Tc1/Mariner")):
            if f.startswith(key):
                return "DNA transposons", name
        return "DNA transposons", "Other DNA"
    if c == "rc":
        return "Rolling-circle (Helitron)", "Helitron"
    if c in ("unknown", "unspecified", "unclassified"):
        return "Unclassified repeats", "Unclassified"
    if c == "simple_repeat":
        return "Tandem and low complexity", "Simple repeats"
    if c == "low_complexity":
        return "Tandem and low complexity", "Low complexity"
    if c == "satellite":
        return "Tandem and low complexity", "Satellites"
    return OTHER_REP[0], cls


def load_rm_out(path: Path):
    by_seq = defaultdict(list)
    with open_any(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit():
                continue
            by_seq[f[4]].append((int(f[0]), int(f[5]) - 1, int(f[6]), classify(f[10])))
    if not by_seq:
        sys.exit(f"ERROR: no annotations parsed from {path}")
    return by_seq


def load_annotation(path: Path):
    """Return per-sequence lists of CDS, exon, and gene-span intervals (0-based, half-open)."""
    cds, exon, span = defaultdict(list), defaultdict(list), defaultdict(list)
    tx_bounds = {}
    have_span = False
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seq, ftype, a, b = f[0], f[2], int(f[3]) - 1, int(f[4])
            if ftype == "CDS":
                cds[seq].append((a, b))
            elif ftype == "exon":
                exon[seq].append((a, b))
            elif ftype in ("gene", "transcript", "mRNA"):
                span[seq].append((a, b))
                have_span = True
            if ftype in ("CDS", "exon"):
                m = re.search(r'transcript_id[ =]"?([^";]+)', f[8]) or re.search(r"Parent=([^;]+)", f[8])
                if m:
                    key = (seq, m.group(1))
                    lo, hi = tx_bounds.get(key, (a, b))
                    tx_bounds[key] = (min(lo, a), max(hi, b))
    if not have_span:  # derive transcript spans from CDS/exon extents
        for (seq, _), iv in tx_bounds.items():
            span[seq].append(iv)
    if not cds:
        sys.exit(f"ERROR: no CDS features in {path}")
    return cds, exon, span


def read_fasta(path: Path):
    name, chunks = None, []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).encode()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).encode()


def exact_partition(fasta: Path, rm_out: Path, annot: Path, min_sub_frac: float):
    import numpy as np

    reps = load_rm_out(rm_out)
    cds, exon, span = load_annotation(annot)

    # label codes
    N_, INTERGENIC, INTRON, UTR, CDS_ = 0, 1, 2, 3, 4
    rep_keys = sorted({k for v in reps.values() for *_, k in v})
    rep_code = {k: 10 + i for i, k in enumerate(rep_keys)}
    if len(rep_code) + 10 > 255:
        sys.exit("ERROR: too many repeat subclasses")

    counts = np.zeros(256, dtype=np.int64)
    total_bp = n_bp = gc_bp = nseq = 0
    rep_union = te_intron = te_utr = te_cds = gene_union = 0
    seen = set()
    for name, seq in read_fasta(fasta):
        nseq += 1
        seen.add(name)
        arr = np.frombuffer(seq, dtype=np.uint8) & 0xDF          # uppercase
        nmask = (arr == ord("N")) | (arr == ord("X"))
        total_bp += arr.size
        n_bp += int(nmask.sum())
        gc_bp += int(((arr == ord("G")) | (arr == ord("C"))).sum())

        L = np.full(arr.size, INTERGENIC, dtype=np.uint8)
        gs = np.zeros(arr.size, dtype=bool)
        for a, b in span.get(name, ()):
            L[a:b] = INTRON
            gs[a:b] = True
        rep = np.zeros(arr.size, dtype=bool)
        for _score, a, b, key in sorted(reps.get(name, ()), key=lambda r: r[0]):
            L[a:b] = rep_code[key]
            rep[a:b] = True
        exm = np.zeros(arr.size, dtype=bool)
        for a, b in exon.get(name, ()):
            L[a:b] = UTR
            exm[a:b] = True
        cm = np.zeros(arr.size, dtype=bool)
        for a, b in cds.get(name, ()):
            L[a:b] = CDS_
            cm[a:b] = True
        L[nmask] = N_

        ok = ~nmask
        rep_union += int((rep & ok).sum())
        gene_union += int((gs & ok).sum())
        te_cds += int((rep & cm & ok).sum())
        te_utr += int((rep & exm & ~cm & ok).sum())
        te_intron += int((rep & gs & ~exm & ~cm & ok).sum())
        counts += np.bincount(L, minlength=256)

    missing = (set(reps) | set(cds)) - seen
    if missing:
        print(f"NOTE: {len(missing)} sequence(s) in .out/annotation absent from FASTA "
              f"were ignored (e.g. {sorted(missing)[0]})", file=sys.stderr)

    denom = total_bp - n_bp
    meta = {"file": fasta.name, "sequences": nseq, "total_bp": total_bp,
            "total_bp_noN": denom, "gc": 100 * gc_bp / denom, "masked_bp": rep_union,
            "overlap_bp": 0, "te_intron": te_intron, "te_utr": te_utr, "te_cds": te_cds,
            "gene_union": gene_union}

    # assemble groups
    bucket = defaultdict(list)
    bucket[GENIC_LABEL] = [("CDS", int(counts[CDS_])), ("UTR", int(counts[UTR])),
                           ("Intron", int(counts[INTRON]))]
    for key, code in rep_code.items():
        bucket[key[0]].append((key[1], int(counts[code])))

    groups = []
    for label, color, other_label in PART_ORDER:
        subs = [(n, b) for n, b in bucket.get(label, []) if b > 0]
        bp = sum(b for _, b in subs)
        if bp == 0:
            continue
        kept, folded = [], []
        for n, b in subs:
            (folded if (b / denom < min_sub_frac and label != GENIC_LABEL) or n == other_label
             else kept).append((n, b))
        if label != GENIC_LABEL:
            kept.sort(key=lambda x: -x[1])
        if len(folded) == 1:
            kept.append(folded[0])
        elif folded:
            kept.append((other_label, sum(b for _, b in folded)))
        groups.append({"label": label, "color": color, "bp": bp,
                       "subs": kept if len(kept) > 1 else []})
    groups.append({"label": INTERGENIC_LABEL, "color": NONREP_COLOR,
                   "bp": int(counts[INTERGENIC]), "subs": []})
    assert sum(g["bp"] for g in groups) == denom, "partition does not sum to assembly"
    for g in groups:
        g["pct"] = 100 * g["bp"] / denom
        g["arc"] = 360 * g["bp"] / denom
        g["subs"] = [(n, b, 100 * b / denom, 360 * b / denom) for n, b in g["subs"]]
    meta["center"] = [f"{100 * rep_union / denom:.1f}% repetitive",
                      f"{100 * gene_union / denom:.1f}% gene space"]
    return groups, meta



def tint(color: str, t: float) -> tuple:
    r, g, b = to_rgb(color)
    return (r + (1 - r) * t, g + (1 - g) * t, b + (1 - b) * t)


def text_color(rgb) -> str:
    r, g, b = to_rgb(rgb)
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.55 else "#1A1A1A"


def polar(r: float, deg: float) -> tuple[float, float]:
    import math
    a = math.radians(deg)
    return r * math.cos(a), r * math.sin(a)


def plot(groups, meta, args, gene=None):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(args.width, args.height))
    ax = fig.add_axes([0.0, 0.0, 0.5, 1.0])
    ax.set_aspect("equal")
    lim = 1.42
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis("off")

    R_IN, W_IN = 1.0, 0.30      # inner ring
    R_OUT, W_OUT = 1.20, 0.17   # outer ring (starts just outside inner ring)
    R_GENE, W_GENE = 1.30, 0.045
    gap = dict(edgecolor="white", linewidth=0.9)

    angle = 90.0  # 12 o'clock, clockwise
    nonrep_start, deg_per_pct = None, None
    for g in groups:
        span = g["arc"]
        ax.add_patch(Wedge((0, 0), R_IN, angle - span, angle, width=W_IN,
                           facecolor=g["color"], **gap))
        mid = angle - span / 2
        if g["pct"] >= args.min_label and not args.no_text:
            x, y = polar(R_IN - W_IN / 2, mid)
            col = text_color(g["color"])
            if g["label"] == GENIC_LABEL:
                ax.text(x, y, f"Genes\n{g['pct']:.1f}%", ha="center", va="center",
                        fontsize=7, fontweight="bold", linespacing=1.3, color=col)
            else:
                ax.text(x, y, f"{g['pct']:.1f}%", ha="center", va="center", fontsize=7,
                        fontweight="bold", color=col)
        if g["label"] == NONREP_LABEL:
            nonrep_start, deg_per_pct = angle, g["arc"] / g["pct"]
        if g["subs"]:
            sub_angle = angle
            n = len(g["subs"])
            for i, (name, _bp, pct, sspan) in enumerate(g["subs"]):
                shade = tint(g["color"], 0.15 + 0.55 * i / max(n - 1, 1))
                ax.add_patch(Wedge((0, 0), R_OUT, sub_angle - sspan, sub_angle,
                                   width=W_OUT, facecolor=shade, **gap))
                arc_len = (R_OUT - W_OUT / 2) * sspan * 3.14159 / 180
                if pct >= args.min_sub_label and 0.045 * len(name) < 0.85 * arc_len:
                    mid_s = sub_angle - sspan / 2
                    x, y = polar(R_OUT - W_OUT / 2, mid_s)
                    rot = (mid_s - 90) % 360          # tangent to the ring
                    if 90 < rot < 270:
                        rot -= 180                    # keep text upright
                    ax.text(x, y, name, ha="center", va="center", fontsize=5.5,
                            rotation=rot, rotation_mode="anchor", color=text_color(shade))
                sub_angle -= sspan
        angle -= span

    if gene is not None:
        # genes are non-repetitive sequence, so the track sits on the outer ring
        # directly over the non-repetitive wedge rather than floating at 12 o'clock
        start = nonrep_start if nonrep_start is not None else 90.0
        gspan = gene["pct"] * (deg_per_pct if deg_per_pct else 3.6)
        if nonrep_start is not None and gene["pct"] > groups[-1]["pct"]:
            print(f"NOTE: gene space ({gene['pct']:.2f}%) exceeds non-repetitive "
                  f"({groups[-1]['pct']:.2f}%); intronic and UTR repeats put part of gene "
                  f"space inside the repeat wedges. Use --rm-out for an exact partition.",
                  file=sys.stderr)
        ax.add_patch(Wedge((0, 0), R_OUT, start - gspan, start, width=W_OUT,
                           facecolor=args.gene_color, **gap))
        mid_g = start - gspan / 2
        x, y = polar(R_OUT - W_OUT / 2, mid_g)
        rot = (mid_g - 90) % 360
        if 90 < rot < 270:
            rot -= 180
        ax.text(x, y, "Genes", ha="center", va="center", fontsize=6,
                fontweight="bold", rotation=rot, rotation_mode="anchor",
                color=text_color(args.gene_color))

    # center summary
    if args.no_text:
        lines, meta = [], dict(meta, gc=None, sequences=None)
    else:
        lines = meta.get("center") or [
        f"{100 * (meta['total_bp_noN'] - groups[-1]['bp']) / meta['total_bp_noN']:.1f}% repetitive"]
    if not args.no_text:
        ax.text(0, 0.22, f"{meta['total_bp'] / 1e6:,.1f} Mb", ha="center", va="center",
                fontsize=11, fontweight="bold", color="#1A1A1A")
    for j, txt in enumerate(lines):
        ax.text(0, 0.04 - 0.13 * j, txt, ha="center", va="center",
                fontsize=7.5, color="#333333")
    sub = []
    if meta["gc"] is not None:
        sub.append(f"GC {meta['gc']:.1f}%")
    if meta["sequences"]:
        sub.append(f"{meta['sequences']:,} scaffolds")
    if sub:
        ax.text(0, 0.04 - 0.13 * len(lines) - 0.03, "  |  ".join(sub), ha="center", va="center",
                fontsize=6, color="#666666")

    # hierarchical legend
    handles = []
    for g in groups:
        handles.append(Patch(facecolor=g["color"], edgecolor="none",
                             label=f"{g['label']}   {g['pct']:.2f}%"))
        n = len(g["subs"])
        for i, (name, _bp, pct, _arc) in enumerate(g["subs"]):
            handles.append(Patch(facecolor=tint(g["color"], 0.15 + 0.55 * i / max(n - 1, 1)),
                                 edgecolor="none", label=f"    {name}   {pct:.2f}%"))
    if gene is not None:
        handles.append(Patch(facecolor=args.gene_color, edgecolor="none",
                             label=f"Gene space (outer track)   {gene['pct']:.2f}%"))
    leg = fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.52, 0.5),
                     frameon=False, fontsize=6.5, handlelength=1.0, handleheight=1.0,
                     handletextpad=0.6, labelspacing=0.45,
                     title="% of assembly (excl. N)", title_fontsize=7)
    leg._legend_box.align = "left"
    leg.get_title().set_fontweight("bold")

    if args.title:
        fig.suptitle(args.title, fontsize=9, fontweight="bold", y=0.98)
    if args.panel:
        fig.text(0.01, 0.97, f"({args.panel})", fontsize=10, fontweight="bold",
                 va="top", ha="left")

    out = Path(args.out)
    for ext in args.formats:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=args.dpi, bbox_inches="tight",
                    pad_inches=0.03, transparent=args.transparent)
        print(f"wrote {out.with_suffix('.' + ext)}")
    plt.close(fig)


def write_table(groups, meta, out: Path, gene=None):
    path = out.with_suffix(".tsv")
    with open(path, "w") as fh:
        fh.write("class\tsubclass\tbp\tpercent_of_assembly_exclN\n")
        for g in groups:
            fh.write(f"{g['label']}\t\t{g['bp']}\t{g['pct']:.4f}\n")
            for name, bp, pct, _arc in g["subs"]:
                fh.write(f"{g['label']}\t{name}\t{bp}\t{pct:.4f}\n")
        if meta.get("overlap_bp"):
            fh.write(f"Overlapping repeat annotations (counted in 2 classes)\t\t"
                     f"{meta['overlap_bp']}\t{100 * meta['overlap_bp'] / meta['total_bp_noN']:.4f}\n")
        for key, lab in (("te_intron", "Repeat bp inside introns (in repeat classes)"),
                         ("te_utr", "Repeat bp inside UTRs (assigned to UTR)"),
                         ("te_cds", "Repeat bp inside CDS (assigned to CDS)"),
                         ("gene_union", "Gene space incl. repeats (union)")):
            if key in meta:
                fh.write(f"{lab}\t\t{meta[key]}\t{100 * meta[key] / meta['total_bp_noN']:.4f}\n")
        if gene is not None:
            fh.write(f"Gene space (overlaps repeats)\t\t{gene['bp']}\t{gene['pct']:.4f}\n")
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tbl", type=Path, nargs="?", help="RepeatMasker .tbl summary (tbl mode)")
    p.add_argument("--rm-out", type=Path, help="RepeatMasker .out (exact partition mode)")
    p.add_argument("--fasta", type=Path, help="unmasked assembly FASTA (exact partition mode)")
    p.add_argument("-o", "--out", default="genome_composition", help="output prefix")
    p.add_argument("--gff", type=Path, help="GFF3/GTF (optionally .gz) for gene space track")
    p.add_argument("--gene-color", default=GENE_COLOR,
                   help="color of the gene space track")
    p.add_argument("--gene-feature", default="gene", help="feature type for gene space")
    p.add_argument("--title", default=None)
    p.add_argument("--panel", default=None, help="panel letter, e.g. b")
    p.add_argument("--min-sub", type=float, default=0.5,
                   help="fold subfamilies below this %% of assembly into 'Other'")
    p.add_argument("--min-label", type=float, default=3.0,
                   help="minimum %% for an in-ring label")
    p.add_argument("--min-sub-label", type=float, default=5.0,
                   help="minimum %% for an outer-ring subfamily label")
    p.add_argument("--width", type=float, default=7.2, help="inches (7.2 = double column)")
    p.add_argument("--height", type=float, default=3.4, help="inches")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--formats", nargs="+", default=["pdf", "svg", "png"])
    p.add_argument("--transparent", action="store_true")
    p.add_argument("--no-text", action="store_true", help="no center summary and no inner-ring percentages "
                   "(default off; output without it is unchanged)")
    args = p.parse_args()

    if args.rm_out:
        if not (args.fasta and args.gff):
            sys.exit("ERROR: exact partition mode needs --rm-out, --fasta and --gff")
        groups, meta = exact_partition(args.fasta, args.rm_out, args.gff, args.min_sub / 100)
        print(f"repeats inside introns: {meta['te_intron']:,} bp; inside UTRs: "
              f"{meta['te_utr']:,} bp; inside CDS: {meta['te_cds']:,} bp", file=sys.stderr)
        if args.tbl:
            tmeta, _ = parse_tbl(args.tbl)
            print(f"check: repeat union {meta['masked_bp']:,} bp vs .tbl bases masked "
                  f"{tmeta['masked_bp']:,} bp", file=sys.stderr)
        plot(groups, meta, args)
        write_table(groups, meta, Path(args.out))
        return
    if not args.tbl:
        sys.exit("ERROR: give a .tbl, or --rm-out with --fasta and --gff")

    meta, roots = parse_tbl(args.tbl)
    groups = build_partition(meta, roots, args.min_sub / 100)

    gene = None
    if args.gff:
        gbp, n = gene_span(args.gff, args.gene_feature)
        gene = {"bp": gbp, "pct": 100 * gbp / meta["total_bp_noN"]}
        print(f"{n:,} '{args.gene_feature}' features, {gbp:,} bp union", file=sys.stderr)

    plot(groups, meta, args, gene)
    write_table(groups, meta, Path(args.out), gene)


if __name__ == "__main__":
    main()
