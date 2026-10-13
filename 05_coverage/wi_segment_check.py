#!/usr/bin/env python3
"""Is each A/B coverage-hit segment present in the West Indian (WI) T2T assembly, and which copy is WI's closest to?

Segments are the deletion-hit runs (contiguous >= 1 kb in low or half; scripts/coverage_sv_hits.py) at the A/B panel
genes, on the reference they were called on (Thille or UoQ Hass). For each segment:
  1. the segment is aligned to WI with minimap2 2.30 -x asm5 -c; WI coverage = union of the query intervals of the
     primary alignments / segment length; WI identity = matching bases / alignment block length (PAF cols 10 / 11)
     over those primary alignments;
  2. the WI sequence spanned by the primary alignments on the best WI chromosome is aligned in the same way to the
     Thille and the UoQ Hass assemblies, giving its coverage and identity against each.
The low fraction of each run (positions in the low class) is reported, since only runs absent in one cultivar can
separate the two states in a single-haplotype assembly. Nothing here calls a flowering type.

Usage (from the repo root):
    python3 scripts/wi_segment_check.py        # writes summaries/wi_segment_check.tsv
"""

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coverage_class_summary import class_index  # noqa: E402
import plot_gene_stacks as pgs  # noqa: E402

MINIMAP2 = "/gpool/bin/minimap2-2.30_x64-linux/minimap2"
IMG = "containers/aligners_v3.sif"
GENOMES = {
    "thille": "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta",
    "hass": "genome/hass_GCA_029852735.1/hass_UoQ_chr_only.fasta",
    "wi": "genome/wi_GCA_051132755.1/wi.chr.fasta",
}
HITS = {"thille": ("summaries/coverage_sv_hits_thille.tsv", ["gwen.clr.thille", "thille.hifi.thille"]),
        "hass": ("summaries/coverage_sv_hits_hass.tsv", ["gwen.clr.hass_uoq", "thille.hifi.hass_uoq"])}
FAMILIES = {"SDMYB", "AP2", "NAC25", "MIK2", "SUD1", "Y5713", "OTU7", "APT3", "SPA1", "KCS11", "PAR1", "CCR1"}
MIK2_CONTIGS = {"chr12", "CM056820.1"}  # the chr12 candidate only
WORK = Path("summaries/wi_segment_work")


def sing(*cmd):
    return subprocess.run(["singularity", "exec", "--bind", "/epool", IMG, *cmd], check=True, text=True, capture_output=True).stdout


def index(name):
    mmi = WORK / f"{name}.asm5.mmi"
    if not mmi.exists():
        subprocess.run([MINIMAP2, "-x", "asm5", "-t", "16", "-d", str(mmi), GENOMES[name]], check=True, capture_output=True)
    return mmi


def align(fasta, target):
    """PAF primary alignments per query: {query: [(qlen, qs, qe, tname, ts, te, matches, block)]}"""
    out = subprocess.run([MINIMAP2, "-c", "-t", "16", str(index(target)), str(fasta)], check=True, text=True, capture_output=True).stdout
    hits = defaultdict(list)
    for line in out.splitlines():
        f = line.split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        hits[f[0]].append((int(f[1]), int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])))
    return hits


def summarize(alns, qlen):
    if not alns:
        return 0.0, 0.0, ""
    iv = sorted((a[1], a[2]) for a in alns)
    covered, cur_s, cur_e = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s > cur_e:
            covered += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    covered += cur_e - cur_s
    ident = sum(a[6] for a in alns) / sum(a[7] for a in alns)
    by_t = defaultdict(int)
    for a in alns:
        by_t[a[3]] += a[2] - a[1]
    return covered / qlen, ident, max(by_t, key=by_t.get)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    segments = {}  # (ref, contig, start, end) -> {genes, called_in}
    for ref, (path, tracks) in HITS.items():
        for r in csv.DictReader(open(path), delimiter="\t"):
            contig = r["locus"].split(":")[0]
            if r["gene"] not in FAMILIES or (r["gene"] == "MIK2" and contig not in MIK2_CONTIGS):
                continue
            for t in tracks:
                for run in filter(None, r[f"{t}_deletion_runs"].split(";")):
                    s, e = map(int, run.split("-"))
                    seg = segments.setdefault((ref, contig, s, e), {"genes": set(), "called_in": set()})
                    seg["genes"].add(r["gene"])
                    seg["called_in"].add("Gwen" if t.startswith("gwen") else "Thille")

    rows, fa = [], WORK / "segments.fa"
    with open(fa, "w") as out:
        for i, ((ref, contig, s, e), seg) in enumerate(sorted(segments.items())):
            name = f"seg{i}"
            seq = "".join(sing("samtools", "faidx", GENOMES[ref], f"{contig}:{s}-{e}").splitlines()[1:])
            out.write(f">{name}\n{seq}\n")
            low = {}
            for t in HITS[ref][1]:
                depth, meta = pgs.load_track("cache", t, contig)
                cls = class_index(np.asarray(depth[s - 1:e]), int(meta["median"]))
                low["Gwen" if t.startswith("gwen") else "Thille"] = (float(np.mean(cls == 0)), float(np.mean(cls == 1)))
            rows.append({"segment": name, "genes": ",".join(sorted(seg["genes"])), "called_on": ref,
                         "region": f"{contig}:{s}-{e}", "length_bp": e - s + 1,
                         "deletion_run_in": ",".join(sorted(seg["called_in"])),
                         "gwen_low_half_frac": "{:.2f}/{:.2f}".format(*low["Gwen"]),
                         "thille_low_half_frac": "{:.2f}/{:.2f}".format(*low["Thille"])})

    wi = align(fa, "wi")
    wi_fa = WORK / "wi_regions.fa"
    with open(wi_fa, "w") as out:
        for row in rows:
            alns = wi.get(row["segment"], [])
            cov, ident, best = summarize(alns, row["length_bp"])
            row.update({"wi_coverage": f"{cov:.3f}", "wi_identity": f"{ident:.4f}" if alns else "", "wi_region": ""})
            on_best = [a for a in alns if a[3] == best]
            if on_best:
                ts, te = min(a[4] for a in on_best), max(a[5] for a in on_best)
                row["wi_region"] = f"{best}:{ts + 1}-{te}"
                seq = "".join(sing("samtools", "faidx", GENOMES["wi"], row["wi_region"]).splitlines()[1:])
                out.write(f">{row['segment']}\n{seq}\n")
    for target in ("thille", "hass"):
        hits = align(wi_fa, target)
        for row in rows:
            if not row["wi_region"]:
                row[f"wi_region_vs_{target}_cov"] = row[f"wi_region_vs_{target}_id"] = ""
                continue
            s, e = row["wi_region"].split(":")[1].split("-")
            cov, ident, best = summarize(hits.get(row["segment"], []), int(e) - int(s) + 1)
            row[f"wi_region_vs_{target}_cov"] = f"{cov:.3f}"
            row[f"wi_region_vs_{target}_id"] = f"{ident:.4f}" if hits.get(row["segment"]) else ""

    with open("summaries/wi_segment_check.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} segments", file=sys.stderr)


if __name__ == "__main__":
    main()
