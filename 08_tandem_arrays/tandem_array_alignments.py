#!/usr/bin/env python3
"""Pairwise alignment of the Tiberius gene models in each Thille tandem group.

For every group in summaries/tandem_copies.tsv (Thille), the gene models overlapping its copies are taken as the
array members. For every pair of members three sequences are aligned globally:
  region  - the whole gene span (genome, gene orientation)          EMBOSS stretcher, DNA
  exon    - the exons of the longest transcript, in transcript order EMBOSS stretcher, DNA
  protein - the Tiberius protein of that transcript                  EMBOSS needle, BLOSUM62
A group is "supported" when at least one pair of members has exon identity >= 70 % over aligned columns
(Edwin, 2026-09-26). Nothing is demoted here; the distances are reported for Edwin's decision.

Tools run inside containers/aligners.sif (EMBOSS 6.6.0, samtools 1.21).

Usage (from the repo root):
    python3 scripts/tandem_array_alignments.py [--jobs 16]
Outputs: summaries/tandem_alignment_pairs.tsv, summaries/tandem_array_evidence.tsv,
         summaries/tandem_alignments/<group>/{region,exon,protein}.fa and raw/ (EMBOSS output).
"""

import argparse
import collections
import csv
import re
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import plot_gene_stacks as p  # noqa: E402

GENOME = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta"
PROTEINS = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.proteins.fasta"
GTF = "annotation/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.gtf"
COPIES = "summaries/tandem_copies.tsv"
IMAGE = "containers/aligners.sif"
OUTDIR = Path("summaries/tandem_alignments")
IDENTITY_SUPPORT = 0.70
MIN_GENE_BP = 1000
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def log(message):
    print(message, file=sys.stderr, flush=True)


def run(args):
    """Run a command inside the aligners image; raise with its stderr on failure."""
    result = subprocess.run(["singularity", "exec", IMAGE, *args], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(args[:2])} failed (exit {result.returncode}): {result.stderr.strip()[-500:]}")
    return result.stdout


def read_fasta_list(text):
    """FASTA text -> [(name, sequence)] in file order; keeps duplicate names (EMBOSS writes both sequences)."""
    records, name, parts = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].split()[0], []
        elif line.strip():
            parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def read_fasta(text):
    records, name, parts = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        elif line.strip():
            parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq):
    return "".join(CODONS.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def fetch(regions):
    """samtools faidx for a list of 'contig:start-end' (1-based inclusive) -> {region: uppercase sequence}."""
    records = {}
    for start in range(0, len(regions), 500):  # keep the command line short
        records.update(read_fasta(run(["samtools", "faidx", GENOME, *regions[start:start + 500]])))
    return {k: v.upper() for k, v in records.items()}


def model_sequences(models):
    """gene model -> {'region', 'exon'} in gene orientation."""
    regions = []
    for m in models:
        regions.append(f"{m.contig}:{m.start}-{m.stop}")
        regions += [f"{m.contig}:{s}-{e}" for s, e in m.exons]
    seqs = fetch(sorted(set(regions)))
    out = {}
    for m in models:
        region = seqs[f"{m.contig}:{m.start}-{m.stop}"]
        exon = "".join(seqs[f"{m.contig}:{s}-{e}"] for s, e in sorted(m.exons))
        if m.strand == "-":
            region, exon = reverse_complement(region), reverse_complement(exon)
        expected = sum(e - s + 1 for s, e in m.exons)
        if len(exon) != expected or len(region) != m.stop - m.start + 1:
            raise ValueError(f"{m.gene_id}: extracted lengths {len(region)}/{len(exon)} != annotation {m.stop - m.start + 1}/{expected}")
        out[m.gene_id] = {"region": region, "exon": exon}
    return out


def alignment_metrics(a, b):
    """Metrics from two gapped aligned sequences of equal length."""
    if len(a) != len(b):
        raise ValueError("aligned sequences differ in length")
    both = [i for i in range(len(a)) if a[i] != "-" and b[i] != "-"]
    identical = sum(1 for i in both if a[i] == b[i])
    len_a, len_b = len(a) - a.count("-"), len(b) - b.count("-")
    return {
        "len_a": len_a, "len_b": len_b, "aln_len": len(a),
        "identity_global": round(identical / len(a), 4) if len(a) else 0.0,
        "identity_aligned": round(identical / len(both), 4) if both else 0.0,
        "coverage_short": round(len(both) / min(len_a, len_b), 4) if min(len_a, len_b) else 0.0,
        "gap_columns": len(a) - len(both),
    }


def align(kind, fa_a, fa_b, raw_path):
    if kind == "protein":
        cmd = ["needle", "-asequence", fa_a, "-bsequence", fa_b, "-gapopen", "10", "-gapextend", "0.5"]
    else:
        cmd = ["stretcher", "-asequence", fa_a, "-bsequence", fa_b, "-snucleotide1", "-snucleotide2"]
    run(cmd + ["-aformat", "fasta", "-outfile", str(raw_path)])
    aligned = [seq for _name, seq in read_fasta_list(raw_path.read_text())]
    if len(aligned) != 2:
        raise RuntimeError(f"{raw_path}: expected 2 aligned sequences, got {len(aligned)}")
    return alignment_metrics(aligned[0].upper(), aligned[1].upper())


def write_distinct_arrays(evidence, pairs, by_id, rank):
    """Merge groups that share gene models into distinct physical arrays; one row per array and per member."""
    multi = [e for e in evidence if e["n_models"] >= 2]
    parent = {e["group"]: e["group"] for e in multi}

    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x
    owner = {}
    for e in multi:
        for m in e["models"].split(","):
            if m in owner:
                parent[root(e["group"])] = root(owner[m])
            else:
                owner[m] = e["group"]
    arrays = collections.defaultdict(list)
    for e in multi:
        arrays[root(e["group"])].append(e)
    lookup = {(r["model_a"], r["model_b"], r["type"]): r for r in pairs}

    def pair(a, b, kind):
        a, b = sorted((a, b), key=lambda i: (by_id[i].contig, by_id[i].start))
        return lookup.get((a, b, kind))

    array_rows, member_rows = [], []
    for groups in sorted(arrays.values(), key=lambda es: (es[0]["contig"], min(int(e["models_span"].split("-")[0]) for e in es))):
        members = sorted({m for e in groups for m in e["models"].split(",")}, key=lambda i: by_id[i].start)
        supported = sorted({m for e in groups for m in e["models_in_ge_70pct_pairs"].split(",") if m},
                           key=lambda i: by_id[i].start)
        contig = by_id[members[0]].contig
        name = f"{contig}:{by_id[members[0]].start}-{max(by_id[m].stop for m in members)}"
        adjacent = [(a, b, by_id[b].start - by_id[a].stop - 1, rank[b] - rank[a] - 1) for a, b in zip(supported, supported[1:])]
        array_rows.append({
            "array": name, "groups": ",".join(e["group"] for e in groups),
            "status": "supported" if supported else "not supported",
            "n_members": len(members), "n_members_in_ge_70pct_pairs": len(supported),
            "n_supported_members_ge_1kb": sum(by_id[m].stop - by_id[m].start + 1 >= MIN_GENE_BP for m in supported),
            "strands": "".join(by_id[m].strand for m in members),
            "adjacent_supported_bp": ",".join(str(x[2]) for x in adjacent),
            "adjacent_supported_genes_between": ",".join(str(x[3]) for x in adjacent),
            "max_adjacent_bp": max((x[2] for x in adjacent), default=""),
            "max_adjacent_genes_between": max((x[3] for x in adjacent), default=""),
        })
        for a in members:
            partners = [(pair(a, b, "exon")["identity_aligned"], pair(a, b, "exon")["coverage_short"], b)
                        for b in members if b != a and pair(a, b, "exon")]
            best = max(partners) if partners else None
            m = by_id[a]
            member_rows.append({
                "array": name, "model": a, "transcript": m.transcript_id, "contig": contig, "start": m.start,
                "end": m.stop, "strand": m.strand, "gene_bp": m.stop - m.start + 1,
                "cds_bp": sum(e - s + 1 for s, e in m.cds), "exons": len(m.exons),
                "best_exon_partner": best[2] if best else "", "exon_identity": best[0] if best else "",
                "exon_coverage": best[1] if best else "",
                "region_identity": pair(a, best[2], "region")["identity_aligned"] if best else "",
                "protein_identity": pair(a, best[2], "protein")["identity_aligned"] if best else "",
                "in_ge_70pct_pair": a in supported, "gene_ge_1kb": m.stop - m.start + 1 >= MIN_GENE_BP})
    for path, rows in (("summaries/tandem_arrays.tsv", array_rows), ("summaries/tandem_array_members.tsv", member_rows)):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    log(f"[arrays] {len(array_rows)} distinct arrays from {len(multi)} groups with >= 2 models")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=16, help="parallel alignments (default: %(default)s)")
    args = parser.parse_args()

    copies = [c for c in csv.DictReader(open(COPIES), delimiter="\t") if c["assembly"] == "Thille"]
    contigs = {c["contig"] for c in copies}
    models_by_contig = p.read_gene_models(GTF, contigs)
    by_id = {m.gene_id: m for v in models_by_contig.values() for m in v}
    rank = {m.gene_id: i for v in models_by_contig.values() for i, m in enumerate(v)}  # genome order per contig
    members = collections.defaultdict(set)
    for c in copies:
        members[c["group"]].update(re.findall(r"([^,(]+)\(", c["gene_models"]))
    groups = sorted({c["group"] for c in copies})

    needed = sorted({i for g in groups for i in members[g]}, key=lambda i: (by_id[i].contig, by_id[i].start))
    log(f"[extract] {len(needed)} gene models from {GENOME}")
    seqs = model_sequences([by_id[i] for i in needed])
    proteins = read_fasta(open(PROTEINS).read())
    for gene_id in needed:
        m = by_id[gene_id]
        protein = proteins.get(m.transcript_id)
        if protein is None:
            raise ValueError(f"{m.transcript_id} not in {PROTEINS}")
        translated = translate(seqs[gene_id]["exon"])
        if translated.rstrip("*") != protein.rstrip("*"):
            raise ValueError(f"{gene_id}: translated CDS does not match protein {m.transcript_id}")
        seqs[gene_id]["protein"] = protein.rstrip("*")
    log(f"[extract] all {len(needed)} exon sequences translate to their Tiberius protein")

    jobs = []
    for group in groups:
        ids = sorted(members[group], key=lambda i: (by_id[i].contig, by_id[i].start))
        if len(ids) < 2:
            continue
        folder = OUTDIR / p.safe_file_name(group.split("|")[0] if group.split("|")[0] == group.split("|")[1] else group)
        (folder / "raw").mkdir(parents=True, exist_ok=True)
        for kind in ("region", "exon", "protein"):
            with open(folder / f"{kind}.fa", "w") as handle:
                for i in ids:
                    handle.write(f">{i} {by_id[i].contig}:{by_id[i].start}-{by_id[i].stop}({by_id[i].strand})\n{seqs[i][kind]}\n")
            for i in ids:
                (folder / "raw" / f"{kind}_{i}.fa").write_text(f">{i}\n{seqs[i][kind]}\n")
        for x in range(len(ids)):
            for y in range(x + 1, len(ids)):
                for kind in ("region", "exon", "protein"):
                    jobs.append((group, ids[x], ids[y], kind, folder))

    log(f"[align] {len(jobs)} alignments with {args.jobs} workers")

    def work(job):
        group, a, b, kind, folder = job
        raw = folder / "raw" / f"{kind}_{a}_vs_{b}.aln.fa"
        metrics = align(kind, str(folder / "raw" / f"{kind}_{a}.fa"), str(folder / "raw" / f"{kind}_{b}.fa"), raw)
        ma, mb = by_id[a], by_id[b]
        return {"group": group, "model_a": a, "model_b": b, "type": kind, **metrics,
                "bp_between": mb.start - ma.stop - 1, "genes_between": rank[b] - rank[a] - 1,
                "gene_bp_a": ma.stop - ma.start + 1, "gene_bp_b": mb.stop - mb.start + 1}

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pairs = list(pool.map(work, jobs))

    with open("summaries/tandem_alignment_pairs.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pairs[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(pairs)

    evidence = []
    for group in groups:
        ids = sorted(members[group], key=lambda i: (by_id[i].contig, by_id[i].start))
        exon_pairs = [r for r in pairs if r["group"] == group and r["type"] == "exon"]
        best = lambda kind: max((r["identity_aligned"] for r in pairs if r["group"] == group and r["type"] == kind), default="")
        support = [r for r in exon_pairs if r["identity_aligned"] >= IDENTITY_SUPPORT]
        supported_models = sorted({r["model_a"] for r in support} | {r["model_b"] for r in support},
                                  key=lambda i: by_id[i].start)
        if len(ids) < 2:
            status = "no gene-model array (0 models)" if not ids else "no gene-model array (1 model)"
        else:
            status = "supported" if support else "not supported"
        top = max(exon_pairs, key=lambda r: r["identity_aligned"], default=None)
        evidence.append({
            "group": group, "contig": by_id[ids[0]].contig if ids else "", "n_models": len(ids),
            "n_models_ge_1kb": sum(by_id[i].stop - by_id[i].start + 1 >= MIN_GENE_BP for i in ids),
            "models": ",".join(ids),
            "models_span": f"{min(by_id[i].start for i in ids)}-{max(by_id[i].stop for i in ids)}" if ids else "",
            "best_exon_identity": best("exon"), "best_exon_pair": f"{top['model_a']}/{top['model_b']}" if top else "",
            "best_exon_pair_coverage": top["coverage_short"] if top else "",
            "median_exon_identity": round(statistics.median(r["identity_aligned"] for r in exon_pairs), 4) if exon_pairs else "",
            "best_region_identity": best("region"), "best_protein_identity": best("protein"),
            "pairs_ge_70pct_exon": len(support), "models_in_ge_70pct_pairs": ",".join(supported_models),
            "max_genes_between_in_supporting_pairs": max((r["genes_between"] for r in support), default=""),
            "max_bp_between_in_supporting_pairs": max((r["bp_between"] for r in support), default=""),
            "proposed_status": status,
        })
    with open("summaries/tandem_array_evidence.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(evidence)
    write_distinct_arrays(evidence, pairs, by_id, rank)
    counts = collections.Counter(e["proposed_status"] for e in evidence)
    log(f"[done] {len(pairs)} alignments; groups: {dict(counts)}")


if __name__ == "__main__":
    main()
