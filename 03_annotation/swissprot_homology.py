#!/usr/bin/env python3
"""Name Thille gene models by Swiss-Prot homology and flag predicted flowering function.

Queries are the Tiberius proteins of: every gene model in a Thille tandem group (summaries/tandem_array_members.tsv
and the single-model groups in summaries/tandem_copies.tsv) and the Thille model of every curated flowering gene
(flowering_beds/genes.tsv). They are searched with blastp against UniProtKB/Swiss-Prot (db/swissprot, release in
db/swissprot/reldate.txt) inside containers/aligners_v2.sif.

A hit qualifies (Edwin, 2026-09-27) when E <= 1e-10, identity >= 30 %, and the best HSP covers >= 50 % of both
the query and the subject. The best qualifying hit (highest bitscore) names the model "<gene>-like".
Predicted flowering function: the hit carries UniProt keyword KW-0287 (Flowering) or a GO biological process
equal to or below GO:0009908, GO:0048573, GO:0010048, GO:0010228 or GO:0009909 (is_a / part_of, go-basic.obo).

Usage (from the repo root):
    python3 scripts/swissprot_homology.py [--threads 16]
Output: summaries/gene_model_functions.tsv, db/swissprot/cc_function_subset.tsv
"""

import argparse
import csv
import gzip
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

IMAGE = "containers/aligners_v2.sif"
SPROT_FASTA_GZ = Path("db/swissprot/uniprot_sprot.fasta.gz")
SPROT_DB = Path("db/swissprot/blast/uniprot_sprot")
SPROT_META = Path("db/swissprot/swissprot_meta.tsv.gz")
GO_OBO = Path("db/go/go-basic.obo")
CC_FUNCTION = Path("db/swissprot/cc_function_subset.tsv")
PROTEINS = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.proteins.fasta"
MAX_EVALUE, MIN_IDENTITY, MIN_COVERAGE = 1e-10, 30.0, 0.50
FLOWERING_KEYWORD = "KW-0287"
FLOWERING_GO_ROOTS = ["GO:0009908", "GO:0048573", "GO:0010048", "GO:0010228", "GO:0009909"]
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"


def log(message):
    print(message, file=sys.stderr, flush=True)


def run(args):
    result = subprocess.run(["singularity", "exec", IMAGE, *args], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{args[0]} failed (exit {result.returncode}): {result.stderr.strip()[-800:]}")
    return result.stdout


def read_fasta(path):
    records, name, parts = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def go_descendants(roots):
    """All GO terms equal to or below the roots through is_a and part_of (go-basic.obo)."""
    children, term = {}, None
    for line in open(GO_OBO):
        line = line.rstrip("\n")
        if line == "[Term]":
            term = None
        elif line.startswith("id: GO:"):
            term = line[4:]
        elif term and line.startswith("is_a: "):
            children.setdefault(line.split()[1], set()).add(term)
        elif term and line.startswith("relationship: part_of "):
            children.setdefault(line.split()[2], set()).add(term)
    seen, stack = set(roots), list(roots)
    while stack:
        for child in children.get(stack.pop(), ()):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


def query_ids():
    """Tiberius transcript IDs to search, with where each came from."""
    sources = {}
    for row in csv.DictReader(open("summaries/tandem_array_members.tsv"), delimiter="\t"):
        sources.setdefault(row["transcript"], set()).add(f"tandem array {row['array']}")
    evidence = {r["group"]: r for r in csv.DictReader(open("summaries/tandem_array_evidence.tsv"), delimiter="\t")}
    for group, row in evidence.items():
        if row["n_models"] == "1":
            sources.setdefault(f"{row['models']}_t1", set()).add(f"single-model tandem group {group}")
    for row in csv.DictReader(open("flowering_beds/genes.tsv"), delimiter="\t"):
        model = row["thille_model"]
        if re.fullmatch(r"g\d+", model):  # skips "unannotated" and "absent (closest paralog ...)"
            sources.setdefault(f"{model}_t1", set()).add(f"flowering gene {row['label']} ({row['category']})")
    return sources


def fetch_cc_function(accessions):
    """UniProt 'Function [CC]' text for the given accessions, cached in db/swissprot/cc_function_subset.tsv."""
    cached = {}
    if CC_FUNCTION.is_file():
        cached = {r["accession"]: r["cc_function"] for r in csv.DictReader(open(CC_FUNCTION), delimiter="\t")}
    todo = sorted(set(accessions) - set(cached))
    for start in range(0, len(todo), 100):
        chunk = todo[start:start + 100]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": query, "fields": "accession,cc_function", "format": "tsv", "size": 500}))
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    lines = response.read().decode().splitlines()
                break
            except OSError as error:
                log(f"[uniprot] attempt {attempt + 1} failed: {error}")
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError("UniProt REST request failed 3 times")
        for line in lines[1:]:
            accession, _, text = line.partition("\t")
            cached[accession] = text
        for accession in chunk:
            cached.setdefault(accession, "")
    with open(CC_FUNCTION, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["accession", "cc_function"])
        writer.writerows(sorted(cached.items()))
    return cached


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()

    if not Path(str(SPROT_DB) + ".pin").is_file():
        SPROT_DB.parent.mkdir(parents=True, exist_ok=True)
        plain = SPROT_DB.parent / "uniprot_sprot.fasta"
        with gzip.open(SPROT_FASTA_GZ, "rt") as src, open(plain, "w") as dst:
            dst.write(src.read())
        run(["makeblastdb", "-in", str(plain), "-dbtype", "prot", "-out", str(SPROT_DB), "-parse_seqids"])
        log(f"[blastdb] built {SPROT_DB}")

    sources = query_ids()
    proteins = read_fasta(PROTEINS)
    missing = sorted(t for t in sources if t not in proteins)
    if missing:
        raise ValueError(f"transcripts not in {PROTEINS}: {missing}")
    query_path = Path("summaries/swissprot_queries.fa")
    with open(query_path, "w") as handle:
        for transcript in sorted(sources):
            handle.write(f">{transcript}\n{proteins[transcript].rstrip('*')}\n")
    log(f"[blastp] {len(sources)} proteins vs Swiss-Prot")
    out = run(["blastp", "-query", str(query_path), "-db", str(SPROT_DB), "-evalue", "1e-5",
               "-max_target_seqs", "25", "-outfmt", OUTFMT, "-num_threads", str(args.threads)])

    best_hsp = {}  # (query, subject) -> best HSP by bitscore
    for line in out.splitlines():
        f = line.split("\t")
        hsp = {"query": f[0], "subject": f[1], "pident": float(f[2]), "qstart": int(f[6]), "qend": int(f[7]),
               "sstart": int(f[8]), "send": int(f[9]), "evalue": float(f[10]), "bitscore": float(f[11]),
               "qlen": int(f[12]), "slen": int(f[13])}
        hsp["qcov"] = (hsp["qend"] - hsp["qstart"] + 1) / hsp["qlen"]
        hsp["scov"] = (abs(hsp["send"] - hsp["sstart"]) + 1) / hsp["slen"]
        key = (hsp["query"], hsp["subject"])
        if key not in best_hsp or hsp["bitscore"] > best_hsp[key]["bitscore"]:
            best_hsp[key] = hsp
    best = {}
    for hsp in best_hsp.values():
        ok = (hsp["evalue"] <= MAX_EVALUE and hsp["pident"] >= MIN_IDENTITY
              and hsp["qcov"] >= MIN_COVERAGE and hsp["scov"] >= MIN_COVERAGE)
        if ok and (hsp["query"] not in best or hsp["bitscore"] > best[hsp["query"]]["bitscore"]):
            best[hsp["query"]] = hsp

    def accession(subject):  # sp|P47927|AP2_ARATH -> P47927 ; also plain accession
        parts = subject.split("|")
        return parts[1] if len(parts) >= 2 else subject.split(".")[0]

    wanted = {accession(h["subject"]) for h in best.values()}
    meta = {}
    with gzip.open(SPROT_META, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Entry"] in wanted:
                meta[row["Entry"]] = row
    flowering_go = go_descendants(FLOWERING_GO_ROOTS)
    cc = fetch_cc_function(wanted)

    rows = []
    for transcript in sorted(sources, key=lambda t: int(re.sub(r"\D", "", t.split("_")[0]) or 0)):
        hit = best.get(transcript)
        row = {"model": transcript.split("_")[0], "transcript": transcript, "sources": "; ".join(sorted(sources[transcript])),
               "protein_aa": len(proteins[transcript].rstrip("*"))}
        if hit is None:
            row.update({"swissprot_accession": "", "swissprot_entry": "", "swissprot_gene": "", "proposed_name": "",
                        "protein_name": "", "organism": "", "identity": "", "qcov": "", "scov": "", "evalue": "",
                        "bitscore": "", "keywords": "", "go_bp": "", "flowering": "no qualifying hit",
                        "flowering_reason": "", "function": ""})
            rows.append(row)
            continue
        acc = accession(hit["subject"])
        m = meta.get(acc, {})
        gene = m.get("Gene Names (primary)", "") or m.get("Entry Name", "").split("_")[0]
        go_ids = set(re.findall(r"GO:\d{7}", m.get("Gene Ontology (biological process)", "")))
        reasons = []
        if FLOWERING_KEYWORD in m.get("Keyword ID", ""):
            reasons.append("KW-0287 Flowering")
        flowering_terms = sorted(go_ids & flowering_go)
        if flowering_terms:
            reasons.append("GO " + ",".join(flowering_terms))
        row.update({"swissprot_accession": acc, "swissprot_entry": m.get("Entry Name", ""), "swissprot_gene": gene,
                    "proposed_name": f"{gene}-like" if gene else "", "protein_name": m.get("Protein names", ""),
                    "organism": m.get("Organism", ""), "identity": round(hit["pident"], 2),
                    "qcov": round(hit["qcov"], 3), "scov": round(hit["scov"], 3), "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"], "keywords": m.get("Keywords", ""),
                    "go_bp": m.get("Gene Ontology (biological process)", ""),
                    "flowering": "yes" if reasons else "no", "flowering_reason": "; ".join(reasons),
                    "function": cc.get(acc, "").removeprefix("FUNCTION: ")})
        rows.append(row)
    with open("summaries/gene_model_functions.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    named = sum(1 for r in rows if r["swissprot_gene"])
    flowering = sum(1 for r in rows if r["flowering"] == "yes")
    log(f"[done] {len(rows)} models: {named} named, {len(rows) - named} without a qualifying hit, {flowering} flowering")


if __name__ == "__main__":
    main()
