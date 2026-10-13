#!/usr/bin/env python3
"""Place the A/B panel genes on the West Indian (WI) T2T genome (GCA_051132755.1).

Same transfer as flowering_beds (README §Method, step 4): the Thille gene with its 1 kb upstream and downstream
flanks are aligned independently to WI with minimap2 2.30 (-a, as in the original transfer); the three primary
alignments are reported with chromosome, strand, MAPQ and the gaps between them, and a locus is marked
"tiled" when all three are on one chromosome and strand, in order ("contiguous" when both gaps are 0 bp, otherwise "gapped" with the gaps listed).
The WI gene model is the overlapping model whose protein has the highest blastp bitscore against the query protein;
the query protein is the chosen Thille model protein, or for KCS11 (no Thille ortholog) the UoQ Hass proteins.
Every query protein is also searched against all WI proteins, as an independent check of the model choice.
Accepted protein hits follow the Swiss-Prot rule used in this study (E <= 1e-10, >= 30 % identity, >= 50 % coverage
of query and subject).

Usage (from the repo root):
    python3 scripts/place_on_wi.py            # writes summaries/wi_placement.tsv and flowering_beds/igv/flowering_genes_wi.bed
    python3 scripts/place_on_wi.py --genes-file fst_genes.txt --label fst   # FAMILY|ver1 keys; *_fst outputs
"""

import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

WI_FASTA = "genome/wi_GCA_051132755.1/wi.chr.fasta"
WI_DIR = "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/westindian/ncbi/ncbi_dataset/data/GCA_051132755.1"
WI_GFF, WI_PROT = f"{WI_DIR}/genomic.gff", f"{WI_DIR}/protein.faa"
THILLE_FASTA = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.fasta"
THILLE_PROT = "genome/thille.hifiv3.asm.hic.scaffold.chr.only.masked.annotation.tiberius.proteins.fasta"
HASS_PROT = "genome/hass_GCA_029852735.1/GCA_029852735.1_ASM2985273v1_protein.faa.gz"
THILLE_BED, HASS_BED = "flowering_beds/igv/flowering_genes_thille.bed", "flowering_beds/igv/flowering_genes_hass.bed"
MINIMAP2 = "/gpool/bin/minimap2-2.30_x64-linux/minimap2"  # the build used for every alignment in this study
IMG = "containers/aligners_v3.sif"
GENES = ["SDMYB", "AP2", "NAC25", "MIK2", "SUD1", "Y5713", "OTU7", "APT3", "SPA1", "KCS11", "PAR1", "CCR1"]
MIK2_CHROM = "chr12"  # the chr12 candidate, not the chr5 annex gene
FLANK = 1_000
WORK = Path("summaries/wi_placement_work")
E_MAX, ID_MIN, COV_MIN = 1e-10, 30.0, 0.5


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kw).stdout


def sing(*cmd):
    return run(["singularity", "exec", "--no-home", "--bind", "/epool/edwin/pamericana/avoflow/assemblies/thille.hic/westindian/ncbi/ncbi_dataset/data/GCA_051132755.1:/epool/edwin/pamericana/avoflow/assemblies/thille.hic/westindian/ncbi/ncbi_dataset/data/GCA_051132755.1:ro", IMG, *cmd])


def read_fasta(text):
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def thille_rows(keys=None):
    rows, seen = [], set()
    for line in open(THILLE_BED):
        if line.startswith("track"):
            continue
        f = line.rstrip("\n").split("\t")
        fam, ver1, model = f[3].split("|")[:3]
        if keys is not None:
            if f"{fam}|{ver1}" not in keys or (f[0], f[1], f[2]) in seen:
                continue
        elif fam not in GENES or (fam == "MIK2" and f[0] != MIK2_CHROM) or (f[0], f[1], f[2]) in seen:
            continue
        seen.add((f[0], f[1], f[2]))  # NAC25 and KCS11 each have two Gwen genes on one Thille locus
        rows.append({"family": fam, "ver1": ver1, "model": model, "chrom": f[0], "start": int(f[1]), "end": int(f[2]), "strand": f[5]})
    return rows


def wi_models():
    """gene locus_tag -> (chrom, start, end, strand, protein_id)"""
    genes, prot = {}, {}
    for line in open(WI_GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        tag = re.search(r"locus_tag=([^;]+)", f[8])
        if f[2] == "gene" and tag:
            genes[tag.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] == "CDS" and tag:
            pid = re.search(r"protein_id=([^;]+)", f[8])
            if pid:
                prot.setdefault(tag.group(1), pid.group(1))
    return {t: (*g, prot.get(t, "")) for t, g in genes.items()}


def blastp(query_fa, subject_fa, out):
    sing("blastp", "-query", query_fa, "-subject", subject_fa, "-evalue", "1e-5", "-max_hsps", "1",
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qlen slen", "-out", out)
    hits = []
    for line in open(out):
        q, s, pid, ln, ev, bits, qs, qe, ss, se, ql, sl = line.split()
        qcov, scov = (int(qe) - int(qs) + 1) / int(ql), (int(se) - int(ss) + 1) / int(sl)
        ok = float(ev) <= E_MAX and float(pid) >= ID_MIN and qcov >= COV_MIN and scov >= COV_MIN
        hits.append({"q": q, "s": s, "pid": float(pid), "evalue": float(ev), "bits": float(bits), "qcov": qcov, "scov": scov, "accepted": ok})
    return hits


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes-file", help="FAMILY|ver1_id keys, one per line (default: the A/B panel families)")
    ap.add_argument("--label", default="", help="suffix for the output files (default: none = the panel outputs)")
    cli = ap.parse_args()
    keys = None
    if cli.genes_file:
        keys = {l.strip() for l in open(cli.genes_file) if l.strip()}
    WORK.mkdir(parents=True, exist_ok=True)
    models = wi_models()
    wi_prot = read_fasta(open(WI_PROT).read())
    thille_prot = read_fasta(open(THILLE_PROT).read())
    with gzip.open(HASS_PROT, "rt") as h:
        hass_prot_raw = h.read()
    hass_by_tag = {}
    for line in hass_prot_raw.splitlines():
        if line.startswith(">"):
            m = re.search(r"(MRB53_\d+)", line)
            if m:
                hass_by_tag[m.group(1)] = line[1:].split()[0]
    hass_prot = read_fasta(hass_prot_raw)

    rows = thille_rows(keys)
    found = {f"{r['family']}|{r['ver1']}" for r in rows} if keys else {r["family"] for r in rows}
    missing = sorted((keys or set(GENES)) - found)
    if missing:
        print(f"not on the Thille BED (no Thille locus to transfer): {missing}", file=sys.stderr)

    # 1. DNA transfer: upstream / gene / downstream of each Thille locus (1-based faidx regions, gene orientation)
    pieces = []
    for r in rows:
        up = (max(1, r["start"] + 1 - FLANK), r["start"]) if r["strand"] == "+" else (r["end"] + 1, r["end"] + FLANK)
        down = (r["end"] + 1, r["end"] + FLANK) if r["strand"] == "+" else (max(1, r["start"] + 1 - FLANK), r["start"])
        for part, (a, b) in (("up", up), ("gene", (r["start"] + 1, r["end"])), ("down", down)):
            pieces.append((f"{r['family']}|{r['ver1']}|{part}", f"{r['chrom']}:{a}-{b}", r["strand"]))
    query = WORK / "thille_loci.fa"
    with open(query, "w") as out:
        for name, region, strand in pieces:
            args = ["samtools", "faidx", THILLE_FASTA, region] + (["-i"] if strand == "-" else [])
            seq = "".join(sing(*args).splitlines()[1:])
            out.write(f">{name}\n{seq}\n")
    sam = run([MINIMAP2, "-a", "-t", "16", "--secondary=no", WI_FASTA, str(query)])
    aln = {}
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x904 or f[0] in aln:  # unmapped, secondary, supplementary
            aln.setdefault(f[0], None)
            continue
        ref_len = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", f[5]))
        aln[f[0]] = (f[2], int(f[3]), int(f[3]) + ref_len - 1, "-" if flag & 16 else "+", int(f[4]))

    # 2. model choice and protein check
    out_rows, bed = [], []
    for r in rows:
        key = f"{r['family']}|{r['ver1']}"
        up, gene, down = (aln.get(f"{key}|{p}") for p in ("up", "gene", "down"))
        tiled, gaps = "no", ""
        if gene and up and down and len({up[0], gene[0], down[0]}) == 1 and len({up[3], gene[3], down[3]}) == 1:
            # in the orientation of the gene on WI, upstream precedes the gene and downstream follows it
            first, last = (up, down) if gene[3] == "+" else (down, up)
            g1, g2 = gene[1] - first[2] - 1, last[1] - gene[2] - 1
            gaps = f"{g1};{g2}"
            in_order = first[2] < last[1] and first[1] <= gene[1] <= last[1]
            tiled = "no" if not in_order else ("contiguous" if g1 == 0 and g2 == 0 else "gapped")

        if r["family"] == "KCS11":  # Thille carries only the paralog g24008; query with the UoQ KCS11 proteins
            tags = [l.split("\t")[3].split("|")[2] for l in open(HASS_BED) if l.startswith("CM") and l.split("\t")[3].startswith("KCS11|")]
            qprot = {f"hass_{t}": hass_prot[hass_by_tag[t]] for t in tags}
        else:
            qprot = {f"thille_{r['model']}": thille_prot[f"{r['model']}_t1"]} if f"{r['model']}_t1" in thille_prot else {}
        qfa = WORK / f"{r['family']}_{r['ver1']}.query.fa"
        qfa.write_text("".join(f">{k}\n{v}\n" for k, v in qprot.items()))

        overlapping = []
        if gene:
            overlapping = [t for t, m in models.items() if m[0] == gene[0] and m[1] <= gene[2] and m[2] >= gene[1] and m[4]]
        best_overlap = None
        if qprot and overlapping:
            sfa = WORK / f"{r['family']}_{r['ver1']}.overlap.fa"
            sfa.write_text("".join(f">{t}\n{wi_prot[models[t][4]]}\n" for t in overlapping if models[t][4] in wi_prot))
            hits = blastp(str(qfa), str(sfa), str(sfa) + ".tsv")
            best_overlap = max(hits, key=lambda h: h["bits"], default=None)
        best_genome = None
        if qprot:
            allfa = WORK / "wi_proteins_by_tag.fa"
            if not allfa.exists():
                allfa.write_text("".join(f">{t}\n{wi_prot[m[4]]}\n" for t, m in models.items() if m[4] in wi_prot))
            hits = blastp(str(qfa), str(allfa), str(WORK / f"{r['family']}_{r['ver1']}.genome.tsv"))
            best_genome = max(hits, key=lambda h: h["bits"], default=None)

        chosen = best_overlap["s"] if best_overlap and best_overlap["accepted"] else ""
        if r["family"] == "KCS11":
            # the DNA transfer starts from the Thille paralog, so its WI locus is the paralog's; the WI KCS11 model
            # is the genome-wide best match of the UoQ KCS11 proteins instead
            chosen = best_genome["s"] if best_genome and best_genome["accepted"] else ""
        m = models.get(chosen)
        row = {
            "family": r["family"], "ver1": r["ver1"], "thille_locus": f"{r['chrom']}:{r['start'] + 1}-{r['end']}({r['strand']})",
            "query_protein": ",".join(qprot) or "none",
            "wi_gene_aln": f"{gene[0]}:{gene[1]}-{gene[2]}({gene[3]}) MAPQ {gene[4]}" if gene else "unmapped",
            "flanks_tiled": tiled, "gaps_up_down_bp": gaps,
            "wi_model": chosen or "none",
            "wi_model_span": f"{m[0]}:{m[1]}-{m[2]}({m[3]})" if m else "",
            "overlap_best": f"{best_overlap['s']} id{best_overlap['pid']:.1f} qcov{best_overlap['qcov']:.2f} scov{best_overlap['scov']:.2f} bits{best_overlap['bits']:.0f} {'accepted' if best_overlap['accepted'] else 'below_thresholds'}" if best_overlap else "none",
            "genome_best": f"{best_genome['s']} id{best_genome['pid']:.1f} qcov{best_genome['qcov']:.2f} scov{best_genome['scov']:.2f} bits{best_genome['bits']:.0f}" if best_genome else "none",
            "model_is_genome_best": "yes" if best_genome and chosen and best_genome["s"] == chosen else "no",
        }
        out_rows.append(row)
        if m:
            src = best_genome if r["family"] == "KCS11" else best_overlap
            bh = f"{src['q']}_id{src['pid']:.2f}_qcov{src['qcov']:.3f}"
            bed.append((m[0], m[1] - 1, m[2], f"{r['family']}|{r['ver1']}|{chosen}|flanks_tiled_{tiled}|{bh}", 0, m[3]))
        elif gene:
            what = f"thille_paralog_{r['model']}_locus" if r["family"] == "KCS11" else "aligned_locus"
            bed.append((gene[0], gene[1] - 1, gene[2], f"{r['family']}|{r['ver1']}|unannotated|{what}", 0, gene[3]))

    sfx = f"_{cli.label}" if cli.label else ""
    with open(f"summaries/wi_placement{sfx}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    with open(f"flowering_beds/igv/flowering_genes_wi{sfx}.bed", "w") as fh:
        fh.write('track name="flowering_genes_wi" description="A/B panel genes on the West Indian T2T (GCA_051132755.1): WI gene models (one line per gene)"\n')
        for b in sorted(bed):
            fh.write("\t".join(map(str, b)) + "\n")
    print(f"{len(out_rows)} loci; {sum(r['wi_model'] != 'none' for r in out_rows)} with a WI model", file=sys.stderr)


if __name__ == "__main__":
    main()
