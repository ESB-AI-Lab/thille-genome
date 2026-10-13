#!/usr/bin/env python3
"""Allele phylogenies and panel genotypes for SDMYB and the five genes with an A-shared allele (Edwin, 2026-10-08;
PLAN.md "SDMYB phylogeny and the missing allele analyses", P2 to P5).

For each gene:
 1. Thille CDS from the v4 genome and the Tiberius exons (gene_models.tsv; minus strand reverse-complemented).
 2. Its locus on the Gwen primary contigs (the reference of the phenotyped Illumina panel of the Gwen study) by
    spliced alignment (minimap2 -ax splice --secondary=no); the aligned exon blocks on the contig.
 3. Each panel sample (GATK-recalibrated BAMs): bcftools mpileup | call -m over the blocks (SNPs only, QUAL >= 20,
    base and mapping quality >= 20); the sample CDS is the contig sequence of the blocks with homozygous ALT applied,
    heterozygous sites as IUPAC codes, and positions with depth < 5 (samtools depth, same quality filters) as N,
    oriented to the gene. Indels are not applied (Illumina, unphased); they are carried by the phased haplotypes only.
 4. The 10 phased haplotype CDS (alleles.fa of the allele analysis) plus the panel and outgroup CDS -> MAFFT (--auto)
    -> IQ-TREE 3 (-m MFP, 1,000 ultrafast bootstraps, rooted on Ocotea_bothanthra).
 5. Diagnostic sites of the A-shared allele: alignment columns at which the A-shared allele carries a base that no
    B-type haplotype (Thille h1/h2, Unnamed Lowland h1/h2) carries. Per panel sample: number of diagnostic sites with the
    A-shared base homozygous, heterozygous (IUPAC containing it), absent, or no data; dosage 2/1/0 by majority.
Outputs: summaries/v4/allele_phylo/<gene>/ (cds_all.fa, aln.fa, iqtree files) and summaries/v4/allele_phylo/panel_dosage.tsv
Usage (from the repo root, on SLURM): python3 scripts/allele_phylo.py [--threads 8]
"""
import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

REPO = Path("/gpool/edwin/avocados/thille_paper")
GENES = ["SDMYB", "AGAMOUS-like", "CSU2", "NPR5-BOP2", "SGRL", "LCYB"]
SHARED = {"SDMYB": "b", "AGAMOUS-like": "c", "CSU2": "d", "NPR5-BOP2": "d", "SGRL": "d", "LCYB": "d"}  # Table S4
MODELS = REPO / "v4work/summaries/fst_alleles/gene_models.tsv"
ALLELES = REPO / "v4work/summaries/fst_alleles/{g}/alleles.fa"
THILLE = REPO / "genome/thille.hifiv4.asm.hic.scaffold.masked.chr.only.fasta"
GWEN = Path("/epool/edwin/pamericana/avoflow/assemblies/thille.hic/pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta")
BAMDIR = Path("/epool/edwin/pamericana/gwen/snpcalling/gatkrecal")
BAMFMT = "pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary_{s}{suffix}_recal.bam"
TYPE_B = ["Bacon", "Fuerte_1", "Ganter", "Linda", "Lyon", "Nabal", "Nimlioh", "Thille", "Zutano"]
TYPE_A = ["Anaheim", "Carlsbad", "Fairchild", "Gwen", "Mendez", "Pinkerton", "Reed", "Simmonds", "Taft", "TopaTopa",
          "VC26", "Waldin"]            # Hass-2 and Hass-R excluded (Edwin, 2026-10-08: UoQ-only Hass rule)
OUTGROUPS = ["Persea_hintonii", "Persea_donnell_smithii", "Ocotea_bothanthra"]
SAMTOOLS = "/gpool/bin/samtools-1.23.1/bin/samtools"
MINIMAP2 = "/gpool/bin/minimap2-2.30_x64-linux/minimap2"
BCF = ["singularity", "exec", "-B", "/gpool,/epool", str(REPO / "containers/bcftools_1.24.sif"), "bcftools"]
MAFFT = ["singularity", "exec", "-B", "/gpool,/epool", str(REPO / "containers/mafft_7.525.sif"), "mafft"]
IQTREE = ["singularity", "exec", "-B", "/gpool,/epool", str(REPO / "containers/iqtree_3.1.4.sif"), "iqtree3"]
OUT = REPO / "summaries/v4/allele_phylo"
COMP = str.maketrans("ACGTNRYKMSWBDHVacgtn", "TGCANYRMKSWVHDBtgcan")
IUPAC = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CG"): "S",
         frozenset("AT"): "W"}
IUPAC_BASES = {v: set(k) for k, v in IUPAC.items()}


def run(cmd, **kw):
    r = subprocess.run(cmd, text=True, capture_output=True, **kw)
    if r.returncode:
        raise RuntimeError(f"{' '.join(map(str, cmd))[:300]}\n{r.stderr[-2000:]}")
    return r.stdout


def faidx(fa, region):
    return "".join(run([SAMTOOLS, "faidx", str(fa), region]).split("\n")[1:])


def rc(s):
    return s.translate(COMP)[::-1]


def bam(sample):
    for suffix in (".merged_sorted", "_sorted"):
        p = BAMDIR / BAMFMT.format(s=sample, suffix=suffix)
        if p.exists():
            return p
    raise FileNotFoundError(f"no recalibrated BAM for {sample} in {BAMDIR}")


def thille_cds(row):
    exons = [tuple(map(int, x.split("-"))) for x in row["cds_exons"].split(",")]
    seq = "".join(faidx(THILLE, f"{row['chrom']}:{s}-{e}") for s, e in exons)
    return rc(seq) if row["strand"] == "-" else seq


def gwen_blocks(gene, cds, d):
    """Spliced alignment of the Thille CDS to the Gwen contigs; returns contig, strand, [(start,end)] 1-based, identity."""
    q = d / "thille_cds.fa"
    q.write_text(f">{gene}\n{cds}\n")
    sam = run([MINIMAP2, "-ax", "splice", "--secondary=no", "-t", "4", str(GWEN), str(q)])
    rec = next(l.split("\t") for l in sam.splitlines() if l and not l.startswith("@") and not int(l.split("\t")[1]) & 0x904)
    flag, contig, pos, cigar = int(rec[1]), rec[2], int(rec[3]), rec[5]
    nm = int(re.search(r"NM:i:(\d+)", "\t".join(rec[11:])).group(1))
    blocks, cur, ref, aln = [], pos, pos, 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(n)
        if op in "M=X":
            ref += n; aln += n
        elif op == "D":
            ref += n
        elif op == "N":
            blocks.append((cur, ref - 1)); ref += n; cur = ref
    blocks.append((cur, ref - 1))
    return contig, "-" if flag & 16 else "+", blocks, 1 - nm / max(aln, 1)


def sample_cds(sample, contig, strand, blocks, d):
    regions = ",".join(f"{contig}:{s}-{e}" for s, e in blocks)
    bed = d / "blocks.bed"
    bed.write_text("".join(f"{contig}\t{s - 1}\t{e}\n" for s, e in blocks))
    b = str(bam(sample))
    vcf = run(BCF[:-1] + ["sh", "-c", f"bcftools mpileup -f {GWEN} -R {bed} -q 20 -Q 20 -a AD,DP {b} 2>/dev/null | "
                                      f"bcftools call -m -V indels 2>/dev/null | bcftools view -i 'QUAL>=20' -H"])
    gt = {}
    for l in vcf.splitlines():
        f = l.split("\t")
        alleles = [f[3]] + f[4].split(",")
        g = re.split(r"[/|]", f[9].split(":")[0])
        if "." in g:
            continue
        bases = {alleles[int(x)] for x in g}
        if any(len(x) != 1 for x in bases):
            continue
        gt[int(f[1])] = bases
    depth = {}
    for l in run([SAMTOOLS, "depth", "-a", "-q", "20", "-Q", "20", "-b", str(bed), b]).splitlines():
        c, p, dp = l.split("\t"); depth[int(p)] = int(dp)
    seq = []
    for s, e in blocks:
        ref = faidx(GWEN, f"{contig}:{s}-{e}").upper()
        for i, base in enumerate(ref):
            p = s + i
            if depth.get(p, 0) < 5:
                seq.append("N")
            elif p in gt:
                bs = gt[p]
                seq.append(next(iter(bs)) if len(bs) == 1 else IUPAC.get(frozenset(bs), "N"))
            else:
                seq.append(base)
    s = "".join(seq)
    return rc(s) if strand == "-" else s


def read_fa(p):
    out, name = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            name = l[1:].split()[0]; out[name] = ""
        elif name:
            out[name] += l
    return out


def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument("--threads", type=int, default=8)
    args = a.parse_args()
    models = {r["label"]: r for r in csv.DictReader(open(MODELS), delimiter="\t")}
    dosage_rows = []
    for g in GENES:
        d = OUT / g
        d.mkdir(parents=True, exist_ok=True)
        cds = thille_cds(models[g])
        contig, strand, blocks, ident = gwen_blocks(g, cds, d)
        print(f"{g}: Thille CDS {len(cds)} bp -> Gwen {contig} {strand} {blocks[0][0]}-{blocks[-1][1]} "
              f"({len(blocks)} blocks, identity {ident:.3f})", flush=True)
        hap = read_fa(ALLELES.format(g=g))
        seqs = {k: v for k, v in hap.items() if re.match(r"^(Thille|WI|Gwen|Puebla|Hass)_h[12]$", k)}
        seqs = {k.replace("WI_", "UnnamedLowland_"): v for k, v in seqs.items()}
        for s in TYPE_B + TYPE_A + OUTGROUPS:
            grp = "B" if s in TYPE_B else "A" if s in TYPE_A else "out"
            seqs[f"panel_{grp}_{s}"] = sample_cds(s, contig, strand, blocks, d)
        fa = d / "cds_all.fa"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        aln = d / "aln.fa"
        aln.write_text(run(MAFFT + ["--auto", "--quiet", str(fa)]))
        run(IQTREE + ["-s", str(aln), "-m", "MFP", "-B", "1000", "-T", str(args.threads), "-o", "panel_out_Ocotea_bothanthra",
                      "--prefix", str(d / "iqtree"), "-redo", "--quiet"])
        # diagnostic columns of the A-shared allele
        al = read_fa(aln)
        allele_of = {k.split()[0]: v.split("=")[-1] for k, v in
                     ((l[1:].split()[0], l.split("allele=")[1].split()[0]) for l in open(ALLELES.format(g=g))
                      if l.startswith(">") and "allele=" in l)}
        shared = f"{g}^{SHARED[g]}"
        a_haps = [k for k in al if not k.startswith("panel_") and allele_of.get(k.replace("UnnamedLowland_", "WI_")) == shared]
        b_haps = [k for k in al if k.startswith(("Thille_", "UnnamedLowland_"))]
        L = len(next(iter(al.values())))
        diag = []
        for i in range(L):
            a_base = {al[k][i].upper() for k in a_haps}
            b_base = {al[k][i].upper() for k in b_haps}
            if len(a_base) == 1 and "-" not in a_base and not (a_base & b_base) and "N" not in b_base:
                diag.append((i, next(iter(a_base))))
        (d / "diagnostic_columns.tsv").write_text("column\tA_shared_base\n" + "".join(f"{i + 1}\t{b}\n" for i, b in diag))
        for k, v in al.items():
            if not k.startswith("panel_"):
                continue
            hom = het = absent = nod = 0
            for i, b in diag:
                x = v[i].upper()
                if x in "N-":
                    nod += 1
                elif x == b:
                    hom += 1
                elif b in IUPAC_BASES.get(x, set()):
                    het += 1
                else:
                    absent += 1
            called = hom + het + absent
            dose = "" if called == 0 else 2 if hom > max(het, absent) else 1 if het >= absent and het > 0 else 0
            dosage_rows.append({"gene": g, "sample": k.split("_", 2)[2], "type": k.split("_")[1], "diagnostic_sites": len(diag),
                                "hom_shared": hom, "het_shared": het, "absent": absent, "no_data": nod, "dosage": dose})
    with open(OUT / "panel_dosage.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dosage_rows[0]), delimiter="\t"); w.writeheader(); w.writerows(dosage_rows)
    print("done", OUT / "panel_dosage.tsv")


if __name__ == "__main__":
    main()
