#!/usr/bin/env python3
"""Variants and CDS of each Fst flowering gene in the assembled haplotypes (Puebla hap1/hap2, Hass UoQ).

1. Locate: all Thille windows (summaries/fst_alleles/gene_models.tsv) are mapped to the assembly with minimap2 2.30
   (-x asm20 --secondary=no -t --threads). Per window, alignments on the same-numbered chromosome are kept, the
   block with the most matching bases is the seed, and blocks within --chain bp of it (same strand) form the locus;
   the assembly region is the locus span +- 2 kb. A window is "placed" when its blocks cover >= --min-cov of it.
2. Align: the assembly region is aligned back to the Thille window (target) with minimap2 -c --cs=long -x asm20, so
   every Thille position gets the assembly's base (or a gap) and insertions, as in scripts/puebla_variants.py.
3. Variants: SNPs, and insertions/deletions >= 2 bp from the cs tags; unaligned Thille stretches between consecutive
   blocks are one deletion, unaligned assembly stretches one insertion (placed at the block junction).
4. CDS: the assembly's bases over the CDS exons in transcript order (1-bp indels included, insertions inside an exon
   added), oriented to the gene strand. An exon position with no aligned base makes the CDS unknown ("N").

Usage (from the repo root):
    python3 scripts/gene_alleles_assembly.py --name Puebla_hap1 --fasta <fna> [--threads 16] [--labels SDMYB,AP2]
Writes summaries/fst_alleles/<label>/<name>.{paf,variants.tsv,cds.fa} and summaries/fst_alleles/placement_<name>.tsv
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fst_common import OUT, load_models, ref_window, sh, translate, cs_target_bases  # noqa: E402

MM = "/gpool/bin/minimap2-2.30_x64-linux/minimap2"
SCR = Path("summaries/fst_alleles/_assembly_tmp")


def chrom_map(fasta):
    """Thille chrN -> assembly sequence name, from the assembly's FASTA headers ('chromosome N')."""
    hdr = subprocess.run(["grep", ">", fasta], capture_output=True, text=True, check=True).stdout.splitlines()
    out = {}
    for h in hdr:
        f = h[1:].split()
        if "chromosome" in f:
            n = f[f.index("chromosome") + 1].rstrip(",")
            if n.isdigit():
                out[f"chr{n}"] = f[0]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--name", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--chain", type=int, default=200000)
    p.add_argument("--min-cov", type=float, default=0.8)
    p.add_argument("--labels", help="comma-separated subset of gene labels")
    a = p.parse_args()
    models = load_models()
    if a.labels:
        models = {k: v for k, v in models.items() if k in a.labels.split(",")}
    cmap = chrom_map(a.fasta)
    SCR.mkdir(parents=True, exist_ok=True)
    wins = SCR / f"windows_{a.name}.fa"
    with open(wins, "w") as fh:
        for lab, g in models.items():
            fh.write(f">{lab}\n{ref_window(g)}\n")
    loc = subprocess.run([MM, "-x", "asm20", "--secondary=no", "-t", str(a.threads), a.fasta, str(wins)],
                         capture_output=True, text=True, check=True).stdout
    hits = {}
    for l in loc.splitlines():
        f = l.split("\t")
        hits.setdefault(f[0], []).append((f[5], f[4], int(f[7]), int(f[8]), int(f[2]), int(f[3]), int(f[9])))
    place_rows = []
    for lab, g in models.items():
        want = cmap.get(g.chrom)
        hs = [h for h in hits.get(lab, []) if h[0] == want]
        status, region = "not placed", ""
        wlen = g.win_end - g.win_start + 1
        if hs:
            seed = max(hs, key=lambda h: h[6])
            chain = [h for h in hs if h[1] == seed[1] and abs(h[2] - seed[2]) <= a.chain]
            cov = sum(h[5] - h[4] for h in chain) / wlen
            s, e = max(1, min(h[2] for h in chain) - 2000), max(h[3] for h in chain) + 2000
            region = f"{want}:{s}-{e}"
            status = "placed" if cov >= a.min_cov else f"partial ({cov:.2f} of window)"
        else:
            cov = 0.0
        place_rows.append((lab, g.chrom, want or "", region, f"{cov:.3f}", status))
        d = OUT / lab; d.mkdir(parents=True, exist_ok=True)
        if not region:
            (d / f"{a.name}.variants.tsv").write_text("status\tnot placed\n")
            continue
        sub = SCR / f"{a.name}_{lab}.fa"
        sub.write_text(sh("samtools", "faidx", a.fasta, region))
        tw = SCR / f"{a.name}_{lab}_thille.fa"
        tw.write_text(f">thille_window\n{ref_window(g)}\n")
        paf = d / f"{a.name}.paf"
        with open(paf, "w") as fh:
            subprocess.run([MM, "-c", "--cs=long", "-x", "asm20", "-t", "4", str(tw), str(sub)], stdout=fh,
                           stderr=subprocess.DEVNULL, check=True)
        base, ins, blocks = cs_target_bases(paf, g.win_start)
        ref = ref_window(g)
        rb = lambda pos: ref[pos - g.win_start]
        rows = []
        for pos in sorted(base):
            b = base[pos]
            if b not in ("-",) and b != rb(pos):
                rows.append(("SNP", pos, rb(pos), b, 1))
        # deletions: runs of '-' (>= 2 bp), anchored at the base before the run
        run = None
        for pos in sorted(base) + [None]:
            if pos is not None and base[pos] == "-" and (run and pos == run[1] + 1):
                run[1] = pos
                continue
            if run and run[1] - run[0] + 1 >= 2:
                rows.append(("DEL", run[0] - 1, "", "", run[1] - run[0] + 1))
            run = [pos, pos] if pos is not None and base[pos] == "-" else None
        for pos, sq in ins.items():
            if len(sq) >= 2:
                rows.append(("INS", pos, "", "", len(sq)))
        for (t1s, t1e, q1s, q1e), (t2s, t2e, q2s, q2e) in zip(blocks, blocks[1:]):
            if t2s - t1e >= 2:
                rows.append(("DEL", g.win_start + t1e - 1, "", "", t2s - t1e))
            qgap = (q2s - q1e) if q2s >= q1e else (q1s - q2e)
            if qgap >= 2:
                rows.append(("INS", g.win_start + t1e - 1, "", "", qgap))
        with open(d / f"{a.name}.variants.tsv", "w") as fh:
            fh.write("type\tpos\tref\talt\tlength\n")
            for r in sorted(rows, key=lambda r: r[1]):
                fh.write("\t".join(map(str, r)) + "\n")
        # CDS
        if g.exons:
            parts = []
            for s, e in g.exons:
                fw = ""
                for pos in range(s, e + 1):
                    b = base.get(pos, "N")
                    fw += "" if b == "-" else b
                    if pos < e and pos in ins:
                        fw += ins[pos]
                parts.append(g.orient(fw))
            cds = "".join(parts)
            with open(d / f"{a.name}.cds.fa", "w") as fh:
                fh.write(f">{a.name}_cds len={len(cds)}\n{cds}\n>{a.name}_protein\n{translate(cds) if 'N' not in cds else 'unknown'}\n")
        print(f"{a.name} {lab}: {status} {region} cov {cov:.2f}; {sum(r[0]=='SNP' for r in rows)} SNPs, "
              f"{sum(r[0]=='DEL' for r in rows)} DEL, {sum(r[0]=='INS' for r in rows)} INS", flush=True)
    with open(OUT / f"placement_{a.name}.tsv", "w") as fh:
        fh.write("label\tthille_chrom\tassembly_chrom\tassembly_region\twindow_cov\tstatus\n")
        for r in place_rows:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
