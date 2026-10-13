#!/usr/bin/env python3
"""Plan v19 Q4: project long-read samples onto a published Gwen-paper PCA.

For one panel: (1) merge each sample's genotype chunks (bcftools concat + sort), (2) import into PLINK 2 keeping only
the published variant IDs (CHROM:POS), (3) score with the PLINK 2 projection recipe (--read-freq <panel>.acount,
--score <panel_pcs>.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15),
(4) score the panel samples with the same command and fit, per PC, published eigenvec = a * score + b by least
squares over the panel; apply a and b to the projected samples, which puts them on the published axes (the PLINK 2
docs note projected PCs differ in scale from .eigenvec). Writes <out>/projection.tsv and <out>/calibration.tsv.

Usage (from assemblies/thille.hic): project.py <panel_dir> <published.eigenvec> <geno_dir> <out_dir> [--samples ...]
"""
import argparse
import glob
import os
import statistics as st
import subprocess

BC = "singularity exec --no-home /gpool/bin/singularityimages/bcftools/bcftools_1.24.sif bcftools"
PL = "singularity exec --no-home -B {b}:{b}:ro /gpool/bin/singularityimages/plink2/plink2_2.0.0a.7.8.sif plink2"


def sh(cmd):
    print("+", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def read_sscore(p):
    rows = {}
    with open(p) as fh:
        hdr = fh.readline().lstrip("#").split()
        for l in fh:
            f = dict(zip(hdr, l.split()))
            rows[f["IID"]] = ([float(f[f"PC{k}_AVG"]) for k in range(1, 11)], int(f["ALLELE_CT"]))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel"); ap.add_argument("eigenvec"); ap.add_argument("geno"); ap.add_argument("out")
    ap.add_argument("--samples", nargs="+", default=["puebla.hifi", "thille.hifi", "gwen.clr"])
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pl = PL.format(b=os.path.dirname(os.path.abspath(a.eigenvec)))
    score = (f"--read-freq {a.panel}/panel_pub.acount --score {a.panel}/panel_pcs.eigenvec.allele 2 5 header-read "
             f"no-mean-imputation variance-standardize --score-col-nums 6-15 --allow-extra-chr --threads {a.threads}")
    # SNPs only: bcftools writes indel alleles with extra flanking bases (e.g. TGTCCGTCCGT/TGTCCGT for the panel's
    # TGTCC/T), so indels cannot be matched by allele; the panel is scored on the same SNP set for the calibration.
    # also drop variants with no information in the panel (no called genotypes, or ALT absent or fixed among the
    # panel calls): their PCA weights are zero, and PLINK 2 cannot variance-standardize a sample dosage there
    informative = set()
    with open(f"{a.panel}/panel_pub.acount") as fh:
        for l in fh:
            if not l.startswith("#"):
                f = l.split()
                if 0 < int(f[4]) < int(f[5]):
                    informative.add(f[1])
    panel_alleles = {}
    with open(f"{a.panel}/panel_pub.pvar") as fh, open(f"{a.out}/panel_snps.txt", "w") as out:
        for l in fh:
            if l.startswith("#"):
                continue
            f = l.split("\t")
            if len(f[3]) == 1 and len(f[4]) == 1 and f[2] in informative:
                panel_alleles[f[2]] = (f[3], f[4])
                out.write(f[2] + "\n")
    print(f"panel SNPs used for projection: {len(panel_alleles)} (informative variants: {len(informative)})", flush=True)
    sh(f"{pl} --pfile {a.panel}/panel_pub --extract {a.out}/panel_snps.txt {score} --out {a.out}/panel_scored")
    for s in a.samples:
        vcf = f"{a.out}/{s}.vcf.gz"
        if not os.path.exists(vcf):
            parts = sorted(glob.glob(f"{a.geno}/{s}.*.vcf.gz"))
            if not parts:
                raise SystemExit(f"no genotype chunks for {s} in {a.geno}")
            # bcftools call -i writes every target site into every chunk (missing outside the chunk's own
            # contigs), so each chunk is first restricted to its own region: the chunk .bed (contig run) or the
            # scaffold named in the file (scaffold run).
            own = []
            for part in parts:
                bed = part[:-len(".vcf.gz")] + ".bed"
                region = f"-T {bed}" if os.path.exists(bed) else f"-t {part.split('/')[-1][len(s) + 1:-len('.vcf.gz')]}"
                bcf = f"{a.out}/{os.path.basename(part)[:-len('.vcf.gz')]}.own.bcf"
                sh(f"{BC} view {region} -Ob -o {bcf} {part}")
                own.append(bcf)
            os.makedirs(f"{a.out}/sorttmp_{s}", exist_ok=True)
            sh(f"{BC} concat -Ou {' '.join(own)} | {BC} sort -m 8G -T {a.out}/sorttmp_{s} -Oz -o {vcf}.tmp && mv {vcf}.tmp {vcf}")
            n = subprocess.run(f"{BC} view -H {vcf} | wc -l", shell=True, check=True, capture_output=True, text=True).stdout.strip()
            print(f"{s}: {n} sites after merging chunks", flush=True)
        sh(f"{pl} --vcf {vcf} --set-all-var-ids @:# --allow-extra-chr --extract {a.panel}/published_variant_ids.txt "
           f"--max-alleles 2 --make-pgen --out {a.out}/{s}")
        # keep SNPs whose REF matches the panel and whose ALT is the panel ALT or "." (sample carries REF only),
        # then set the panel ALT so PLINK 2 sees the same biallelic variant as in the panel counts
        kept = 0
        with open(f"{a.out}/{s}.pvar") as fh, open(f"{a.out}/{s}_matched.txt", "w") as out, \
                open(f"{a.out}/{s}_alts.tsv", "w") as alts:
            for l in fh:
                if l.startswith("#"):
                    continue
                f = l.rstrip("\n").split("\t")
                pa = panel_alleles.get(f[2])
                if pa and f[3] == pa[0] and f[4] in (pa[1], "."):
                    out.write(f[2] + "\n"); alts.write(f"{f[2]}\t{pa[1]}\n"); kept += 1
        print(f"{s}: {kept} panel SNPs with matching alleles", flush=True)
        sh(f"{pl} --pfile {a.out}/{s} --allow-extra-chr --extract {a.out}/{s}_matched.txt "
           f"--alt-allele {a.out}/{s}_alts.tsv 2 1 --make-pgen --out {a.out}/{s}_snps")
        sh(f"{pl} --pfile {a.out}/{s}_snps {score} --out {a.out}/{s}_scored")
    pub = {l.split()[1]: [float(x) for x in l.split()[2:12]] for l in open(a.eigenvec) if not l.startswith("#")}
    panel = read_sscore(f"{a.out}/panel_scored.sscore")
    ids = [i for i in panel if i in pub]
    cal = []
    with open(f"{a.out}/calibration.tsv", "w") as fh:
        fh.write("PC\tslope\tintercept\tr\tn_panel\n")
        for k in range(10):
            x = [panel[i][0][k] for i in ids]; y = [pub[i][k] for i in ids]
            mx, my = st.mean(x), st.mean(y)
            sxy = sum((u - mx) * (v - my) for u, v in zip(x, y)); sxx = sum((u - mx) ** 2 for u in x)
            syy = sum((v - my) ** 2 for v in y)
            slope = sxy / sxx; icpt = my - slope * mx; r = sxy / (sxx * syy) ** 0.5
            cal.append((slope, icpt)); fh.write(f"PC{k + 1}\t{slope:.6g}\t{icpt:.6g}\t{r:.6f}\t{len(ids)}\n")
    with open(f"{a.out}/projection.tsv", "w") as fh:
        fh.write("sample\tsource\tallele_ct\t" + "\t".join(f"PC{k}" for k in range(1, 11)) + "\n")
        for i in ids:
            fh.write(f"{i}\tpanel_published\t\t" + "\t".join(f"{v:.6f}" for v in pub[i]) + "\n")
        for s in a.samples:
            sc = read_sscore(f"{a.out}/{s}_scored.sscore")
            for iid, (v, n) in sc.items():
                fh.write(f"{s}\tprojected\t{n}\t" + "\t".join(f"{cal[k][0] * v[k] + cal[k][1]:.6f}" for k in range(10)) + "\n")
    print(open(f"{a.out}/calibration.tsv").read())


if __name__ == "__main__":
    main()
