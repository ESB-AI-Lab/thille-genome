# 12_population: Puebla projected onto the Gwen-paper PCAs

Copies (thille-hic session; origins in MANIFEST_12_population.tsv). Run from
/epool/edwin/pamericana/avoflow/assemblies/thille.hic/. Published panels: 2021 contig PCA
(gwen/snpcalling/plink, noOG, 30 samples, 34,472,447 variants) and 2022 scaffold PCAs
(gwen/scaffolding/scaffolds/plink/amcfilter/newanalysis, noOVL and noOVL.onlypsch); each is reproduced in PLINK 2
(r = 1.000000, identical eigenvalues) before projection. Admixture: K = 3 only (Edwin, 2026-10-11).

| Script | Makes | Inputs | Command |
|---|---|---|---|
| sbatch_panel_pca.sh, sbatch_panel_pca2.sh | contig panel pgen, allele counts, PCA allele weights | panel VCF and .map | `sbatch --mem=150G puebla_pca/sbatch_panel_pca2.sh` |
| sbatch_panel_pca_scf.sh | same for the two scaffold panels | newanalysis VCFs and .map | `sbatch --mem=150G puebla_pca/sbatch_panel_pca_scf.sh` |
| make_sites_pub.sh, sbatch_sites_union_scf.sh | bcftools targets files of the published sites | panel pvar | `bash puebla_pca/make_sites_pub.sh`; `sbatch puebla_pca/sbatch_sites_union_scf.sh` |
| sbatch_align_puebla_gwen.sh, sbatch_align_controls.sh, sbatch_align_scaffolds.sh | Puebla HiFi, Thille HiFi, Gwen CLR on the Gwen contigs and scaffolds | reads, Gwen contigs/scaffolds | `sbatch puebla_pca/<script>` |
| sbatch_genotype.sh, sbatch_genotype_scf.sh | genotypes at the published sites (bcftools 1.24 call -m -C alleles, DP >= 5) | BAMs, targets | `sbatch puebla_pca/<script>` |
| project.py, sbatch_project.sh | projection onto the published axes (PLINK 2 --score, variance-standardize; per-PC calibration on the panel) | panel dir, published eigenvec, genotypes | `sbatch puebla_pca/sbatch_project.sh <panel_dir> <published.eigenvec> <geno_dir> <out_dir>` |
