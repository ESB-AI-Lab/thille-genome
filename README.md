# Scripts for the Thille genome manuscript

Copies of the scripts behind the results, figures, and tables of the manuscript, grouped by analysis (Edwin,
2026-10-11: copied, not moved; the working copies stay where they are). Copied from thille_paper commit f0cfbc6; each
copy is byte-identical to its source listed below.

Scripts are written to run from the root of the working repository (paths such as `summaries/`, `cache/`, `genome/`
are relative to it), mostly under SLURM, with the tool versions given in the manuscript Methods and in
`containers/VERSIONS*.txt`. Helpers that other scripts import (`fst_common.py`, `coverage_class_summary.py`,
`gene_alleles_reads.py`, `gene_allele_summary.py`, `sdmyb_allele_tables.py`, `figure_paths_v4.py`,
`figure_legends.py`) sit in the same folder as the scripts that import them. Usage is in each script's docstring
(`python3 <script> --help`).

Assembly, Omni-C, comparative genomics (SyRI, plotsr, dot plots, ModDotPlot, inversions, karyotype) and the Puebla
PCA/admixture (K = 3) were run in the assembly working directory and are added by that work stream to `02_assembly/`,
`09_hic/`, `11_comparative/`, and `12_population/`; each of those folders has a `README_<folder>.md` (script, output, inputs, command) and a `MANIFEST_<folder>.tsv` (origin paths).


## 01_read_qc: Read quality

| Script | Used for | Source |
|---|---|---|
| `read_qv_hist.py` | Figure 1B (planned) / Methods: per-read QV histogram | `scripts/read_qv_hist.py` |
| `read_qv_fraction.py` | Results: fraction of reads at Q20/Q30 | `scripts/read_qv_fraction.py` |

## 02_assembly: Assembly

| Script | Used for | Source |
|---|---|---|
| `plot_cumulative_length.py` | Figure 1 (HapSolo2 assembly_stats.py, no title) | `scripts/plot_cumulative_length.py` |
| `sbatch_hapsolo_thille.sh` | HapSolo2 on the Thille hifiasm Hi-C contigs (haplotig check) | `scripts/sbatch_hapsolo_thille.sh` |

## 03_annotation: Annotation and repeats

| Script | Used for | Source |
|---|---|---|
| `annotation_table.py` | Table 2 (annotation and repeat summary) | `scripts/annotation_table.py` |
| `genome_composition_donut.py` | Figure 2 (genome composition) | `scripts/genome_composition_donut.py` |
| `swissprot_homology.py` | Gene naming by Swiss-Prot homology (Tables 4-6) | `scripts/swissprot_homology.py` |
| `gene_id_map.py` | Thille gene IDs v3 to v4 | `scripts/gene_id_map.py` |
| `relabel_v4_ids.py` | Relabel tables with v4 gene IDs | `scripts/relabel_v4_ids.py` |

## 04_read_mapping: Read mapping (Table 3)

| Script | Used for | Source |
|---|---|---|
| `map_v4.sh` | Read alignment to Thille v4 (Table 3) | `scripts/map_v4.sh` |
| `mapping_rerun.sh` | Thille HiFi re-mapping (Table 3) | `scripts/mapping_rerun.sh` |
| `mapping_stats_full.py` | Table 3 mapping statistics | `scripts/mapping_stats_full.py` |
| `identity_stats.py` | Table 3 alignment identity | `scripts/identity_stats.py` |
| `identity_all.sh` | Table 3 alignment identity, all BAMs | `scripts/identity_all.sh` |
| `stats_all.sh` | Table 3 samtools stats check | `scripts/stats_all.sh` |
| `stats_existing.sh` | Table 3 samtools stats check, existing BAMs | `scripts/stats_existing.sh` |
| `sbatch_mapping_stats_v4.sh` | Table 3 (SLURM driver, v4) | `scripts/sbatch_mapping_stats_v4.sh` |
| `aggregate_stats.py` | Table 3 aggregation | `summaries/mapping_stats/aggregate_stats.py` |
| `run_stats_pairs.sh` | Table 3 driver | `summaries/mapping_stats/run_stats_pairs.sh` |
| `run_stats_parallel.sh` | Table 3 driver | `summaries/mapping_stats/run_stats_parallel.sh` |
| `run_stats_wi.sh` | Table 3 driver (Unnamed Lowland) | `summaries/mapping_stats/run_stats_wi.sh` |

## 05_coverage: Read coverage at candidate genes

| Script | Used for | Source |
|---|---|---|
| `plot_gene_stacks.py` | Figure 7, Supplementary Figures S6, S7, S13, S14 (coverage stacks) | `plot_gene_stacks.py` |
| `make_coverage_figures_v4_notext.sh` | Driver for the coverage figures (no-text renders) | `scripts/make_coverage_figures_v4_notext.sh` |
| `sbatch_coverage_figures_v4_notext.sh` | SLURM wrapper for the coverage figures | `scripts/sbatch_coverage_figures_v4_notext.sh` |
| `genomecov_track.sh` | Per-position depth tracks (BEDTools genomecov) | `scripts/genomecov_track.sh` |
| `uoq_genomecov.sh` | Per-position depth on Hass (UoQ) | `scripts/uoq_genomecov.sh` |
| `sbatch_genomecov_puebla_v4.sh` | Per-position depth, Puebla on Thille v4 | `scripts/sbatch_genomecov_puebla_v4.sh` |
| `coverage_class_summary.py` | Coverage classes genome-wide and per gene; Supplementary Table S3 | `scripts/coverage_class_summary.py` |
| `coverage_class_runs.py` | Coverage class runs as BED | `scripts/coverage_class_runs.py` |
| `coverage_sv_hits.py` | Coverage SV hits at candidate genes (Table 6) | `scripts/coverage_sv_hits.py` |
| `genome_wide_coverage.py` | Supplementary Figure S8 and Table S3 | `scripts/genome_wide_coverage.py` |
| `ab_window_classes.py` | Coverage classes over the figure windows | `scripts/ab_window_classes.py` |
| `ab_panel_comparisons.py` | Coverage comparisons for the A/B panel genes | `scripts/ab_panel_comparisons.py` |
| `wi_segment_check.py` | Coverage-hit segments in the Unnamed Lowland assembly | `scripts/wi_segment_check.py` |
| `place_on_wi.py` | Gene placement on the Unnamed Lowland assembly | `scripts/place_on_wi.py` |
| `sv_near_regions.py` | SVs near the candidate genes | `scripts/sv_near_regions.py` |
| `sv_region_map.py` | SV maps of the chr3/chr10 regions | `scripts/sv_region_map.py` |

## 06_candidate_genes: Candidate genes and clusters

| Script | Used for | Source |
|---|---|---|
| `make_gene_tables.py` | Tables 4-6 (candidate gene tables) | `scripts/make_gene_tables.py` |
| `florid_screen/gather_gene_evidence.py` | FLOR-ID screen input: per-gene evidence table (results/final_candidates/evidence_long.tsv) | `/epool/edwin/pamericana/gwen/annotation/fstpeaksearch/scripts/gather_gene_evidence.py` (fstpeaksearch 763b601) |
| `florid_screen/flowering_screen.py` | FLOR-ID screen of the top 1% Fst genes (T1 FLOR-ID, T2 GO incl. light signalling, T3 keywords); results/final_candidates/screen.tsv, screen_v2.tsv; FLOR-ID list reference/florid/florid_flowering_genes.tsv; methods docs/methods_flowering_genes_v1.md | `/epool/edwin/pamericana/gwen/annotation/fstpeaksearch/scripts/flowering_screen.py` (commits 44430c1, 00dee42) |
| `cluster_te_sv.py` | TE and SV density at the candidate clusters | `scripts/cluster_te_sv.py` |
| `snp_markers/sbatch_blastn_snps_thille.sh` | Flowering SNP markers (SHRSPa) placed on Thille v4 by BLASTn; assembly changed from `thille.hifiv3.asm.p_ctg.fasta` to `thille.hifiv4.asm.hic.scaffold.chr.only.fasta` (Edwin, 2026-10-13), otherwise identical | `old/SNP_markers/sbatch_blastn_snps_thille.sh` |
| `snp_markers/sbatch_blastn_snps_UoQ.sh` | Flowering SNP markers placed on Hass (GCA_029852735.1) by BLASTn | `old/SNP_markers/sbatch_blastn_snps_UoQ.sh` |

## 07_alleles: Alleles of the candidate genes

| Script | Used for | Source |
|---|---|---|
| `fst_common.py` | Shared helpers for the allele pipeline | `scripts/fst_common.py` |
| `gene_models.py` | Gene models and windows of the candidate genes | `scripts/gene_models.py` |
| `gene_alleles_reads.py` | Phased read-based alleles (Thille, Unnamed Lowland, Gwen) | `scripts/gene_alleles_reads.py` |
| `gene_alleles_assembly.py` | Alleles in the Puebla and Hass haplotype assemblies | `scripts/gene_alleles_assembly.py` |
| `gene_allele_summary.py` | Allele letters and AF (Supplementary Table S4 input) | `scripts/gene_allele_summary.py` |
| `allele_class_pattern.py` | Allele classes under the SDMYB model | `scripts/allele_class_pattern.py` |
| `af_skew_check.py` | Supplementary Table S5 (AF-skew structural check) | `scripts/af_skew_check.py` |
| `plot_gene_matrix.py` | Figure 6B-D, Supplementary Figures S5, S11 (allele matrices) | `scripts/plot_gene_matrix.py` |
| `plot_a_shared_alleles.py` | Figure 6A and Supplementary Table S4 | `scripts/plot_a_shared_alleles.py` |
| `run_alleles_v4.sh` | Driver for the allele pipeline on v4 | `scripts/run_alleles_v4.sh` |
| `indel_clusters.py` | Indel clustering across reads | `scripts/indel_clusters.py` |
| `breakpoint_scan.py` | Clipped-read breakpoints (AF-skew check) | `scripts/breakpoint_scan.py` |
| `sdmyb_alleles_from_fst.py` | SDMYB allele files from the general pipeline | `scripts/sdmyb_alleles_from_fst.py` |
| `sdmyb_allele_tables.py` | SDMYB allele tables | `scripts/sdmyb_allele_tables.py` |
| `puebla_sdmyb_alleles.py` | SDMYB alleles of the Puebla haplotypes | `scripts/puebla_sdmyb_alleles.py` |
| `puebla_variants.py` | Puebla variants at SDMYB | `scripts/puebla_variants.py` |
| `sbatch_kcs11_deletion.sh` | KCS11 deletion across the haplotype assemblies | `scripts/sbatch_kcs11_deletion.sh` |
| `allele_phylo.py` | Allele phylogenies and panel dosage (in progress) | `scripts/allele_phylo.py` |
| `sdmyb_kmer_test.py` | SDMYB raw-read 21-mer test (Methods, Results): rebuilds the k-mer design and summarizes counts | `scripts/sdmyb_kmer_test.py` |
| `sbatch_sdmyb_kmer_test.sh` | SDMYB 21-mer counting in the raw reads of Thille, Unnamed Lowland, Gwen, and Puebla | `scripts/sbatch_sdmyb_kmer_test.sh` |

## 08_tandem_arrays: Tandem arrays

| Script | Used for | Source |
|---|---|---|
| `summarize_tandem_groups.py` | Table 7 (tandem arrays) | `scripts/summarize_tandem_groups.py` |
| `tandem_array_alignments.py` | Pairwise alignment of tandem copies | `scripts/tandem_array_alignments.py` |
| `tandem_member_te.py` | TE overlap of tandem copies | `scripts/tandem_member_te.py` |

## 09_hic: Omni-C (Hi-C)

| Script | Used for | Source |
|---|---|---|
| `hic_trans_contacts.py` | Inter-chromosomal Omni-C contacts | `scripts/hic_trans_contacts.py` |

## 10_figures_manuscript: Figures and manuscript build

| Script | Used for | Source |
|---|---|---|
| `figure_paths_v4.py` | Paths of every figure panel and table | `scripts/figure_paths_v4.py` |
| `sync_figures.py` | figures/ folder sync | `scripts/sync_figures.py` |
| `write_figure_legends.py` | Figure text md files | `scripts/write_figure_legends.py` |
| `figure_legends.py` | Figure text into the manuscript | `manuscript/build_scripts/figure_legends.py` |
| `build_latex.py` | LaTeX main and supplemental PDFs | `manuscript/build_scripts/build_latex.py` |


