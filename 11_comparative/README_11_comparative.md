# 11_comparative: comparisons with Hass (UoQ), Unnamed Lowland, Puebla and Gwen

Copies (thille-hic session; origins in MANIFEST_11_comparative.tsv). Run from
/epool/edwin/pamericana/avoflow/assemblies/thille.hic/ unless noted; SLURM partitions gcbigmem,gcluster.

| Script | Makes (paper) | Inputs | Command |
|---|---|---|---|
| lib/fastaio.py, concat_chroms.py, split_fasta.py, filter_fasta_minlen.py, wg_orientation.py | helpers (FASTA I/O, chromosome concatenation for ModDotPlot, orientation table) | | imported or called by the jobs below |
| alignment/sbatch_prep_v4.sh, sbatch_cmp_westindian_v4.sh, sbatch_puebla_prep_v4.sh | minimap2 asm5 alignments, chromosome matching and same-strand FASTAs | v4 chromosomes, Hass UoQ, GCA_051132755.1, GCA_060308345.1 | `sbatch compare/<script>` |
| alignment/sbatch_cmp_gwen_v4.sh, sbatch_prep_gwen_cigar_v4.sh, build_gwen_pseudochr.py | Gwen contig placement and Figure S2 input | Gwen primary contigs | `sbatch compare/<script>` |
| alignment/sbatch_cmp_hass_wg_v4.sh, sbatch_cmp_westindian_wg_v4.sh, sbatch_cmp_puebla_wg_v4.sh | whole-genome ModDotPlot (Figure S15) | concatenated chromosomes | `sbatch compare/<script>` |
| syri/sbatch_syri_v4.sh, sbatch_syri_threeway_v4.sh, sbatch_puebla_syri_v4.sh | SyRI 1.7.0 and plotsr (Figure 5B, Table S1 input) | minimap2 -ax asm5 --eqx BAMs | `sbatch compare/<script>` |
| syri/sv_to_bed_v4.py, sv_events_table_v4.py | SV BEDs and Supplementary Table S1 | SyRI outputs, Gwen SV candidates | `python3 compare/sv_to_bed_v4.py compare/v4/sv_bed`; `python3 compare/sv_events_table_v4.py compare/v4/sv_bed/sv_events.ge100kb.tsv` |
| syri/gwen_sv_candidates.py | Gwen SV candidates with read-depth classes (Results) | Gwen placement, Gwen CLR depth | see script usage |
| syri/syri_counts.py, shared_svs.py, thille_only_stats.py | SyRI totals, shared Puebla/Unnamed Lowland SVs, Thille-only sequence (Results) | SyRI outputs, centromere calls | `python3 compare/v4/<script>` |
| figures/plot_syri_overlay.py | Figure 5A | SyRI out, FASTA indexes | see figures/sbatch_notext.sh |
| figures/plot_genome_dotplot.py | Figure S2 | Gwen PAF, order files | see figures/sbatch_notext_part2.sh |
| figures/plot_wg_compare.py, moddotplot_stats.py | Figure S15 and ModDotPlot window statistics (Results) | ModDotPlot bedpe, boundaries | see figures/sbatch_notext.sh; `python3 compare/v4/moddotplot_stats.py <bedpe> compare/v4/wg/thille_wg.boundaries.tsv <y boundaries>` |
| figures/plot_region_tracks_v4.py, plot_region_threeway_v4.py, bin_genomecov.sh, sbatch_bin_puebla.sh, sbatch_repeatmasker_windows.sh | Figure S9 | coverage bins, SyRI, annotation, RepeatMasker windows | see figures/sbatch_notext.sh |
| figures/run_region_plotsr_v4.sh, sbatch_relabel2.sh | Figure S10 ribbons and the Unnamed Lowland-labelled plotsr panels | SyRI chain | `sbatch compare/v4/relabel2/sbatch_relabel2.sh` |
| figures/sbatch_notext.sh, sbatch_notext_part2.sh, sbatch_notext_karyotype.sh | final no-text renders of Figures 3B, 4, 5A, 5C, S2, S9, S12, S15 | as above | `sbatch compare/v4/notext/<script>` |
| karyotype/render_karyotype.py, sbatch_karyotype_dots_v4.sh, make_karyotype_dots.py | Figure 5C | chromosome sizes, GTF, TE GFF, centromere calls, flowering gene positions | `sbatch karyotype/sbatch_karyotype_dots_v4.sh` |
| inversions/make_windows.py, sbatch_inv_pileup.sh, sbatch_inv_pileup_puebla.sh, inv_pileup.py | read support at the 152 inversions (Results) | inversions.tsv, read BAMs | `sbatch compare/v4/inversions/<script>` |
| inversions/summarize_support.py, plot_inversions.py, sbatch_plot_inversions.sh, te/sbatch_repeatmasker_invwin.sh | Figure S12 and support calls | pileup bins, clips, annotation | see figures/sbatch_notext_part2.sh |
| inversions/make_inversion_beds.py | inversion BEDs with read-support calls | inversions.tsv, support_summary_v2.tsv | `python3 make_inversion_beds.py` (in compare/v4/inversions/) |
| ragtag/split_at_gaps.py, sbatch_ragtag.sh, finalize_ragtag.py, sbatch_validate.sh | Thille-guided scaffolds of Gwen and Hass (UoQ) (Results) | Gwen contigs, Hass UoQ | `sbatch sbatch_ragtag.sh` (in assemblies/ragtag_thille_v4/) |
