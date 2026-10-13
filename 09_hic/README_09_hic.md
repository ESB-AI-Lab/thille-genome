# 09_hic: Omni-C contact maps, centromeres and telomeres

Copies (thille-hic session; origins in MANIFEST_09_hic.tsv). Run from /epool/edwin/pamericana/avoflow/assemblies/thille.hic/.

| Script | Makes (paper) | Inputs | Command |
|---|---|---|---|
| hic/sbatch_hic_v4.sh | juicer pre pairs on the curated AGP, Table 1A Hi-C rows, .hic for Juicebox | YaHS bin, v4 AGP | `sbatch hic/sbatch_hic_v4.sh` |
| hic/make_chromsizes.py | chromosome sizes for cooler | FASTA index | called by sbatch_cool_v4.sh |
| hic/sbatch_cool_v4.sh | cooler .cool/.mcool, ICE balancing | pairs | `sbatch hic/sbatch_cool_v4.sh` |
| hic/sbatch_balance_variants_v4.sh | strict-filter ICE weights (--min-nnz 50 --mad-max 3) used for Figure 4 | all.mcool | `sbatch hic/sbatch_balance_variants_v4.sh` |
| hic/plot_contacts.py | Figure 4 (with --no-text) | all.variants.mcool | `python3 hic/plot_contacts.py hic/v4/all/all.variants.mcool all <outdir> hic/v4/all/all.pairstats.tsv --genome-only --weight-name weight_strict --no-text` |
| hic/sbatch_chr12_profile.sh | chromosome 12 start contact partners (Results) | all.pairs.sorted.txt | `sbatch hic/v4/chr12_start/sbatch_chr12_profile.sh` |
| centromere/sbatch_trf_v4.sh, sbatch_trf_windows.sh, merge_trf_windows.py | tandem repeat calls (TRF) per window | v4 chromosomes | `sbatch centromere/sbatch_trf_v4.sh` |
| centromere/resolve_centromeres.py | centromere and telomere calls, arm ratios (Results, Figure 5C); satellite-first and Omni-C-first rankings | chr12.mcool, Tiberius GTF, TRF calls | `for m in satellite hic; do python3 centromere/resolve_centromeres.py hic/v4/chr12/chr12.mcool thille.hifiv4.asm.hic.scaffold.masked.chr.only.annotation.tiberius.gtf centromere/v4/trf centromere/v4/thille_centromeres.${m}first --rank $m; done` |
| centromere/call_centromeres.py | earlier single-ranking caller (v3), kept for provenance; superseded by resolve_centromeres.py | | `call_centromeres.py <chr12.mcool> <genes.gtf> <repeats.gff> <trf_dir> <out_prefix>` |
