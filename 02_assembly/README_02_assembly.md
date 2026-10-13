# 02_assembly: genome assembly, Hi-C scaffolding, chromosome 8 curation, QC

Copies of the scripts behind the assembly results (thille-hic session; origins in MANIFEST_02_assembly.tsv). Run from
the origin directory given; SLURM jobs ran on partitions gcbigmem,gcluster.

| Script | Makes (paper) | Inputs | Command |
|---|---|---|---|
| hifiasm/sbatch_hifiasm-hic.sh | hifiasm Hi-C mode primary contigs (Methods: Genome Assembly; Table 1B contigs) | Thille HiFi, Omni-C reads | `sbatch sbatch_hifiasm-hic.sh` (in hifi.hic.gcluster/) |
| hifiasm/bash_gfa2fasta.sh | contig FASTA from the hifiasm GFA | *.p_ctg.gfa | `bash bash_gfa2fasta.sh` |
| yahs/sbatch_bwamem.sh | Omni-C alignments for YaHS | contigs, Omni-C reads | `sbatch sbatch_bwamem.sh` |
| yahs/sbatch_yahs.sh | YaHS scaffolds (Table 1B scaffolds) | contigs, Omni-C BAM | `sbatch sbatch_yahs.sh` |
| chr8_curation/sbatch_map_hifi_yahs.sh | HiFi depth on the YaHS scaffolds (chr8 heterotype duplication evidence, Results) | YaHS FASTA, HiFi reads | `sbatch curation/chr8/sbatch_map_hifi_yahs.sh` (in thille.hic/) |
| chr8_curation/sbatch_jbat.sh | Juicebox Assembly Tools map (juicer pre -a, juicer_tools pre -n) and juicer post of the edited assembly (v4) | YaHS bin/AGP, Omni-C BAM, review.assembly | `sbatch sbatch_jbat.sh` (in yahs_thille.hifiv3.asm.hic.p_ctg.curation/) |
| naming/sbatch_align_lowland.sh, naming/name_by_lowland.py | chromosome names by alignment to Unnamed Lowland (Methods) | v4 scaffolds, GCA_051132755.1 | `sbatch curation/naming/sbatch_align_lowland.sh`; `python3 curation/naming/name_by_lowland.py <paf>` |
| stats/sbatch_hic_seqkit.sh | assembly statistics (Table 1B) | contig, scaffold, chromosome FASTA | `sbatch stats/sbatch_hic_seqkit.sh` |
| merqury/sbatch_merqury.sh, merqury/run_merqury.sh | QV, k-mer completeness, Figure 3A spectrum | HiFi meryl db (k = 21), v4 chromosomes | `sbatch merqury/sbatch_merqury.sh` |
| busco/sbatch_busco.sh, busco/run_busco_local_2026-10-04.sh | BUSCO genome and protein mode (Table 1B, Figure 3B) | v4 chromosomes, Tiberius proteins, embryophyta_odb12.2 | `sbatch sbatch_busco.sh` (in annotation/tiberius/thille.hic/) |
| busco/busco_plot.py | Figure 3B | BUSCO short summaries | `python3 busco_plot.py -i <proteins summary> <genome summary> --labels "Thille Proteins" "Thille HiC Scaffolds" --no-title --annotation none -o Fig3B_busco` |
| annotation_stats/density_stats.py | gene and TE density at chromosome ends and centromeres (Results) | Tiberius GTF, TE GFF, centromere BED, FASTA index | `python3 compare/v4/density_stats.py <gtf> <TE.gff> centromere/v4/thille_centromeres.satellitefirst.bed <fai>` |
