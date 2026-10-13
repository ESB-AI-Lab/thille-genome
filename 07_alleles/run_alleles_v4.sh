#!/usr/bin/env bash
# Allele pipeline (methods.md §6-7) on the Thille v4 assembly, run inside v4work/ (v3 file names linked to v4 files) with
# THILLE_ASSEMBLY=v4 so the read sets come from mappingback/v4/. Same commands as the v3 run; outputs in
# v4work/summaries/fst_alleles/ and v4work/Results_fst_alleles/.
# Usage (from the repo root): bash scripts/run_alleles_v4.sh
set -euo pipefail
cd v4work
export THILLE_ASSEMBLY=v4 SINGULARITY_BIND=/gpool,/epool
mkdir -p summaries/fst_alleles/logs
G=genome; N=ncbi_dataset/data
LABELS=$(awk -F'\t' 'NR>1{print $1}' summaries/fst_alleles/gene_models.tsv)
for l in $LABELS; do for s in Thille WI Gwen; do echo "$l $s"; done; done \
  | xargs -P 6 -L 1 bash -c 'python3 scripts/gene_alleles_reads.py "$0" "$1" > summaries/fst_alleles/logs/v4_${0}_${1}.log 2>&1 || { echo "FAILED reads $0 $1"; exit 255; }'
echo "reads done $(date)"
python3 scripts/gene_alleles_assembly.py --name Puebla_hap1 --fasta $G/puebla/GCA_060308345.1/$N/GCA_060308345.1/GCA_060308345.1_ASM6030834v1_genomic.fna --threads 8 > summaries/fst_alleles/logs/v4_Puebla_hap1.log 2>&1 &
python3 scripts/gene_alleles_assembly.py --name Puebla_hap2 --fasta $G/puebla/GCA_060308325.1/$N/GCA_060308325.1/GCA_060308325.1_ASM6030832v1_genomic.fna --threads 8 > summaries/fst_alleles/logs/v4_Puebla_hap2.log 2>&1 &
wait
python3 scripts/gene_alleles_assembly.py --name Hass_pri --fasta $G/hass_pseudohap/GCA_031208685.1/$N/GCA_031208685.1/GCA_031208685.1_ASM3120868v1_genomic.fna --threads 8 > summaries/fst_alleles/logs/v4_Hass_pri.log 2>&1 &
python3 scripts/gene_alleles_assembly.py --name Hass_alt --fasta $G/hass_pseudohap/GCA_031208655.1/$N/GCA_031208655.1/GCA_031208655.1_ASM3120865v1_genomic.fna --threads 8 > summaries/fst_alleles/logs/v4_Hass_alt.log 2>&1 &
wait
python3 scripts/gene_alleles_assembly.py --name Hass --fasta $G/hass_GCA_029852735.1/hass_UoQ_chr_only.fasta --threads 8 > summaries/fst_alleles/logs/v4_Hass.log 2>&1
echo "assemblies done $(date)"
python3 scripts/gene_allele_summary.py
for l in $LABELS; do python3 scripts/plot_gene_matrix.py "$l" > /dev/null; done
python3 scripts/allele_class_pattern.py
python3 scripts/plot_a_shared_alleles.py
python3 scripts/af_skew_check.py
echo ALLELES_DONE
