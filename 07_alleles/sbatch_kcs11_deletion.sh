#!/bin/bash
#SBATCH --job-name=kcs11.del
#SBATCH --partition=gcluster,gcbigmem
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/kcs11_deletion.%j.log
# KCS11 deletion (Edwin 2026-10-09): the Gwen-study KCS11 (ver1 ChrU.g440620 + g440630 = Hass UoQ MRB53_020451/020452,
# CM056814.1:21,128,267-21,133,279) has no Thille ortholog; Thille and Unnamed Lowland HiFi reads have zero depth on Hass
# over 21,107,947-21,199,536 and 21,094,382-21,199,536, Gwen CLR half depth over 21,120,839-21,162,493
# (summaries/coverage_sv_hits_hass.tsv). This job places the Hass region (+-150 kb flanks) and the KCS11 gene in every
# assembly to reconstruct the deletion: which haplotypes carry KCS11 and where the breakpoints are.
# Usage (repo root): sbatch scripts/sbatch_kcs11_deletion.sh   -> summaries/v4/kcs11_deletion/
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node $(hostname)  cpus ${SLURM_CPUS_PER_TASK:-?}  start $(date)"
ST=/gpool/bin/samtools-1.24/samtools
MM=/gpool/bin/minimap2-2.30_x64-linux/minimap2
O=summaries/v4/kcs11_deletion; mkdir -p $O
HASS=genome/hass_GCA_029852735.1/hass_UoQ_chr_only.fasta
$ST faidx $HASS CM056814.1:20944382-21349536 > $O/hass_region.fa     # zero-depth runs +-150 kb
$ST faidx $HASS CM056814.1:21128267-21133279 > $O/hass_kcs11_gene.fa # MRB53_020451..020452
declare -A ASM=(
  [thille_v4]=genome/thille.hifiv4.asm.hic.scaffold.masked.chr.only.fasta
  [hass_pri]=genome/hass_pseudohap/GCA_031208685.1/ncbi_dataset/data/GCA_031208685.1/GCA_031208685.1_ASM3120868v1_genomic.fna
  [hass_alt]=genome/hass_pseudohap/GCA_031208655.1/ncbi_dataset/data/GCA_031208655.1/GCA_031208655.1_ASM3120865v1_genomic.fna
  [puebla_hap1]=genome/puebla/GCA_060308345.1/ncbi_dataset/data/GCA_060308345.1/GCA_060308345.1_ASM6030834v1_genomic.fna
  [puebla_hap2]=genome/puebla/GCA_060308325.1/ncbi_dataset/data/GCA_060308325.1/GCA_060308325.1_ASM6030832v1_genomic.fna
  [lowland]=genome/wi_GCA_051132755.1/wi.chr.fasta
  [gwen_primary]=/epool/edwin/pamericana/avoflow/assemblies/thille.hic/pamer.contigs.c21.consensus.consensus_pilon_pilon_new_1000_0.5157_0.5426to2.4158_0.2004_primary.fasta
)
for n in "${!ASM[@]}"; do
  echo "== $n $(date)"
  $MM -c -x asm20 -t ${SLURM_CPUS_PER_TASK:-16} --secondary=yes -N 5 ${ASM[$n]} $O/hass_region.fa > $O/region.$n.paf
  $MM -c -x asm20 -t ${SLURM_CPUS_PER_TASK:-16} --secondary=yes -N 5 ${ASM[$n]} $O/hass_kcs11_gene.fa > $O/gene.$n.paf
done
echo "done $(date)"
