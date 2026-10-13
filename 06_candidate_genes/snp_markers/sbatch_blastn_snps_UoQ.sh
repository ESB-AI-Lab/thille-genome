#!/bin/bash
#SBATCH -J blastn_snps_UoQ       # Job name
#SBATCH -o blastn_snps_UoQ.o%j   # Output file name (%j expands to jobId)
#SBATCH -e blastn_snps_UoQ.e%j   # Error file name
#SBATCH -A YOURACCOUNTID            # Project code
#SBATCH --nodes=1                # Number of nodes
#SBATCH -p RM-shared             # Queue (partition)
#SBATCH -t 20:00:00              # Runtime (hh:mm:ss)
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=ALL          # Email notifications

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

# Map the extracted flowering SNP markers onto UoQ genome

# Input files
SHRSPa_FASTA=SHRSPa_sequence.fasta
FLOWERING_SNPS_TXT=vanessafloweringsnps.new.txt
UOQ_GENOME=GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta

# Output files
FLOWERING_MARKERS_FASTA=flowering_snp_markers.fasta
UOQ_OUT=vanessaMSL_UoQ_blastn_5.outfmt6

# Extract flowering SNP markers' IDs from the TXT file and Extract sequences from SHRSPa_sequence.fasta based on these IDs
cut -f1 $FLOWERING_SNPS_TXT > flowering_snp_ids.txt
awk 'NR==FNR{ids[$1]; next} /^>/{f=substr($1,2) in ids} f' flowering_snp_ids.txt $SHRSPa_FASTA > $FLOWERING_MARKERS_FASTA
echo "Extracted flowering SNP markers into $FLOWERING_MARKERS_FASTA"

# Load BLAST module
module load BLAST
makeblastdb -in $UOQ_GENOME -input_type fasta -dbtype nucl -parse_seqids -out uoq_db
blastn -query $FLOWERING_MARKERS_FASTA -db uoq_db -out $UOQ_OUT -outfmt 6 -max_target_seqs 10 -evalue 1e-5 -num_threads 4
echo "BLASTn against UoQ genome completed."
