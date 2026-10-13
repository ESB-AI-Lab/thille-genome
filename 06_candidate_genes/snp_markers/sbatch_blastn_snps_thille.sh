#!/bin/bash
#SBATCH -J blastn_snps_thille       # Job name
#SBATCH -o blastn_snps_thille.o%j   # Output file name (%j expands to jobId)
#SBATCH -e blastn_snps_thille.e%j   # Error file name
#SBATCH -A YOURACCOUNTID            # Project code
#SBATCH --nodes=1                # Number of nodes
#SBATCH -p RM-shared             # Queue (partition)
#SBATCH -t 20:00:00              # Runtime (hh:mm:ss)
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=ALL          # Email notifications

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

# Map the extracted flowering SNP markers onto Thille genome

# Input files
SHRSPa_FASTA=SHRSPa_sequence.fasta
FLOWERING_SNPS_TXT=vanessafloweringsnps.new.txt
THILLE_GENOME=thille.hifiv4.asm.hic.scaffold.chr.only.fasta

# Output files
FLOWERING_MARKERS_FASTA=flowering_snp_markers.fasta
THILLE_OUT=vanessaMSL_thille_blastn_5.outfmt6

# Extract flowering SNP markers' IDs from the TXT file and Extract sequences from SHRSPa_sequence.fasta based on these IDs
cut -f1 $FLOWERING_SNPS_TXT > flowering_snp_ids.txt
awk 'NR==FNR{ids[$1]; next} /^>/{f=substr($1,2) in ids} f' flowering_snp_ids.txt $SHRSPa_FASTA > $FLOWERING_MARKERS_FASTA
echo "Extracted flowering SNP markers into $FLOWERING_MARKERS_FASTA"

# Load BLAST module
module load BLAST
makeblastdb -in $THILLE_GENOME -input_type fasta -dbtype nucl -parse_seqids -out thille_db
blastn -query $FLOWERING_MARKERS_FASTA -db thille_db -out $THILLE_OUT -outfmt 6 -max_target_seqs 10 -evalue 1e-5 -num_threads 4
echo "BLASTn against Thille genome completed."
