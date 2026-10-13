#!/bin/bash
#SBATCH --job-name=sdmyb.kmers
#SBATCH --partition=gcluster,gcbigmem
#SBATCH --cpus-per-task=4
#SBATCH --array=0-3
#SBATCH --output=logs/sdmyb_kmers.%A_%a.log
# SDMYB 21-mer test in the raw reads (scripts/sdmyb_kmer_test.py): rerun of summaries/sdmyb_kmers/ with the saved
# design (kmers.txt), plus Puebla HiFi reads (Edwin, 2026-10-11). One array task per read set; writes
# summaries/sdmyb_kmers_rerun/counts_<set>.txt (the original folder is not touched).
# Usage (repo root): sbatch scripts/sbatch_sdmyb_kmer_test.sh
#   then: python3 scripts/sdmyb_kmer_test.py summarize --counts-dir summaries/sdmyb_kmers_rerun \
#           --out summaries/sdmyb_kmers_rerun/allele_kmer_counts.tsv
set -euo pipefail
export PYTHONUNBUFFERED=1
NAMES=(Thille_hifi WI_hifi Gwen_clr Puebla_hifi)
READS=(/epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz
       /epool/edwin/pamericana/avoflow/rawdata/westindian/HiFi/pame.hifi.fq.gz
       /epool/edwin/pamericana/gwen/rawdata/pacbio/gwen/gwen_pacb_subreads.fastq.gz
       /epool/edwin/pamericana/avoflow/rawdata/hifi/puebla/SRR36605778_subreads.fastq.gz)
i=${SLURM_ARRAY_TASK_ID}
O=summaries/sdmyb_kmers_rerun; mkdir -p $O
out=$O/counts_${NAMES[$i]}.txt
[ -e "$out" ] && { echo "refusing to overwrite $out"; exit 1; }
echo "node $(hostname) task $i ${NAMES[$i]} ${READS[$i]} start $(date)"; grep --version | head -1
zcat "${READS[$i]}" | LC_ALL=C grep -o -F -f summaries/sdmyb_kmers/kmers.txt | sort | uniq -c > "$out.tmp"
mv "$out.tmp" "$out"
echo "done $(date): $(wc -l < "$out") distinct k-mers matched"
