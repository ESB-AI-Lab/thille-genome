#!/bin/bash
#SBATCH -J thille.hic.seqkit
#SBATCH -o logs/thille.hic.seqkit.o%j
#SBATCH -e logs/thille.hic.seqkit.e%j
#SBATCH -c 4
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --export=ALL
# Read stats for the Thille Dovetail Omni-C reads (Table 1A): reads, total bases, length distribution, Q20/Q30, GC.
# Submit from assemblies/thille.hic: sbatch stats/sbatch_hic_seqkit.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} job=${SLURM_JOB_ID}"
D=/epool/edwin/pamericana/avoflow/rawdata/dovetail-hic/thille
/gpool/bin/seqkit/seqkit version 2>/dev/null || /gpool/bin/seqkit version
SK=$( [[ -x /gpool/bin/seqkit/seqkit ]] && echo /gpool/bin/seqkit/seqkit || echo /gpool/bin/seqkit )
$SK stats -a -T -j 2 $D/DTG_OmniC_839_R1.fastq.gz $D/DTG_OmniC_839_R2.fastq.gz > stats/thille_hic_seqkit.tsv.tmp
mv stats/thille_hic_seqkit.tsv.tmp stats/thille_hic_seqkit.tsv
cat stats/thille_hic_seqkit.tsv
