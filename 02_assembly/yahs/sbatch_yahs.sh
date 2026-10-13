#!/bin/bash
#SBATCH -J yahshic.thille                # jobname
#SBATCH -o yahshic.thille.o%A #.%a         # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e yahshic.thille.e%A #.%a         # error file name A is the jobid and a is the arraytaskid
#SBATCH -c 8                    # start and stop of the array start-end
###SBATCH -A YOURACCOUNTID
###SBATCH --nodes=1
###SBATCH -N 1
#SBATCH -p gcbigmem,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --mem-per-cpu=6000
###SBATCH -t 48:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
#export MV2_SHOW_CPU_BINDING=1
#  $SLURM_CPUS_ON_NODE
#  $SLURM_CPUS_PER_TASK

###yahs 6c46061
PATH=/gpool/bin/yahs:$PATH

REF=thille.hifiv3.asm.hic.p_ctg.fasta
REF=thille.hifiv3.asm.p_ctg.fasta
BAM=$(basename $REF .fasta)_omnic_mapped.bam

yahs $REF $BAM
juicer pre yahs.out.bin yahs.out_scaffolds_final.agp ${REF}.fai > ${REF}_out_jbmap.txt
