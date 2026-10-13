#!/bin/bash
#SBATCH -J bwahic.thille                # jobname
#SBATCH -o bwahic.thille.o%A #.%a         # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bwahic.thille.e%A #.%a         # error file name A is the jobid and a is the arraytaskid
#SBATCH -c 24                    # start and stop of the array start-end
###SBATCH -A YOURACCOUNTID
###SBATCH --nodes=1
###SBATCH -N 1
#SBATCH -p gcbigmem,gcluster,debug           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --mem-per-cpu=5000
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

### bwa mem 0.7.19 & samtools 1.23.1
PATH=/gpool/bin/bwa-0.7.19:/gpool/bin/samtools-1.23.1:$PATH

REF=thille.hifiv3.asm.hic.p_ctg.fasta
REF=thille.hifiv3.asm.p_ctg.fasta
H1=DTG_OmniC_839_R1.fastq.gz
H2=DTG_OmniC_839_R2.fastq.gz

bwa index $REF
bwa mem -5SP -t $SLURM_CPUS_PER_TASK $REF $H1 $H2 | samtools view -bhS - > $(basename $REF .fasta)_omnic_mapped.bam
