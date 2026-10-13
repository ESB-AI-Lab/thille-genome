#!/bin/bash
#SBATCH -J hifiasmhic.thille                # jobname
#SBATCH -o hifiasmhic.thille.o%A #.%a         # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e hifiasmhic.thille.e%A #.%a         # error file name A is the jobid and a is the arraytaskid
#SBATCH -c 24                    # start and stop of the array start-end
###SBATCH -A YOURACCOUNTID
###SBATCH --nodes=1
###SBATCH -N 1
#SBATCH -p gcbigmem,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --mem-per-cpu=20000
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

###hifiasm 0.25.1
PATH=/gpool/bin/hifiasm-0.25.0:$PATH

#FASTQ=thille.hifi.fastq.gz
FASTQ=thille.hifiv3.fastq.gz
H1=DTG_OmniC_839_R1.fastq.gz
H2=DTG_OmniC_839_R2.fastq.gz

hifiasm -o $(basename ${FASTQ} .fastq.gz).asm -t ${SLURM_CPUS_PER_TASK} --h1 $H1 --h2 $H2 ${FASTQ}
#hifiasm -o $(basename ${FASTQ} .fastq.gz).asm --primary -t ${SLURM_CPUS_PER_TASK} ${FASTQ}

#bash bash_gfa2fasta.sh $(basename ${FASTQ} .fastq.gz).asm.p_ctg.gfa
#bash bash_gfa2fasta.sh $(basename ${FASTQ} .fastq.gz).asm.a_ctg.gfa
