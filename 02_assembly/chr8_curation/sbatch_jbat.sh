#!/bin/bash
#SBATCH -J jbat.thille
#SBATCH -o jbat.thille.o%j
#SBATCH -e jbat.thille.e%j
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -p gcbigmem,gcluster
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Manual-curation inputs for Juicebox Assembly Tools (chr8 = scaffold_12 + scaffold_13 join), following the YaHS 1.2.2
# README "Manual curation with Juicebox (JBAT)" recipe. Inputs are md5-verified copies of the original YaHS run
# (../yahs_thille.hifiv3.asm.hic.p_ctg, contigs thille.hifiv3.asm.hic.p_ctg.fasta), so the original run is untouched.
# --mem is set because the JVM (-Xmx32G) exceeds the 2.5 GB/CPU default (20 GB for 8 CPUs).
# Outputs: out_JBAT.{txt,liftover.agp,assembly,assembly.agp,log}, out_JBAT.hic (load .hic + .assembly in Juicebox).
# After curation: juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp thille.hifiv3.asm.hic.p_ctg.fasta
# Submit from this directory: sbatch sbatch_jbat.sh
set -euo pipefail

PATH=/gpool/bin/yahs:$PATH
JUICER_TOOLS=/gpool/bin/juicebox/juicer_tools.2.20.00.jar
CTG=thille.hifiv3.asm.hic.p_ctg.fasta

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-NA} mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID:-NA}"
juicer --version
java -version 2>&1 | head -1
free -g

if [[ -f out_JBAT.pre.done ]]; then
  echo "[$(date)] juicer pre -a outputs exist with .done, skipping"
else
  echo "[$(date)] juicer pre -a start"
  juicer pre -a -o out_JBAT yahs.out.bin yahs.out_scaffolds_final.agp "$CTG.fai" > out_JBAT.log 2>&1
  date > out_JBAT.pre.done
fi
grep -E "scale factor|JUICER_PRE CMD" out_JBAT.log || true

echo "[$(date)] juicer_tools pre start"
mkdir -p tmp_jbat_${SLURM_JOB_ID:-local}
# -n: job 12990 finished the matrix but juicer_tools 2.20.00 then crashed computing genome-wide norms
# (NullPointerException in GenomeWideNormalizationVectorUpdater); the partial file is out_JBAT.hic.part_pre_norm_npe.
java -Xmx32G -jar "$JUICER_TOOLS" pre -n -j "${SLURM_CPUS_PER_TASK:-8}" -t tmp_jbat_${SLURM_JOB_ID:-local} \
    out_JBAT.txt out_JBAT.hic.part <(grep PRE_C_SIZE out_JBAT.log | awk '{print $2" "$3}')
mv out_JBAT.hic.part out_JBAT.hic
echo "[$(date)] done: out_JBAT.hic + out_JBAT.assembly"
ls -la out_JBAT.*
