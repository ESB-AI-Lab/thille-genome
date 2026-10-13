#!/bin/bash
#SBATCH -J chr8.hifi.yahs
#SBATCH -o logs/chr8.hifi.yahs.o%j
#SBATCH -e logs/chr8.hifi.yahs.e%j
#SBATCH -c 32
#SBATCH --mem=100G
#SBATCH -p gcbigmem,gcluster,debug
#SBATCH --mail-user=youraddress@domain.org
#SBATCH --mail-type=end
#SBATCH --export=ALL
# Step 1 of the chr8 curation plan: Thille HiFi read depth over the scaffold_13 / chr8 (scaffold_12) overlap.
# The existing Thille HiFi BAMs were mapped to the 12-chromosome FASTA, which lacks scaffold_13, so reads from that
# arm would pile onto chr8. Here the reads are mapped to the complete YaHS assembly (all 1,491 scaffolds), with the
# same settings as the mappingback alignments (minimap2 2.30 -ax map-hifi --secondary=no), then depth is computed with
# samtools depth -a over scaffold_12 and scaffold_13. Half the median depth on both copies of the overlap would mark a
# haplotype duplication (false duplication) rather than real sequence.
# Submit from assemblies/thille.hic: sbatch curation/chr8/sbatch_map_hifi_yahs.sh
set -euo pipefail
export PYTHONUNBUFFERED=1
MINIMAP2=/gpool/bin/minimap2-2.30_x64-linux/minimap2
SAMTOOLS=/gpool/bin/samtools-1.24/samtools
REF=/epool/edwin/pamericana/avoflow/hifiasm/thille/hifi.hic.gcluster/yahs_thille.hifiv3.asm.hic.p_ctg/yahs.out_scaffolds_final.fa
READS=/epool/edwin/pamericana/avoflow/rawdata/hifi/thille/thille.hifiv3.fastq.gz
OUT=curation/chr8
BAM=$OUT/thille.hifi.yahs_final.bam
T=${SLURM_CPUS_PER_TASK}
echo "node=$(hostname) cpus=$T mem=${SLURM_MEM_PER_NODE:-NA}MB job=${SLURM_JOB_ID}"
echo "minimap2 $($MINIMAP2 --version)  $($SAMTOOLS --version | head -1)"
free -g
if [[ ! -s $BAM.done ]]; then
  echo "[$(date)] minimap2 -ax map-hifi start"
  $MINIMAP2 -ax map-hifi --secondary=no -t "$T" "$REF" "$READS" | $SAMTOOLS sort -@ 8 -m 3G -o "$BAM.tmp.bam" -
  mv "$BAM.tmp.bam" "$BAM"; $SAMTOOLS index -@ 8 "$BAM"; date > "$BAM.done"
fi
for s in scaffold_12 scaffold_13; do
  [[ -s $OUT/depth.$s.tsv.gz ]] || $SAMTOOLS depth -a -r "$s" "$BAM" | gzip > "$OUT/depth.$s.tsv.gz"
done
echo "[$(date)] done"
