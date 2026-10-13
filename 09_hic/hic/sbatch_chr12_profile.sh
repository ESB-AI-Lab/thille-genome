#!/bin/bash
#SBATCH -J chr12_hic_profile
#SBATCH -p gcbigmem,gcluster
#SBATCH -c 1
#SBATCH -o hic/v4/chr12_start/chr12_profile.%j.log
# Where do the Hi-C contacts of each 1 Mb bin of Thille v4 chr12 (0-20 Mb) go? Partner classes: chr12 0-8 Mb,
# chr12 8-20 Mb, chr12 >= 20 Mb, other chromosomes (by name). Same pairs as the Figure 4 maps (juicer pre, MAPQ >= 10).
# Submit from assemblies/thille.hic.
set -euo pipefail
echo "node=$(hostname) job=$SLURM_JOB_ID"
awk -F'\t' '
function cls(c,p){ if(c!="chr12") return "other:"c; if(p<8000000) return "chr12_0-8"; if(p<20000000) return "chr12_8-20"; return "chr12_ge20" }
function add(c,p,c2,p2){ if(c=="chr12" && p<20000000){ b=int(p/1000000); k=cls(c2,p2); n[b,k]++; tot[b]++; keys[k]=1 } }
{ add($2,$3,$6,$7); if(!($2==$6 && $3==$7)) add($6,$7,$2,$3) }
END{ for(b=0;b<20;b++) for(k in keys) if(n[b,k]) printf "%d\t%s\t%d\t%d\n", b, k, n[b,k], tot[b] }' hic/v4/all/all.pairs.sorted.txt | sort -k1,1n -k3,3nr > hic/v4/chr12_start/chr12_bin_partners.tsv
echo done
