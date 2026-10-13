#!/bin/bash
# Mean/median/pooled alignment identity for the read-mapping BAMs (scripts/identity_stats.py).
set -uo pipefail
export PATH="/gpool/bin/samtools-1.23.1/bin:$PATH"
O=summaries/mapping_stats/rerun_20261002; A=/epool/edwin/pamericana/avoflow/mappingback; B=$A/rerun_20261002
run() { samtools view -@ 2 -F 0x904 $2 | python3 scripts/identity_stats.py $1 > $O/$1.identity.tsv.tmp && mv $O/$1.identity.tsv.tmp $O/$1.identity.tsv; echo "done $1 $(date)"; }
run thille_12chr $A/minimap2/thille.hifi.thille.bam &
run thille_13seq $B/thille.hifi.thille_chr_scf13.bam &
run hass_uoq_12chr $A/minimap2/thille.hifi.hass_uoq.bam &
run hass_uoq_full $B/thille.hifi.hass_uoq_full.bam &
run gwen_clr_thille_12chr $A/gwen.clr.thille.bam &
run gwen_clr_hass_uoq_12chr $A/gwen.clr.hass_uoq.bam &
wait; echo ALL_DONE
