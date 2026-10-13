#!/bin/bash
# v4 copy (2026-10-05): v4 SyRI chain, genome label "Lowland accession". Four-genome plotsr ribbons (Lowland accession, Thille, Hass, Puebla; SyRI chain wi_ref -> hass -> hass_puebla) zoomed on the
# chr3 (R3) and chr10 (R10) SV loci. --reg is given on Thille (the chain's second genome) as the locus window used by
# the track figures; --rtr keeps all rearrangements within the homologous regions. Run from assemblies/thille.hic.
set -euo pipefail
SIF=/gpool/bin/singularityimages/graphgenomics.sif
O=compare/v4/regions/ribbons
singularity exec --no-home "$SIF" bash -c "
set -euo pipefail
export MPLCONFIGDIR=\$PWD/$O/.mpl
for r in R3:chr3:7598588-11696806 R10:chr10:10061230-15097080; do
  id=\${r%%:*}; reg=\${r#*:}
  for ext in pdf png; do
    plotsr --sr compare/v4/syri/wi_ref/wi_refsyri.out --sr compare/v4/syri/hass/hasssyri.out --sr compare/v4/syri/hass_puebla/hass_pueblasyri.out \
      --genomes $O/genomes.txt --reg Thille:\$reg --rtr -o $O/\${id}_fourway_ribbon.\$ext -s 10000 -H 6 -W 9
  done
done
"
