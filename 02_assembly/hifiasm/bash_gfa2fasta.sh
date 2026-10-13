#!/bin/bash

PATH=/gpool/bin/bioawk:$PATH

i=$1
awk '/^S/{print ">"$2"\n"$3}' ${i} > $(basename ${i} .gfa).fasta
bioawk -cfastx '{print $name"\t"length($seq)}' $(basename ${i} .gfa).fasta | sort -rnk2,2 > $(basename ${i} .gfa).headers
