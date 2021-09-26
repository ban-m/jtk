#!/bin/bash
REFERENCE=~/work/sandbox/result/DBB_QBL_without_gaps.fasta
for UNIT in $@
do
    cat ${PWD}/units.fa | rg -A 1 ">${UNIT}$" | \
        minimap2 -p0.2 -c ${REFERENCE} - | awk '($2*9/10 < $4-$3)' 
done

