#!/bin/bash
REFERENCE=~/work/sandbox/result/DBB_QBL_without_gaps.fasta
for UNIT in $@
do
    cat ${PWD}/units.fa | rg -A 1 ">${UNIT}" | \
        minimap2 -c ${REFERENCE} - | awk '(1900 < $4-$3)' 
done

