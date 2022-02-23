#!/bin/bash
REFERENCE=${PWD}/data/COX_PGF.fa
for UNIT in $@
do
    cat ${PWD}/units.fa | rg -A 1 ">${UNIT}$" | \
        minimap2 --eqx -p0.2 -c ${REFERENCE} - | awk '($2*9/10 < $4-$3)' 
done

