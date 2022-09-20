#!/bin/bash
NAME=$1
REGION=$2
REFERENCE=$3
REF_HAP1=$4
REF_HAP2=$5
ASM=$6
ASM_HAP1=$7
ASM_HAP2=$8

if [[ "$ASM" == *.gfa ]]; then
    TEMP=$RANDOM.fa
    awk -f ./script/gfa_to_fa.awk "$ASM" >"$TEMP"
    ASM=$TEMP
fi

echo -ne "$REGION & $NAME"
if [ $# -ge 9 ]; then
    cargo run --release --bin compare_haplotypes -- "$REFERENCE" "$ASM" "$REF_HAP1" "$REF_HAP2" "$ASM_HAP1" "$ASM_HAP2" "$9"
else
    cargo run --release --bin compare_haplotypes -- "$REFERENCE" "$ASM" "$REF_HAP1" "$REF_HAP2" "$ASM_HAP1" "$ASM_HAP2"
fi
