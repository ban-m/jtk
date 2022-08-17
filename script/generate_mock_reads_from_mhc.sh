#!/bin/bash
# Require: Nanosim
set -ue
REFERENCE=${PWD}/data/COX_PGF/COX_PGF.fa
REANME=${PWD}/script/rename_fastq.awk
DATA=${PWD}/data/COX_PGF/
mkdir -p ${DATA}
TARGET=${DATA}/COX_PGF_ONT_30x
simulator.py genome --seed 100 -n 50000 -rg ${REFERENCE} -o ${TARGET} -t 56 -c ${PWD}/../../tools/nanosim/training