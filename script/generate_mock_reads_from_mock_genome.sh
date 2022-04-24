#!/bin/bash
set -ue 
OUTPUT=${PWD}/data/mock_genome/
mkdir -p ${OUTPUT}
ONT_READS=${OUTPUT}/haps.ont
HAPA=${OUTPUT}/hapA.fa
HAPB=${OUTPUT}/hapB.fa
REFERENCE=${OUTPUT}/hap_ref.fa
MERGED=${OUTPUT}/haps.fa
SEED=3490283
SEED_BR=3910
RENAME=${PWD}/script/rename_fastq.awk


# Gen genome
# cargo run --release --bin gen_sim_genome -- ${REFERENCE} ${HAPA} ${HAPB} ${SEED} 
# cat ${HAPA} ${HAPB} > ${MERGED}

# simulator.py genome --seed 100 -n 10000 -rg ${MERGED} -o ${ONT_READS}.30x -t 56 -c ${PWD}/../../tools/nanosim/training


OUTPUT=${PWD}/data/mock_genome_segdup/
mkdir -p ${OUTPUT}
HAPS=${OUTPUT}/haps.fa
ONT_READS=${OUTPUT}/haps.ont
SEED=3490282393
cargo run --release --bin gen_sim_genome_segdup -- ${HAPS} ${SEED} 
simulator.py genome --seed 100 -n 23000 -rg ${HAPS} -o ${ONT_READS}.30x -t 56 -c ${PWD}/../../tools/nanosim/training
