#!/bin/bash
CANU=${PWD}/script/canu.sh
FLYE=${PWD}/script/flye.sh
HIFIASM=${PWD}/script/hifiasm.sh

ROOT=${PWD}
set -ue
for mut_rate in 300 600 1500 3000 6000 15000
do
    read=/data/hacone/randseq/mut${mut_rate}/CCS_reads.${mut_rate}.1M.fa
    output=${ROOT}/result/mut${mut_rate}
    mkdir -p ${output}
    mkdir -p ${output}/flye
    mkdir -p ${output}/canu
    mkdir -p ${output}/hifiasm
    bash ${CANU} ${read} 1M ${output}/canu 2> ${output}/canu/log
    qsub -o ${output}/flye/log -j y ${FLYE} ${read} 1M ${output}/flye
    cd ${output}/hifiasm
    qsub -o ${output}/hifiasm -j y ${HIFIASM} ${read} 1M hifiasm
    cd ../
done
