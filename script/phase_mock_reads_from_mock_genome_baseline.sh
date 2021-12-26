#!/bin/bash
PHASE=${PWD}/script/phasing_w_reference.sh
ROOT=${PWD}
set -ue

for coverage in 30x 50x
do
    READS=${ROOT}/data/mock_genome/haps.clr.${coverage}.fq
    REFERENCE=${ROOT}/data/mock_genome/hap_ref.fa
    OUTPUT=${ROOT}/result/mock_genome/clr.${coverage}/phasing
    rm -r ${OUTPUT}
    mkdir -p ${OUTPUT}
    qsub -q centos7.q -o ${OUTPUT}/log -j y \
         ${PHASE} ${READS} ${REFERENCE} ${OUTPUT}
done

for coverage in 15x 30x
do
    READS=${ROOT}/data/mock_genome/haps.ccs.${coverage}.fa
    REFERENCE=${ROOT}/data/mock_genome/hap_ref.fa
    OUTPUT=${ROOT}/result/mock_genome/ccs.${coverage}/phasing
    mkdir -p ${OUTPUT}
    qsub -q centos7.q -o ${OUTPUT}/log -j y \
         ${PHASE} ${READS} ${REFERENCE} ${OUTPUT}
done


