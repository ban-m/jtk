#!/bin/bash
PHASE=${PWD}/script/phasing_w_reference.sh
FLYE=${PWD}/script/flye.sh
FLYE_HIFI=${PWD}/script/flye_hifi.sh
ROOT=${PWD}
set -ue
for read_type in clr.30x clr.50x 
do
    # with-reference
    READS=${ROOT}/data/mock_genome/haps.${read_type}.fq
    REFERENCE=${ROOT}/data/mock_genome/hap_ref.fa
    OUT_ROOT=${ROOT}/result/mock_genome/${read_type}
    OUTPUT=${OUT_ROOT}/phasing
    rm -rf ${OUTPUT}
    mkdir -p ${OUTPUT}
    qsub -q centos7.q -o ${OUTPUT}/log -j y \
         ${PHASE} ${READS} ${REFERENCE} ${OUTPUT}
    # Flye
    OUTPUT=${OUT_ROOT}/flye
    rm -rf ${OUTPUT}
    mkdir -p ${OUTPUT}
    qsub -o ${OUTPUT}/log -j y ${FLYE} ${READS} 2M $OUTPUT
done


for read_type in ccs.15x ccs.30x
do
    # with-reference
    READS=${ROOT}/data/mock_genome/haps.${read_type}.fq
    REFERENCE=${ROOT}/data/mock_genome/hap_ref.fa
    OUT_ROOT=${ROOT}/result/mock_genome/${read_type}
    OUTPUT=${OUT_ROOT}/phasing
    rm -rf ${OUTPUT}
    mkdir -p ${OUTPUT}
    qsub -q centos7.q -o ${OUTPUT}/log -j y \
         ${PHASE} ${READS} ${REFERENCE} ${OUTPUT}
    # Flye
    OUTPUT=${OUT_ROOT}/flye
    rm -rf ${OUTPUT}
    mkdir -p ${OUTPUT}
    qsub -o ${OUTPUT}/log -j y ${FLYE_HIFI} ${READS} 2M $OUTPUT
done
