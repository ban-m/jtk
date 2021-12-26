#!/bin/bash
CANU=${PWD}/script/canu.sh
FLYE=${PWD}/script/flye.sh
HIFIASM=${PWD}/script/hifiasm.sh
ROOT=${PWD}
PHASE=${PWD}/script/phasing_w_reference.sh
set -ue
for NAME in COX_PGF_CLR_25Kbp_30x COX_PGF_CLR_50Kbp_30x COX_PGF_ONT_75Kbp_30x COX_PGF_CCS_25Kbp_30x
do
    output=${ROOT}/result/${NAME}/
    read=${PWD}/data/COX_PGF/${NAME}.fastq
    # Phasing
    GENOME="/grid/ban-m/human_genome/chm13.draft_v1.1.fasta"
    REGION="chr6:28573494-33272792"
    rm -fr ${output}/phasing
    mkdir -p ${output}/phasing
    REFERENCE=${output}/phasing/chm13_mhc.fasta
    samtools faidx ${GENOME} ${REGION} > ${REFERENCE}
    qsub -q centos7.q -o ${output}/phasing/log -j y \
         ${PHASE} ${read} ${REFERENCE} ${output}/phasing/
    
    # Canu 
    rm -fr ${output}/canu
    mkdir -p ${output}/canu
    bash ${CANU} ${read} 5M ${output}/canu 2> ${output}/canu/log
    
    # Flye
    rm -fr ${output}/flye
    mkdir -p ${output}/flye
    qsub -o ${output}/flye/log -j y ${FLYE} ${read} 5M ${output}/flye
    
    # FlyeHap
    # rm -fr ${output}/flyehap
    # mkdir -p ${output}/flyehap
    # qsub -q centos7.q -o ${output}/flyehap/log -j y ${PHASE} ${read} 2M ${output}/flyehap

    # Ra
    # rm -rf ${output}/ra
    # mkdir -p ${output}/ra
    # cd ${output}/ra/
    # ra -t 24 -x pb ${read} > ${output}/ra/${NAME}.fa
done

