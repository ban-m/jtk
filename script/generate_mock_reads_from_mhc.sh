#!/bin/bash
set -ue
# Update 2021/10/28 DBB_QBL->COX_PGF because of the previous study.
# REFERENCE=~/work/sandbox/result/DBB_QBL_without_gaps.fasta
REFERENCE=${PWD}/data/COX_PGF.fa
REANME=${PWD}/script/rename_fastq.awk
DATA=${PWD}/data/COX_PGF/
mkdir -p ${DATA}

TARGET=${DATA}/COX_PGF_CLR_25Kbp_30x
badread simulate --reference ${REFERENCE} --quantity 30x \
        --error_model pacbio2016 --qscore_model pacbio2016 \
        --length 25000,2000 \
        --identity 85,95,5 \
        --seed 10 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" |\
    awk -f ${REANME} > ${TARGET}.fastq
cat ${TARGET}.fastq | paste - - - - |\
    cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${TARGET}.fa


TARGET=${DATA}/COX_PGF_CLR_50Kbp_30x
badread simulate --reference ${REFERENCE} --quantity 30x \
        --error_model pacbio2016 --qscore_model pacbio2016 \
        --length 50000,2000 \
        --identity 85,95,5 \
        --seed 10 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" |\
    awk -f ${REANME} > ${TARGET}.fastq
cat ${TARGET}.fastq | paste - - - - |\
    cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${TARGET}.fa


TARGET=${DATA}/COX_PGF_CCS_25Kbp_30x
badread simulate --reference ${REFERENCE} --quantity 30x \
        --error_model pacbio2016 --qscore_model pacbio2016 \
        --length 25000,2000 \
        --identity 99.9,100,0.2 \
        --seed 10 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" |\
    awk -f ${REANME} > ${TARGET}.fastq
cat ${TARGET}.fastq | paste - - - - |\
    cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${TARGET}.fa

TARGET=${DATA}/COX_PGF_ONT_75Kbp_30x
badread simulate --reference ${REFERENCE} --quantity 30x \
        --error_model nanopore2020 --qscore_model nanopore2020 \
        --length 75000,5000 \
        --identity 90,95,3 \
        --seed 100 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" |\
    awk -f ${REANME} > ${TARGET}.fastq
cat ${TARGET}.fastq | paste - - - - |\
    cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${TARGET}.fa


