#!/bin/bash
OUTPUT=${PWD}/result/DBB_QBL/
mkdir -p ${OUTPUT}
for NAME in DBB_QBL_CCS_15Kbp_10x DBB_QBL_CCS_15Kbp_50x
do
    data=${PWD}/result/${NAME}.fastq
    cat ${data} | \
        paste - - - - |\
        cut -f1-2 |\
        tr '\t' '\n' |\
        sed -e "s/@/>/g" > ${OUTPUT}/${NAME}.fa
    data=${OUTPUT}/${NAME}.fa
    qsub -o ${OUTPUT}/${NAME}.pipeline.out -e ${OUTPUT}/${NAME}.pipeline.log\
         ./script/pipeline.sh \
         ${data} \
         ${OUTPUT}/${NAME} \
         6000
done

