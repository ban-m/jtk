#!/bin/bash
for NAME in DBB_QBL_CCS_15Kbp_10x
do
    data=~/work/hla_haplotyper/result/${NAME}.fastq
    cat ${data} | \
        paste - - - - |\
        cut -f1-2 |\
        tr '\t' '\n' |\
        sed -e "s/@/>/g" > ${PWD}/result/${NAME}.fa
    data=${PWD}/result/${NAME}.fa
    qsub -o ./logfiles/${NAME}.out -e ./logfiles/${NAME}.log\
         ./script/pipeline.sh \
         ${data} \
         ${PWD}/result/${NAME}
done

