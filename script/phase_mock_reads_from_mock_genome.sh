#!/bin/bash
set -ue
ROOT=${PWD}
for diploid in dip-A-CS dip-A-D dip-A-DS dip-A-C 
do
    output=${ROOT}/result/${diploid}/jtk
    mkdir -p ${output}
    reads=/data/hacone/randseq/synthetic_diploids/diploids/${diploid}.CLR.fa
    if [ -e ${output}/${diploid}.CLR.log ]
    then
        rm -f ${output}/${diploid}.CLR.log
    fi
    qsub -j y -o ${output}/${diploid}.CLR.log \
         ./script/pipeline_clr.sh \
         ${reads} \
         ${output}/${diploid}.CLR
done
