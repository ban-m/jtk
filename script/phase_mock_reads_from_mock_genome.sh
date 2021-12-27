#!/bin/bash
set -ue
ROOT=${PWD}

readtype=clr.30x
reads=${PWD}/data/mock_genome/haps.${readtype}.fa
output=${ROOT}/result/mock_genome/${readtype}/jtk
mkdir -p ${output}
if [ -e ${output}/log ]
then
    rm -f ${output}/log
fi
qsub -q centos7.q -j y -o ${output}/log \
     ./script/pipeline.sh \
     ${reads} \
     ${output}/${readtype}.CLR

# for diploid in dip-A-CS dip-A-D dip-A-DS dip-A-C 
# do
#     output=${ROOT}/result/${diploid}/jtk
#     mkdir -p ${output}
#     reads=/data/hacone/randseq/synthetic_diploids/diploids/${diploid}.CLR.fa
#     if [ -e ${output}/${diploid}.CLR.log ]
#     then
#         rm -f ${output}/${diploid}.CLR.log
#     fi
#     qsub -j y -o ${output}/${diploid}.CLR.log \
#          ./script/pipeline_clr.sh \
#          ${reads} \
#          ${output}/${diploid}.CLR
# done
