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
qsub -q centos7.q -o ${output}/log \
     -S /bin/bash -N Workflow -cwd -pe smp 24 -j y -V \
     ./script/pipeline.sh \
     ${reads} ${output}/${readtype}

# for diploid in dip-A-CS dip-A-D dip-A-DS dip-A-C 
# do
#     output=${ROOT}/result/${diploid}/jtk
#     mkdir -p ${output}
#     reads=/data/hacone/randseq/synthetic_diploids/diploids/${diploid}.CLR.fa
#     if [ -e ${output}/${diploid}.CLR.log ]
#     then
#         rm -f ${output}/${diploid}.CLR.log
#     fi
#     qsub -q centos7.q -S /bin/bash -N Workflow -cwd -pe smp 24 -V \
#          -j y -o ${output}/${diploid}.CLR.log \
#          ./script/pipeline.sh \
#          ${reads}  ${output}/${diploid}.CLR
# done
