#!/bin/bash
set -ue
OUTPUT=${PWD}/result/
mkdir -p ${OUTPUT}
# for FASTA in ${PWD}/data/COX_PGF/*.fa
for FASTA in ${PWD}/data/COX_PGF/COX_PGF_CLR_25Kbp_30x.fa
do
    NAME=${FASTA#${PWD}/data/COX_PGF/}
    NAME=${NAME%.fa}
    mkdir -p ${OUTPUT}/${NAME}/jtk
    if [ -f ${OUTPUT}/${NAME}/jtk/${NAME}.out ]
    then
        echo "Log file found. Moving."
        mv ${OUTPUT}/${NAME}/jtk/${NAME}.out ${OUTPUT}/${NAME}/jtk/${NAME}.old.${RANDOM}.out
    fi
    qsub -q centos7.q -o ${OUTPUT}/${NAME}/jtk/${NAME}.out -j y\
         -V -S /bin/bash -cwd -pe smp 24 -V \
         ./script/pipeline.sh \
         ${FASTA} \
         ${OUTPUT}/${NAME}/jtk/${NAME}
done

