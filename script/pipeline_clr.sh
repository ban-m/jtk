#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 23
#$ -j y
#$ -m e 
#$ -V
set -ue
PATH="${PATH}:${PWD}/target/release/"
TARGET=$1
CLUSTERED=${2}.entry.units.encode.clustered.json
RESULT=${2}.json
GFA=${2}.gfa
DRAFT_GFA=${2}.draft.gfa
LOG=${2}.log
STAT=${2}.stat
THREADS=23

### Take number.
if [ $# -gt 3 ];then
    echo "Unit guess is changed to" $3
    UNIT_GUESS=$3
else
    UNIT_GUESS=10000
fi

if [ -f ${2}.entry.json ]
then
    echo "Entry file found. Skip entry proc."
else
    jtk entry --input ${TARGET} --read_type CLR |\
        jtk mask_repeats -k 15 -t ${THREADS} -vv |\
        jtk select_unit -vv -t ${THREADS} --take_num ${UNIT_GUESS} |\
        jtk encode -vv --threads ${THREADS} >  ${2}.entry.json
fi
if [ -f ${CLUSTERED} ]
then
    echo "Clustered file found.Skip clustering proc."
else
    cat ${2}.entry.json |\
        jtk polish_encoding --threads ${THREADS} -vv |\
        jtk estimate_multiplicity -vv --threads ${THREADS} \
            --draft_assembly ${DRAFT_GFA} --max_cluster_size 6 |\
        jtk partition_local -vv --threads ${THREADS} >  ${CLUSTERED}
fi
if [ -f ${RESULT} ]
then
    echo "Global clustering seems to be done. Just assemble these files."
else    
    cat ${CLUSTERED} |\
        jtk correct_clustering -vv  --threads ${THREADS} |\
        # jtk partition_global -vv --threads ${THREADS} --graph |\
        jtk stats -vv -f ${STAT} > ${RESULT}
fi
cat ${RESULT} | jtk assemble -t ${THREADS} -vv --output ${GFA} --no_polish > /dev/null
cat ${GFA} | awk '($1 ~ /S/)' | awk 'BEGIN{OFS="\n"}{print ">" $2, $4}' > ${GFA%.gfa}.fa
