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
mkdir -p ${2}
if [ -f ${2}.entry.json ]
then
    echo "Entry file found. Skip entry proc."
else
    jtk entry --input ${TARGET} --read_type CLR |\
        jtk repeat_masking -k 15 -t ${THREADS} -vv |\
        jtk select_unit -vv -t ${THREADS} --take_num 10000 |\
        jtk encode -vv --threads ${THREADS} >  ${2}.entry.json
fi
if [ -f ${CLUSTERED} ]
then
    echo "Clustered file found.Skip clustering proc."
else
    cat ${2}.entry.json |\
        jtk multiplicity_estimation -vv --threads ${THREADS} \
            --draft_assembly ${DRAFT_GFA} --max_cluster_size 6 |\
        jtk local_clustering -vv --threads ${THREADS} >  ${CLUSTERED}
fi
if [ -f ${RESULT} ]
then
    echo "Global clustering seems to be done. Just assemble these files."
else    
    cat ${CLUSTERED} |\
        jtk clustering_correction -vv --graph --threads ${THREADS} |\
        jtk global_clustering -vv --threads ${THREADS} |\
        jtk stats -vv -f ${STAT} > ${RESULT}
fi
cat ${RESULT} | jtk assemble -t ${THREADS} -vv --output ${GFA} > /dev/null
cat ${GFA} | awk '($1 ~ /S/)' | awk 'BEGIN{OFS="\n"}{print ">" $2, $4}' > ${GFA%.gfa}.fa
