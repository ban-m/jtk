#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 23
#$ -j y
#$ -m e 
#$ -V
set -ue
PATH="${PWD}/target/release/:${PATH}"
TARGET=$1
CLUSTERED=${2}.entry.units.encode.clustered.json
GFA=${2}.gfa
LOG=${2}.log
STAT=${2}.stat
RESULT=${2}.json
THREADS=23
jtk entry --input ${TARGET} |\
    jtk select_unit -vv --threads ${THREADS} |\
    jtk encode -vv --threads ${THREADS} |\
    jtk local_clustering -vv --threads ${THREADS} --cluster_num 2 |\
    tee ${CLUSTERED} |\
    jtk global_clustering -vv --threads 2 |\
    jtk stats -vv -f ${STAT} > ${RESULT} 
cat ${RESULT} | jtk assemble -t ${THREADS} -vv --output ${GFA} > /dev/null
