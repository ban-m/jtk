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
RESULT=${2}.json
GFA=${2}.gfa
DRAFT_GFA=${2}.draft.gfa
LOG=${2}.log
STAT=${2}.stat
THREADS=23

cat ${TARGET} |\
    jtk entry |\
    jtk select_unit -vv -t ${THREADS} --read_type CLR |\
    jtk encode -vv --threads ${THREADS} |\
    jtk multiplicity_estimation -vv --threads ${THREADS} \
        --draft_assembly ${DRAFT_GFA} --max_multiplicity 3 |\
    jtk local_clustering -vv --threads ${THREADS} --read_type CLR|\
    tee ${CLUSTERED} |\
    jtk clustering_correction -vv --threads ${THREADS} |\
    jtk local_clustering -vv --threads ${THREADS} --read_type CLR \
        --retain_current_clustering |\
    jtk clustering_correction -vv --threads ${THREADS} |\
    jtk global_clustering -vv --threads ${THREADS} |\
    jtk stats -vv -f ${STAT} |\
    jtk assemble -t2 -vv --output ${GFA} > ${RESULT}
