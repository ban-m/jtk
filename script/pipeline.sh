#!/bin/bash
set -ue
PATH="${PATH}:${PWD}/target/release/"
TARGET=$1
CLUSTERED=${2}.entry.units.encode.clustered.json
PURGED=${2}.entry.units.encode.clustered.purged.json
RESOLVED=${2}.entry.units.encode.clustered.purged.resolved.json
RESULT=${2}.json
GFA=${2}.gfa
DRAFT_GFA=${2}.draft.gfa
DRAFT_GFA_2=${2}.draft2.gfa
LOG=${2}.log
STAT=${2}.stat
READTYPE=${3}
THREADS=32

source /bio/package/gcc/setup8.sh

### Take number.
if [ $# -ge 4 ]; then
    echo "Unit guess is changed to" $4
    UNIT_GUESS=$4
else
    UNIT_GUESS=500
fi

if [ -f ${2}.entry.json ]
then
    echo "Entry file found. Skip entry proc."
else
    jtk entry --input ${TARGET} --read_type $READTYPE |\
        jtk mask_repeats -k 15 -t ${THREADS} -vv |\
        jtk select_unit -vv -t ${THREADS} --take_num ${UNIT_GUESS} |\
        jtk pick_components -vv -c1 -t${THREADS}|\
        jtk correct_deletion -vv --threads ${THREADS} > ${2}.entry.json
fi


UPPER_COPY_NUM=6
if [ -f ${CLUSTERED} ]
then
    echo "Clustered file found.Skip clustering proc."
else
    cat ${2}.entry.json |\
        jtk polish_encoding --threads ${THREADS} -vv |\
        jtk estimate_multiplicity -vv --threads ${THREADS} --draft_assembly ${DRAFT_GFA} --purge_copy_num ${UPPER_COPY_NUM} |\
        tee ${2}.entry.units.encode.json |\
        jtk partition_local -vv --threads ${THREADS} >  ${CLUSTERED}
fi


if [ -f ${PURGED} ]
then
    echo "Suspicious encodings are already purged."
else
    cat ${CLUSTERED} | jtk purge_diverged --threads ${THREADS} -vv > ${PURGED}
fi

if [ -f ${RESOLVED} ]
then
    echo "Tangle resolved. Skip resolving proc".
else
    cat ${PURGED} |\
        jtk correct_deletion -vv --threads ${THREADS} --re_cluster |\
        jtk encode_densely -vv --threads ${THREADS} |\
        jtk correct_deletion -vv --threads ${THREADS} --re_cluster >${RESOLVED}
fi

if [ -f ${RESULT} ]
then
    echo "Global clustering seems to be done. Just assemble these files."
else
    cat ${RESOLVED} |\
        # jtk correct_clustering -vv  --threads ${THREADS} |\
        jtk stats -vv -f ${STAT} > ${RESULT}
fi

cat ${RESULT} | jtk assemble -t ${THREADS} -vv --output ${GFA} --no_polish > /dev/null
cat ${GFA} | awk '($1 ~ /S/)' | awk 'BEGIN{OFS="\n"}{print ">" $2, $4}' > ${GFA%.gfa}.fa
