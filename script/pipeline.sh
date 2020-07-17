#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 23
# -o ./logfiles/pipeline.out
# -e ./logfiles/pipeline.log
#$ -j y
#$ -m e 
#$ -V
set -ue
JTK=${PWD}/target/release/jtk
TARGET=$1
ENTRY=${2}.entry.units.json
UNITS=${2}.units.fa
ENCODED=${2}.entry.units.encode.json
CLUSTERED=${2}.entry.units.encode.clustered.json
GFA=${2}.gfa
LOG=${2}.log

if [ $# -gt 2 ];then
CHUNK_NUM=$3
else
CHUNK_NUM=2000
fi

# cargo build --release
cat ${TARGET} | ${JTK} entry |\
    ${JTK} select_unit -vv -t 23 --chunk_num ${CHUNK_NUM} |\
    tee ${ENTRY} |\
    ${JTK} extract -t units > ${UNITS}

cd ${PWD}/result
REFNAME=${RANDOM}
echo ${REFNAME} 1>&2
lastdb -R00 -Q0 ${REFNAME} ${UNITS}
last-train -P23 -Q0 ${REFNAME} ${TARGET} > ${REFNAME}.matrix
lastal -f maf -P23 -R00 -Q0 -p ${REFNAME}.matrix ${REFNAME} ${TARGET}|\
    maf-convert tab --join 500 > ${2}.alignment.tab
rm ${REFNAME}*
cd ../

cat ${ENTRY} |\
    ${JTK} encode -vv -a ${2}.alignment.tab |\
    ${JTK} local_clustering -vv --threads 23 --cluster_num 3 |\
    ${JTK} stats -vv -f ${LOG} |\
    ${JTK} global_clustering -vv --threads 2 > ${CLUSTERED}
cat ${CLUSTERED} | \
    ${JTK} polish_clustering -t23 -vv |\
    ${JTK} assemble -t2 -vv > ${GFA}
