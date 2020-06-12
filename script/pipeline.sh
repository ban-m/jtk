#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 24
#$ -o ./logfiles/pipeline.out
#$ -e ./logfiles/pipeline.log
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
LOG=${2}.log

cargo build --release
cat ${TARGET} | ${JTK} entry |\
    ${JTK} select_unit -vv |\
    tee ${ENTRY} |\
    ${JTK} extract -f fasta -t units > ${UNITS}

cd ${PWD}/result
REFNAME=${RANDOM}
echo ${REFNAME} 1>&2
lastdb -R00 -Q0 ${REFNAME} ${UNITS}
last-train -P12 -Q0 ${REFNAME} ${TARGET} > ${REFNAME}.matrix
lastal -f maf -P12 -R00 -Q0 -p ${REFNAME}.matrix ${REFNAME} ${TARGET}|\
    maf-convert tab --join 500 > ${REFNAME}_alignments.tab
cd ../


cat ${ENTRY} |\
    ${JTK} encode -vv -a ${PWD}/result/${REFNAME}_alignments.tab |\
    ${JTK} stats -vv -f ${LOG} |\
    ${JTK} clustering -vv --threads 24 --cluster_num 3 > ${CLUSTERED}
# ./target/release/debug ${ENCODED} 
