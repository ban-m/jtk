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
TARGET=${PWD}/result/CCS_reads.15000.1M.fa
ENTRY=${PWD}/result/CCS_reads.15000.1M.entry.units.json
UNITS=${PWD}/result/CCS_reads.15000.1M.units.fa
JTK=${PWD}/target/release/jtk
CLUSTERED=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.clustered.json


cargo build --release
cat ${TARGET} | ${JTK} entry |\
    ${JTK} select_unit -vv |\
    tee ${ENTRY} |\
    ${JTK} extract -f fasta -t units > ${UNITS}

cd ${PWD}/result
lastdb -R00 -Q0 units ${UNITS}
last-train -P12 -Q0 units ${TARGET} > score.matrix
lastal -f maf -P12 -R00 -Q0 -p score.matrix units ${TARGET}|\
    maf-convert tab --join 500 > alignments.tab
cd ../

ENCODED=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.json
cat ${ENTRY} |\
    ${JTK} encode -vv -a ${PWD}/result/alignments.tab |\
    ${JTK} stats -vv -f ${PWD}/result/encode.log > ${ENCODED}
cat ${ENCODED} | ${JTK} clustering -vv --threads 24 --cluster_num 4 > ${CLUSTERED}
