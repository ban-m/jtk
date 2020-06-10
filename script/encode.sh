#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 12
#$ -o ./logfiles/encode.log
#$ -j y
#$ -m e 
#$ -V
set -ue 
TARGET=${PWD}/result/CCS_reads.15000.1M.fa
ENTRY=${PWD}/result/CCS_reads.15000.1M.entry.units.json
UNITS=${PWD}/result/CCS_reads.15000.1M.units.fa
JTK=${PWD}/target/release/jtk
cargo build --release
cat ${TARGET} | ${JTK} entry |\
    ${JTK} select_unit -vv |\
    ${JTK} stats -vv -f encode.log |\
    tee ${ENTRY} |\
    ${JTK} extract -f fasta -t units > ${UNITS}

# ALIGNMENT=${PWD}/result/CCS_reads.15000.1M.units.sam
# minimap2 -a -t 12 -x map-pb ${UNITS} ${TARGET} > ${ALIGNMENT}
cd ${PWD}/result
lastdb -R00 -Q0 units ${UNITS}
last-train -P12 -Q0 units ${TARGET} > score.matrix
lastal -f maf -P12 -R00 -Q0 -p score.matrix units ${TARGET}|\
    maf-convert tab --join 500 > alignments.tab
cd ../

ENCODE=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.json
cat ${ENTRY} |\
    ${JTK} encode -vv -a ${PWD}/result/alignments.tab > ${ENCODE}

