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
TARGET=/data/hacone/randseq/mut15000/CCS_reads.15000.1M.fa
ENCODE=${PWD}/result/CCS_reads.15000.1M.encode.units.json
UNITS=${PWD}/result/CCS_reads.15000.1M.units.fa
JTK=${PWD}/target/release/jtk
cargo build --release
cat ${TARGET} | ${JTK} entry | ${JTK} select_unit -vv \
    | ${JTK} stats -vv -f encode.log > ${ENCODE}
cat ${ENCODE} | ${JTK} extract -f fasta -t units > ${UNITS}


ALIGNMENT=${PWD}/result/CCS_reads.15000.1M.units.sam
minimap2 -a -t 12 -x map-pb ${UNITS} ${TARGET} > ${ALIGNMENT}

