#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 12
#$ -o ./logfiles/debug.out
#$ -e ./logfiles/debug.log
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
ENCODED=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.json
cat ${ENCODED} | RUST_BACKTRACE=full ${JTK} clustering -vv --threads 12 --cluster_num 4 > ${CLUSTERED}

## ./target/release/debug ${ENCODED} 
