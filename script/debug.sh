#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 12
#$ -o ./logfiles/debug_chunk.out
#$ -e ./logfiles/debug_chunk.log
#$ -j y
#$ -m e 
#$ -V
set -ue 
TARGET=${PWD}/result/CCS_reads.15000.1M.fa
ENTRY=${PWD}/result/CCS_reads.15000.1M.entry.units.json
UNITS=${PWD}/result/CCS_reads.15000.1M.units.fa
JTK=${PWD}/target/release/jtk
CLUSTERED=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.clustered.json
ENCODED=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.json

cargo build --release

# cat ${ENCODED} | RUST_BACKTRACE=full ${JTK} clustering -vv --threads 12 --cluster_num 4 > ${CLUSTERED}

./target/release/debug ${ENCODED} 
