#!/bin/bash
#$ -S /bin/bash
#$ -N Workflow
#$ -cwd
#$ -pe smp 24
#$ -o ./logfiles/debug_chunk.out
#$ -e ./logfiles/debug_chunk.log
#$ -j y
#$ -m e 
#$ -V
set -ue 
UNITS=${PWD}/result/CCS_reads.15000.1M.entry.units.encode.clustered.json
JTK=${PWD}/target/release/jtk
cargo build --release 
cat ${UNITS} |\
    ${JTK} global_clustering -vv --threads 24 > ./logfiles/temp.log



