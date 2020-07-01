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
JTK=${PWD}/target/release/jtk
cargo build --release
for err in 300 600 1500 3000 6000 15000
do
    UNITS=${PWD}/result/CCS_reads.${err}.1M.entry.units.encode.clustered.json    
    cat ${UNITS} |\
        ${JTK} global_clustering -vv --threads 24 > ./logfiles/${err}.out 2> ./logfiles/${err}.log
done
