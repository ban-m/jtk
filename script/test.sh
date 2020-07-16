#!/bin/bash
set -ue
cargo build --release 
INPUT=${PWD}/result/randseq/CCS_reads.1500.1M.entry.units.encode.clustered.new.json
cat ${INPUT} | \
    jtk polish_clustering -t 24 -vv 2> logfiles/log.polish |\
    jtk assemble -t 2 -vv > logfiles/assemble.gfa 2> logfiles/assemble.log
