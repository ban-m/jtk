#!/bin/bash
TARGET=/data/hacone/randseq/mut15000/CCS_reads.15000.1M.fa
ENCODE=${PWD}/result/CCS_reads.15000.1M.encode.json
JTK=${PWD}/target/release/jtk
cargo build --release
cat ${TARGET} | ${JTK} entry | ${JTK} stats -f encode.log > ${ENCODE}
