#!/bin/bash
# Synopsis: `bash ./script/benchmark_mcmc.sh $OUTPUT_DIR $THREADS`
BENCH="$PWD"/target/release/benchmark_mcmc
THREADS=${2:-1}
mkdir -p "$1"
NUM_LOOP=250
LOOP=$((NUM_LOOP / THREADS))
cargo build --release

RESULT="$1"/mcmc_compare.tsv
rm -f "$RESULT"
touch "$RESULT"
for cov in 20; do
    for var_num in 2 4 6 8; do
        for cluster_num in 2; do
            for loop_count in $(seq 0 $((LOOP - 1))); do
                for inner_count in $(seq 0 $((THREADS - 1))); do
                    seed=$((loop_count * THREADS + inner_count))
                    $BENCH --coverage $cov --variant-num "$var_num" --seed "$seed" --cluster-num "$cluster_num" >>"$RESULT" &
                done
                wait
            done
        done
    done
done