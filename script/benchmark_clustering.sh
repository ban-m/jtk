#!/bin/bash
# Synopsis: `bash ./script/benchmark_clustering.sh $OUTPUT_DIR $THREADS`
BENCH="$PWD"/target/release/benchmark_clustering
THREADS=${2:-1}
mkdir -p "$1"
NUM_LOOP=1000
LOOP=$((NUM_LOOP / THREADS))
cargo build --release

RESULT="$1"/two_cluster.tsv
rm -f "$RESULT"
touch "$RESULT"
for cov in 5 10 20 30 40; do
    for len in 1000 2000 4000 8000; do
        for error_rate in 0.01 0.05 0.10 0.15; do
            for loop_count in $(seq 0 $((LOOP - 1))); do
                for inner_count in $(seq 0 $((THREADS - 1))); do
                    seed=$((loop_count * THREADS + inner_count))
                    $BENCH --coverage $cov --seed $seed --error-rate $error_rate --template-len $len >>"$RESULT" &
                done
                wait
            done
        done
    done
done

RESULT="$1"/four_cluster.tsv
rm -f "$RESULT"
touch "$RESULT"
for cov in 10 20 30; do
    for error_rate in 0.01 0.05 0.10 0.15; do
        for len in 1000 2000 4000 8000; do
            for loop_count in $(seq 0 $((LOOP - 1))); do
                for inner_count in $(seq 0 $((THREADS - 1))); do
                    seed=$((loop_count * THREADS + inner_count))
                    $BENCH --coverage $cov --seed $seed --error-rate $error_rate --cluster-num 4 --template-len $len >>"$RESULT" &
                done
                wait
            done
        done
    done
done
