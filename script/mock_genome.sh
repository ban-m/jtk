#!/bin/bash
REFERENCE=~/work/sandbox/result/DBB_QBL_without_gaps.fasta
badread simulate --reference ${REFERENCE} --quantity 50x \
        --error_model pacbio --qscore_model pacbio \
        --identity 99,100,0.5 \
        --seed 10 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" > ./result/DBB_QBL_CCS_15Kbp_50x.fastq

badread simulate --reference ${REFERENCE} --quantity 10x \
        --error_model pacbio --qscore_model pacbio \
        --length 15000 1000 \
        --identity 99,100,0.5 \
        --seed 10 \
        --junk_reads 0 --random_reads 0 --chimeras 0\
        --start_adapter_seq "" \
        --end_adapter_seq "" > ./result/DBB_QBL_CCS_15Kbp_10x.fastq

