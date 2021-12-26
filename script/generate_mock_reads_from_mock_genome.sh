#!/bin/bash
set -ue 
OUTPUT=${PWD}/data/mock_genome/
mkdir -p ${OUTPUT}
CLR_READS=${OUTPUT}/haps.clr
CCS_READS=${OUTPUT}/haps.ccs
HAPA=${OUTPUT}/hapA.fa
HAPB=${OUTPUT}/hapB.fa
REFERENCE=${OUTPUT}/hap_ref.fa
MERGED=${OUTPUT}/haps.fa
SEED=3490283
SEED_BR=39120
RENAME=${PWD}/script/rename_fastq.awk
# Gen genome
cargo run --release --bin gen_sim_genome -- ${REFERENCE} ${HAPA} ${HAPB} ${SEED} 
cat ${HAPA} ${HAPB} > ${MERGED}

for quantity in 30x 50x
do
    CLR_READS_QU=${CLR_READS}.${quantity}
    badread simulate --reference ${MERGED} --quantity ${quantity} \
            --error_model pacbio2016 --qscore_model pacbio2016 \
            --length 25000,2000 \
            --identity 85,95,5 \
            --seed ${SEED_BR} \
            --junk_reads 0 --random_reads 0 --chimeras 0\
            --start_adapter_seq "" \
            --end_adapter_seq "" |\
        awk -f ${RENAME} > ${CLR_READS_QU}.fq
    cat ${CLR_READS_QU}.fq | paste - - - - | cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${CLR_READS_QU}.fa
done

for quantity in 15x 30x
do
    CCS_READS_QU=${CCS_READS}.${quantity}
    badread simulate --reference ${MERGED} --quantity ${quantity} \
            --error_model pacbio2016 --qscore_model pacbio2016 \
            --length 25000,2000 \
            --identity 99.9,100,0.01 \
            --seed ${SEED_BR} \
            --junk_reads 0 --random_reads 0 --chimeras 0\
            --start_adapter_seq "" \
            --end_adapter_seq "" |\
        awk -f ${RENAME} > ${CCS_READS_QU}.fq
    cat ${CCS_READS_QU}.fq | paste - - - - | cut -f1,2 | sed -e 's/@/>/g' | tr '\t' '\n' > ${CCS_READS_QU}.fa
done
