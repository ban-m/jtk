#!/bin/bash
## Phasing module.
## Usage: bash phasing.sh ${FAST[Q|A]} ${REFERENCE} ${OUTDIR} $PEPPER_SIF_FILE $THREADS
## Requirements: Longphase, minimap2, cuteSV, samtools, seqtk, flye, and Pepper-margin-deepvariant.
set -ue
READS=$1
REFERENCE=$2
OUTDIR=$3
PMD=$4
THREADS=$5
## Variables.
BAM=${OUTDIR}/aln.bam
PHASED_BAM=${OUTDIR}/aln.phased.bam
SNV_PREFIX=pmd
SV_PREFIX=cutesv


## 0. alignments
mkdir -p ${OUTDIR}
minimap2 --MD -a -x map-ont -t ${THREADS} ${REFERENCE} ${READS} | samtools view -OBAM | samtools sort -@ $THREADS -m1G -OBAM > $BAM
samtools index ${BAM}

## 1. call variants.
REF_STEM=${REFERENCE%/*.fa}
singularity exec --bind /usr/lib/locale/,${REF_STEM},${OUTDIR},${PWD} $PMD\
    run_pepper_margin_deepvariant call_variant -b ${BAM} -f ${REFERENCE} -o ${OUTDIR} -p $SNV_PREFIX -t${THREADS} --ont_r9_guppy5_sup 2> log   
cuteSV ${BAM} ${REFERENCE} ${OUTDIR}/${SV_PREFIX}.vcf ${OUTDIR} --report_readid --genotype \
    --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3

## 2. Phase 
longphase_linux-x64 phase -s ${OUTDIR}/${SNV_PREFIX}.vcf.gz --sv-file ${OUTDIR}/${SV_PREFIX}.vcf -b $BAM -r ${REFERENCE} -t ${THREADS} -o ${OUTDIR}/phased  --ont 
longphase_linux-x64 haplotag -s ${OUTDIR}/phased.vcf --sv-file ${OUTDIR}/phased_SV.vcf --log -b ${BAM} -t $THREADS -o ${OUTDIR}/haplotag

## 3. Assemble each haplotype
cat ${OUTDIR}/haplotag.out | awk '(3 < NF && $1 ~ !/^#/ && $5 != "1")' | cut -f1 > ${OUTDIR}/h2.ids
cat ${OUTDIR}/haplotag.out | awk '(3 < NF && $1 ~ !/^#/ && $5 != "2")' | cut -f1 > ${OUTDIR}/h1.ids
seqtk subseq ${READS} ${OUTDIR}/h1.ids > ${OUTDIR}/h1.fq
seqtk subseq ${READS} ${OUTDIR}/h2.ids > ${OUTDIR}/h2.fq
flye --nano-raw ${OUTDIR}/h1.fq --out-dir ${OUTDIR}/h1 --genome-size 5M -t$THREADS 
flye --nano-raw ${OUTDIR}/h2.fq --out-dir ${OUTDIR}/h2 --genome-size 5M -t$THREADS 