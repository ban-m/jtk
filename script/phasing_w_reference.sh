#!/bin/bash
## Phasing module.
## Usage: bash phasing.sh ${FASTQ} ${REFERENCE} ${OUT_DIR}
## Requirements: Whatshap, minimap2, longshot, samtools, cargo.
set -ue
FASTQ=$1
REFERENCE=$2
OUTDIR=$3

## Variables.
### Alignment
BAM=${OUTDIR}/aln.bam

### Long-shot
UNPHASED_RAW=${OUTDIR}/unphased.raw.vcf
UNPHASED=${OUTDIR}/unphased.vcf

### WhatsHap phasing list.
PHASED=${OUTDIR}/phased.vcf.gz
PHASED_LIST=${OUTDIR}/whatshap_phase.tsv
PHASED_BAM=${OUTDIR}/phased_aln.bam

## 1. Long-shot to detect the contig with variants.
minimap2 -a -x map-pb --secondary=no ${REFERENCE} ${FASTQ} |\
    samtools view -O BAM -@ 24 |\
    samtools sort -O BAM -@ 24 -m 1G > ${BAM}
samtools index ${BAM}
samtools faidx ${REFERENCE}
longshot --min_alt_count 5 -F --no_haps --out ${UNPHASED_RAW} --ref ${REFERENCE} --bam ${BAM}

## 3. Run whatshap.
cargo run --release --bin longshot_to_whatshap -- ${UNPHASED_RAW} > ${UNPHASED}

whatshap phase --distrust-genotypes -o ${PHASED} --reference ${REFERENCE} \
         --ignore-read-groups ${UNPHASED} ${BAM}
tabix ${PHASED}
whatshap haplotag -o ${PHASED_BAM} --ignore-read-groups \
         --reference ${REFERENCE} ${PHASED} ${BAM} \
         --output-haplotag-list ${PHASED_LIST}
whatshap stats --gtf=${OUTDIR}/phased_blocks.gtf ${PHASED}
samtools index ${PHASED_BAM}
