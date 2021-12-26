#!/bin/bash
#$ -S /bin/bash
#$ -q centos7.q
#$ -pe smp 24
#$ -cwd
#$ -V
## Phasing module.
## Usage: bash phasing.sh ${FASTA} ${REFERENCE} ${OUT_DIR}
# >>> Conda
__conda_setup="$('/home/ban-m/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/ban-m/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/ban-m/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/ban-m/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< Conda
# set -ue
FASTA=$1
REFERENCE=$2
OUTDIR=$3
conda activate hla

## Variables.
## Alignment
BAM=${OUTDIR}/aln.bam

## Long-shot
UNPHASED_RAW=${OUTDIR}/unphased.raw.vcf
UNPHASED=${OUTDIR}/unphased.vcf

## WhatsHap phasing list.
PHASED=${OUTDIR}/phased.vcf.gz
PHASED_LIST=${OUTDIR}/whatshap_phase.tsv
PHASED_BAM=${OUTDIR}/phased_aln.bam

set -e
## 1. Long-shot to detect the contig with variants.
minimap2 -a -x map-pb --secondary=no ${REFERENCE} ${FASTA} |\
    samtools view -O BAM -@ 24 |\
    samtools sort -O BAM -@ 24 -m 10G > ${BAM}
samtools index ${BAM}
samtools faidx ${REFERENCE}
longshot -F --no_haps --out ${UNPHASED_RAW} --ref ${REFERENCE} --bam ${BAM}

## 3. Run whatshap.
cargo run --release --bin longshot_to_whatshap -- ${UNPHASED_RAW} > ${UNPHASED}

whatshap phase --indels -o ${PHASED} --reference ${REFERENCE} \
         --ignore-read-groups ${UNPHASED} ${BAM}
tabix ${PHASED}
whatshap haplotag -o ${PHASED_BAM} --ignore-read-groups \
         --reference ${REFERENCE} ${PHASED} ${BAM} \
         --output-haplotag-list ${PHASED_LIST}
whatshap stats --gtf=${OUTDIR}/phased_blocks.gtf ${PHASED}
samtools index ${PHASED_BAM}
