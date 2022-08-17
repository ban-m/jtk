#!/bin/bash
#$ -S /bin/bash
#$ -N Map
#$ -cwd
#$ -pe smp 24
#$ -m e 
#$ -V
#$ -j y
#$ -o ./data/maphg02.2.log
set -ue
## Alignment
THREADS=24
INDEX=/grid/ban-m/human_genome/chm13.draft_v1.1.mmi
REGION="chr6:28573494-33272792"
mkdir -p ${PWD}/data/hg002_filtered/
for accession in `seq 9972588 9972596`
do
    for READ in /grid/ban-m/human_reads/clr_hg002/SRR${accession}/*.fastq.gz
    do
        BAM=${READ%.fastq.gz}.sorted.bam
        minimap2 -a -t ${THREADS} ${INDEX} ${READ} |\
            samtools view -O BAM -@ ${THREADS} |\
            samtools sort -O BAM -@ ${THREADS} -m 10G > ${BAM}
        samtools index ${BAM}
        # remove secondary alignment.
        samtools view -O BAM -@ ${THREADS} -F 256 ${BAM} ${REGION}|\
            samtools sort -O BAM -@ ${THREADS} -m 10G > ${PWD}/data/hg002_filtered/${accession}.bam
        samtools index ${PWD}/data/hg002_filtered/${accession}.bam
        samtools fasta ${PWD}/data/hg002_filtered/${accession}.bam > ${PWD}/data/hg002_filtered/${accession}.fa
    done
done
