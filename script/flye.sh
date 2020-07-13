#!/bin/bash
#$ -S /bin/bash
#$ -q big.q
#$ -pe smp 12
#$ -cwd
#$ -V

FASTA=$1

if [ $# -gt 1 ];then
SIZE=$2
else
SIZE=3500000
fi

source /bio/package/miniconda/miniconda3/bin/activate /grid/yoshimura/conda/flye

/grid/yoshimura/conda/flye/bin/flye --threads 12 \
	--out-dir ${3} \
	--genome-size ${SIZE} \
	--keep-haplotypes \
	--pacbio-hifi $FASTA
