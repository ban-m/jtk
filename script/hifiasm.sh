#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -pe smp 12
#$ -cwd
#$ -V

FASTA=$1

if [ $# -gt 1 ];then
SIZE=$2
else
SIZE=3500000
fi

/grid/yoshimura/conda/hifiasm/hifiasm \
	-t 12 \
	-o ${3} \
	$FASTA
