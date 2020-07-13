#!/bin/bash
# bash canu.sh ${FASTA} ${GENOME_SIZE} ${DIRECTRY}
FASTA=$1
JOBN=`basename $FASTA`

if [ $# -gt 1 ];then
SIZE=$2
else
SIZE=3500000
fi

unset PERL5LIB
unset LD_LIBRARY_PATH
source /bio/package/miniconda/miniconda3/bin/activate /grid/yoshimura/conda/canu
export PATH=/usr/bin:$PATH
source /grid/sgeadmin/default/common/settings.sh
export PERL5LIB=/grid/yoshimura/conda/canu/lib
/grid/yoshimura/conda/canu/bin/canu --version 
/grid/yoshimura/conda/canu/bin/canu -d ${3} -p ${JOBN%.fa*}.canu \
	genomeSize=$SIZE \
	gridOptions="-S /bin/bash -V -q all.q" \
	gridEngineResourceOption="-pe smp THREADS -l mem_free=MEMORY" \
	minInputCoverage=1 \
	stopOnLowCoverage=1 \
	-pacbio-hifi $FASTA
