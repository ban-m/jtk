#!/bin/bash
# See https://www.ncbi.nlm.nih.gov/nuccore/GL000251.2 to confirm this is the haplotype in the COX cell line
echo ">chr6_GL000251v2_alt" > ${PWD}/data/COX.fa
samtools faidx /grid/ban-m/human_genome/hg38.fa chr6_GL000251v2_alt | tail -n+2 |\
    tr -d '\nN'| awk '{print(toupper($0))}' | fold -w 80  >> ${PWD}/data/COX.fa
    

# Do `minimap2 -c <( samtools faidx /grid/ban-m/human_genome/hg38.fa chr6 ) cox.fa > aln.paf` to convince yourself
# that the MHC region in the hg38 is chr6:28734407-33411973
echo ">chr6:28734407-33411973" > ${PWD}/data/PGF.fa
samtools faidx /grid/ban-m/human_genome/hg38.fa chr6:28734407-33411973 | tail -n+2 |\
    tr -d '\nN'| awk '{print(toupper($0))}' | fold -w 80 >> ${PWD}/data/PGF.fa

cat ${PWD}/data/COX.fa ${PWD}/data/PGF.fa > ${PWD}/data/COX_PGF.fa
