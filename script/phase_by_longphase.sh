#!/bin/bash
## Phasing module.
## Usage: bash phasing.sh ${FAST[Q|A]} ${REFERENCE} "$OUTDIR" $PEPPER_SIF_FILE $THREADS
## Requirements: Longphase, minimap2, cuteSV, samtools, seqtk, flye, and Pepper-margin-deepvariant.
set -ue
READS=$1
REFERENCE=$2
OUTDIR=$3
PMD=$4
THREADS=$5
## Variables.
BAM="$OUTDIR"/aln.bam
SNV_PREFIX=pmd
SV_PREFIX=cutesv

## 0. alignments
mkdir -p "$OUTDIR"
minimap2 --MD -a -x map-ont -t "$THREADS" "$REFERENCE" "$READS" | samtools view -OBAM | samtools sort -@ "$THREADS" -m1G -OBAM >"$BAM"
samtools index "$BAM"

## 1. call variants.
REF_STEM=${REFERENCE%/*.fa}
singularity exec \
    --bind /usr/lib/locale/ --bind "$REF_STEM" \
    --bind "$OUTDIR" \
    --bind "$PWD" \
    "$PMD" run_pepper_margin_deepvariant call_variant \
    -b "$BAM" -f "$REFERENCE" -o "$OUTDIR" -p "$SNV_PREFIX" -t "$THREADS" --ont_r9_guppy5_sup
cuteSV "$BAM" "$REFERENCE" "$OUTDIR"/"$SV_PREFIX".vcf "$OUTDIR" --report_readid --genotype \
    --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3

## 2. Phase
longphase_linux-x64 phase -s "$OUTDIR"/"$SNV_PREFIX".vcf.gz --sv-file "$OUTDIR"/${SV_PREFIX}.vcf -b "$BAM" -r "$REFERENCE" -t "$THREADS" -o "$OUTDIR"/phased --ont
longphase_linux-x64 haplotag -s "$OUTDIR"/phased.vcf --sv-file "$OUTDIR"/phased_SV.vcf --log -b "$BAM" -t "$THREADS" -o "$OUTDIR"/haplotag

## 3. Assemble each haplotype
awk '(3 < NF && $1 ~ !/^#/ && $5 != "1")' "$OUTDIR"/haplotag.out | cut -f1 >"$OUTDIR"/h2.ids
awk '(3 < NF && $1 ~ !/^#/ && $5 != "2")' "$OUTDIR"/haplotag.out | cut -f1 >"$OUTDIR"/h1.ids
seqtk subseq "$READS" "$OUTDIR"/h1.ids >"$OUTDIR"/h1.fq
seqtk subseq "$READS" "$OUTDIR"/h2.ids >"$OUTDIR"/h2.fq
if [ -e "/usr/bin/time" ]; then
    /usr/bin/time -v flye --nano-raw "$OUTDIR"/h1.fq --out-dir "$OUTDIR"/h1 --genome-size 5M -t"$THREADS"
    /usr/bin/time -v flye --nano-raw "$OUTDIR"/h2.fq --out-dir "$OUTDIR"/h2 --genome-size 5M -t"$THREADS"
else
    flye --nano-raw "$OUTDIR"/h1.fq --out-dir "$OUTDIR"/h1 --genome-size 5M -t"$THREADS"
    flye --nano-raw "$OUTDIR"/h2.fq --out-dir "$OUTDIR"/h2 --genome-size 5M -t"$THREADS"
fi
awk 'BEGIN{num=0}($0 ~ />/){if(num != 0){printf("\n")} printf(">contig_h1_%d\n", num); num +=1;next}{printf("%s", $0)}END{printf("\n")}' "$OUTDIR"/h1/assembly.fasta > "$OUTDIR"/final_assembly.fasta
awk 'BEGIN{num=0}($0 ~ />/){if(num != 0){printf("\n")} printf(">contig_h2_%d\n", num); num +=1;next}{printf("%s", $0)}END{printf("\n")}' "$OUTDIR"/h2/assembly.fasta >> "$OUTDIR"/final_assembly.fasta