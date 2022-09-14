#!/bin/bash
set -uex
## Synopsis:
REFERENCE=$1
ANNOTATION=$2
DIPLO_ASM=$3
DIPLO_ASM_H1=$4
DIPLO_ASM_H2=$5
DIPLO_NAME=$6
JTK_ASM=$7
JTK_ASM_H1=$8
JTK_ASM_H2=$9
JTK_NAME=${10}
CHR=${11}
START=${12}
END=${13}
PREFIX=${14}

mkdir -p "$PREFIX"

### 0. extract diploid assemblies into four parts.
function cut_and_rename (){
    ASM=$1
    CTG_NAMES=$2
    HAP_NAME=$3
    HAP_NUM=$4
    TEMP=$RANDOM
    for ctgname in $CTG_NAMES
    do
        samtools faidx "$ASM" "$ctgname" >> "$TEMP".fa 
    done
    awk --assign name="$HAP_NAME" --assign hnum="$HAP_NUM"\
        '($0 ~ />/){sum+=1; print(">" name "_" hnum "_" sum);next}{print $0}'\
        "$TEMP".fa > "$PREFIX"/"$HAP_NAME"_h"$HAP_NUM".fa
    rm "$TEMP".fa
}

cut_and_rename "$DIPLO_ASM" "$DIPLO_ASM_H1" "$DIPLO_NAME" 1
cut_and_rename "$DIPLO_ASM" "$DIPLO_ASM_H2" "$DIPLO_NAME" 2
cut_and_rename "$JTK_ASM" "$JTK_ASM_H1" "$JTK_NAME" 1
cut_and_rename "$JTK_ASM" "$JTK_ASM_H2" "$JTK_NAME" 2

# samtools faidx "$JTK_ASM" "$JTK_ASM_H1" |\
#     awk --assign name="$JTK_NAME" \
#     '($0 ~ />/){sum+=1;print(">" name "_1_" sum);next}{print $0}' > "$PREFIX"/"$JTK_NAME"_h1.fa
# samtools faidx "$JTK_ASM" "$JTK_ASM_H2" |\
#     awk --assign name="$JTK_NAME" \
#     '($0 ~ />/){sum+=1;print(">" name "_2_" sum);next}{print $0}' > "$PREFIX"/"$JTK_NAME"_h2.fa

## TODO: Check there is only one reference sequence.
awk '($0 ~ />/){print(">ref");next}{print $0}' "$REFERENCE" > "$PREFIX"/reference.fa

### 1. liftover alignments
function liftover_by_liftoff () {
    GFFFILE=$1
    TARGET_CTG=$2
    REF_CTG=$3
    OUT_GFF=${TARGET_CTG%.fa}.gff3
    OUT_TSV=${TARGET_CTG%.fa}.tsv
    TEMP="$RANDOM"
    liftoff -o "$OUT_GFF" -dir "$TEMP" -g "$GFFFILE" "$TARGET_CTG" "$REF_CTG"
    grep -v "^#" "$OUT_GFF" |\
        awk 'BEGIN{OFS="\t"}($3 ~/gene/){if ($7 == "+"){ dir = "1"} else { dir = "-1"};print($1,$4,$5,dir,"gene")}' > "$OUT_TSV"
    rm -r "$TEMP"
}

cargo run --release --bin filter_gff3 -- "$ANNOTATION" "$CHR" "$START" "$END" ref > "$PREFIX"/filtered.gff3
grep -v "^#" "$PREFIX"/filtered.gff3 |\
    awk 'BEGIN{OFS="\t"}($3 ~ /gene/){if ($7 == "+"){ dir = "1"} else { dir = "-1"};print($1,$4,$5,dir,"gene")}' > "$PREFIX"/filtered.tsv
for CTG in "$DIPLO_NAME"_h1.fa "$DIPLO_NAME"_h2.fa "$JTK_NAME"_h1.fa "$JTK_NAME"_h2.fa
do
    TARGET="$PREFIX"/"$CTG"
    liftover_by_liftoff "$PREFIX"/filtered.gff3 "$TARGET" "$PREFIX"/reference.fa
done

### 2. All vs All alignments
cat "$PREFIX"/"$DIPLO_NAME"_h1.fa "$PREFIX"/"$DIPLO_NAME"_h2.fa\
    "$PREFIX"/"$JTK_NAME"_h1.fa "$PREFIX"/"$JTK_NAME"_h2.fa \
    "$PREFIX"/reference.fa |
    sed -e 's/>\(.*\)/>db_\1/g'> "$PREFIX"/seqs.fa
lastdb "$PREFIX"/seqs.db "$PREFIX"/seqs.fa
for seq in "$DIPLO_NAME"_h1.fa "$DIPLO_NAME"_h2.fa "$JTK_NAME"_h1.fa "$JTK_NAME"_h2.fa reference.fa
do
    outfile="$PREFIX"/${seq%.fa}.maf
    echo "##maf version=1 scoring=lastz.v1.03.73" > "$outfile"
    lastal -E0.0001 "$PREFIX"/seqs.db "$PREFIX"/"$seq" | last-split -r >> "$outfile"
done
