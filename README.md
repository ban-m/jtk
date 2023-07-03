# JKT -- Targeted Diploid Genome Assembler

JTK is a targeted diploid genome assembler.

## Requirements 

- [minimap2](https://github.com/lh3/minimap2) with version >= 2.23.
- [Rust](https://www.rust-lang.org/) with version >= 1.72.0 **nightly**.

## Installation

First, check the version of the Rust language.

```
cargo --version
```

The minimum required version is 1.72.0-nightly (otherwise, `rustup update` so that the version is greater than or equal to 1.72.0-nightly). Then,


```
git clone https://github.com/ban-m/jtk.git
cd jtk
cargo build --release
./target/release/jtk --version 
minimap2 --verion # minimap2 is needed.
```

Move `./target/release/jtk` to the location in $PATH.

## Usage

Suppose we have the following variables -- 

- `REFERENCE`: a reference sequence in FASTA format.
- `READS`: reads in FASTA format (recommended: ONT reads)
- `REGION`: target region in [CHR]:[start]-[end] format (e.g., chr1:10000000-15000000). This region should be smaller than 10Mbp and should not end/starts in a segmental duplication.



0. Prepare the reads for assembly.
   - Map your reads to the reference and index the bam file with `minimap2`:  
   ```
   minimap2 -x map-ont -t $THREADS --secondary=no -a $REFERENCE $READS |\
      samtools sort -@$THREADS -OBAM > aln.bam && samtools index aln.bam
   ```
   - `bash ./script/extract_region.sh $REFERENCE $BAM $REGION > reads.fastq`  
   This is equivalent to `samtools view -OBAM $BAM $REGION | samtools fastq > reads.fastq`.
1. Modify `example.toml`
   - See `example.toml` for an explanation of the parameters.
   - I recommend using absolute paths for the input and the output directory.
   - `JTK` would create the temporary file at the *current directory*. Please exec at the location where you have a write permission.
   - `sed` is useful. For example, `cat example.toml | sed -e "/^input_file/c input_file = \""$DATA"\"" ... > profile.toml` would replace the input file with `$DATA`.
2. Run `jtk pipeline -p example.toml`
   - This would create several JSON files and assmbly graphs.
   - In addition, `prefix.sam` is the alignment between the reads and the assembly, and `prefix.coverage.tsv` is the coverage trace on the assembly.

## Test Data Set

```
minimap2 --version # Should be greater than 2.23 
jtk --version # Should be greater than 0.1
wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/jtk/COX_PGF.fastq.gz
gunzip COX_PGF.fastq.gz
wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/jtk/COX_PGF.toml
jtk pipeline -p COX_PGF.toml 2> log
```

Then, the assembly graph is available at `cox_pgf/temp.gfa`.


## *Caution*

Please do not enter reads coming from a region longer than 10M bp long. It has not been tested.

## Info 

Contact: Bansho Masutani banmasutani@gmail.com