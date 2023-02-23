# JKT -- Regional Diploid Genome Assembler

JTK is a regional diploid genome assembler.

## Requirements 

- [minimap2](https://github.com/lh3/minimap2) with version >= 2.23.
- [Rust](https://www.rust-lang.org/) with version >= 1.66.0 **nightly**.

## Installation

First, check the version of the Rust language.

```
cargo --version
```

If it is not up to date, please `rustup update` so that the version is greater than or equal to 1.66.0-nightly. Then run,


```
git clone https://github.com/ban-m/jtk.git
cd jtk
cargo build --release
./target/release/jtk --version 
minimap2 --verion
```

Move `./target/release/jtk` to the location in $PATH.

## Usage


0. Prepare the reads for assembly.
   - Map your reads to the reference and index the bam file with `minimap2`. For example, `minimap2 -x map-ont -t $THREADS --secondary=no -a $REFERENCE $READS | samtools sort -@$THREADS -OBAM > aln.bam && samtools index aln.bam`.
   - Use `bash ./script/extract_region.sh $REFERENCE $BAM $REGION > reads.fastq`, where `$REGION` should be a region specification, such as `chr1:10000000-15000000`.
   - This is equivalent to `samtools view -OBAM $BAM $REGION | samtools fastq > reads.fastq`
   - A region should be smaller than 10Mbp and should not end in a segmental duplication (TODO: this should be fixed).
1. Modify `example.toml` as you wish.
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

Contact: Bansho Masutani ban-m@g.ecc.u-tokyo.ac.jp