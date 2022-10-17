# JKT -- regional diploid genome assembler

JTK is a regional diploid genome assembler.

## Requirements 

- [minimap2](https://github.com/lh3/minimap2) with version >= 2.23
- [Rust](https://www.rust-lang.org/) with version >= 1.66.0 **nightly**

## Installation

First, check the version of the Rust language.

```
cargo --version
```

If it is not up to date, please `rustup update` so that the version would be more than or equal to 1.66.0-nightly.

After installing minimap2 to the location included by `$PATH`, do

```
git clone https://github.com/ban-m/jtk.git
cd jtk
cargo build --release 
```

Then, move `./target/release/jtk` at the location you want.

## Usage

0. Prepare the reads to be assembled.
    - Align your reads to the reference and index the bam file (`--secondary=no` if available). 
    - Use `bash ./script/extract_region.sh $REFERENCE $BAM $REGION > reads.fastq`
    - `$REFERENCE` should be a fasta file, `$BAM` should be an indexed bam file, and `$REGION` should be a region specification, such as `chr1:10000000-15000000`.
    - It is the same as `samtools view -OBAM $BAM $REGION | samtools fastq `
    - A region should be smaller than 10Mbp and should not end with a segmental duplication (TODO: this should be fixed).
1. Modify `example.toml` as you want.
    - See `example.toml` for the explanation of the parameters.
    - I recommend to use absolute path for the input and the output directory.
    - `JTK` would create the temporary file at the *current directory*. Please exec at the location where you have a write permission.
    - `sed` is useful. For example, `cat example.toml | sed -e "/^input_file/c input_file = \""$DATA"\"" ... > profile.toml` would replace the input file with `$DATA`.
2. Run `jtk pipeline -p example.toml`
    - Several JSON files and assmbly graphs would be created.


Also, there is an agnostic subcommend, `jtk stats -f out.stat $JSON_FILE > /dev/null` and `out.stat` contains some useful information. 
(Caution: This command dumps the input file to the output file as-is. So do not forget to redirect to the `/dev/null`. )

## Test dataset

```
minimap2 --version # Should be larger than 2.23 
wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/jtk/COX_PGF.fastq.tar.gz
gunzip COX_PGF.fastq.tar.gz
wget https://mlab.cb.k.u-tokyo.ac.jp/~ban-m/jtk/COX_PGF.toml
jtk pipeline -p COX_PGF.toml 2> log # This commend uses 20 threads.
```

Then, the assembly graph is available at `cox_pgf/temp.gfa`.


## *Caution*

Please do not input reads comming from a region more than 10M bp long. It has not been tested.

## Info 

Contact: Bansho Masutani ban-m@g.ecc.u-tokyo.ac.jp 
