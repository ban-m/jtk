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
    - Align your reads to the reference and index the bam file. 
    - Use `bash ./script/extract_region.sh $REFERENCE $BAM $REGION > reads.fastq`
    - `$REFERENCE` should be a fasta file, `$BAM` should be an indexed bam file, and `$REGION` should be a region specification, such as `chr1:10000000-15000000`.
    - It is the same as `samtools view -OBAM $BAM $REGION | samtools fastq `
    - A region should be smaller than 10Mbp. We have not tested regions longer than 10Mbp.
1. Modify `example.toml` as you want.
    - See `example.toml` for the explanation of the parameters.
    - I recommend to use absolute path for the input and the output directory.
    - `JTK` would create the temporary file at the *current directory*. Please exec at the location where you have a write permission.
2. Run `jtk pipeline -p example.toml`
    - Several JSON files and assmbly graphs would be created.

If you stoped or a panic occured in `JTK`, you can resume the execution by 

1. Edit the TOML file by `sed -i -e "/resume/c resume = true"`
2. Run `jtk pipeline -p foo.toml`

## Reproduce/Test dataset

```
wget 
./target/release/jtk pipeline -p example.toml
```


## *Caution*

Please do not input reads comming from a region more than 10M bp long. It has not been tested.

## Info 

Contact: Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>
