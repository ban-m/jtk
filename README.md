# JKT -- HLA decomposer

Author: Bansho Masutani

Mail: ban-m@g.ecc.u-tokyo.ac.jp

# Desc

It has ONLY ONE command line interface, `jtk.` `jtk` (short for Jittoku or Japanese HAL analysis toolkit) is THE entry point of entire analysis. It provides ways to modify dataset, inspect intermidiate file, and generate filnal product.

To see the detail, run `cargo build --release` and then `./target/release/jtk --help`.

## Implementation details


### Unit selection

1. Sort reads in descending order w.r.t the length of the read.
2. Until we collect sufficient number of chunks or there is no read, do:
   1. Pop the longest read.
   2. Take the L-length chunk every K-base skipping interval



### Unit encoding 

# Memos


## TODO

### Documentation

Currently, I've written documents on only a few types. We should write more and more ducumentation, before we forget the meaning of fragments of codes.


### Sanity check

At each stage of the pipeline, we should check whether the input data has enough information; In the clustering step, the dataset should have, at least, an encoded read set and selected units.

Also, we should have some "invariant checker" for sanity check. For example, we need to check the dataset is consistent. For example, all encoded reads should have corresponding raw reads. This checker is required, as in the future we filter out raw reads and encoded reads at some stage among the pipiline, or there would be some interplay between Python <-> Rust or Javascript <-> Rust.

### Compact serialization(BinCode) vs Readable serialization(JSON)

- It is preferrable to switch between them easily. Currently I think the message passing should be based on JSON by default, as sometimes we do need some help of other language or script such as Python, Ruby, or even Javascript.
However, as the final product or 'publishable state', it would be more elegant to impelemt them entirely in Rust language, where the use of JSON does not make any sense.


# Reference

The name, Jittoku, is taken from [a person in Zen literature](https://en.wikipedia.org/wiki/Hanshan_and_Shide).