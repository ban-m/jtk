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

`rg TODO` would show where to check.

### Unit encoding
LAST -> Something else. Hopefully, it is written in Rust. The priority of this TODO is extremely low, as LAST would serve a good aligner as far.

###  Documentation
Currently, I've written documents on only a few types. We should write more and more documentation before we forget the meaning of fragments of codes.

### Determine the number of the cluster

We should predict the number of the cluster for local clustering for each determiend units.
Although we can predict the number of the cluster as well as clusterings simultaneously (Variational Byes'?), it might be though to tune hyper parameters or find an approprieate model.


### Clustering (Global)
The entire algorithm should be more mature and sophisticated. We need to find some foundation of our clustering algorithm. Maybe a theory on Markov's walk on graphs serves a good guide.


### Sanity check
At each stage of the pipeline, we should check whether the input data has enough information. In the clustering step, the dataset should have, at least, an encoded read set and selected units.
Also, we should have some "invariant checker" for a sanity check. For example, we need to check the dataset is consistent. For instance, all encoded reads should have corresponding raw reads. This checker is required, as we filter out raw reads and encoded reads at some stage among the pipeline in the future. Also, there would be some interplay between Python <-> Rust or Javascript <-> Rust.



### Compact serialization(BinCode) vs Readable serialization(JSON)

It is preferable to switch between them easily. Currently, message passing is based on JSON, as we need some help from other scripts such as Python and Javascript.
However, when we run the pipeline only with `jtk`, which is entirely written in Rust, the use of JSON does not make any sense.
The priority of this TODO is extremely low. The pipeline is not broken, and we do not have to fix something not broken.

# Reference

The name, Jittoku, is taken from [a person in Zen literature](https://en.wikipedia.org/wiki/Hanshan_and_Shide).