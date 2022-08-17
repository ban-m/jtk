# JKT -- regional diploid genome assembler

## Requirements 

- [minimap2](https://github.com/lh3/minimap2) with version >= 2.23
- [Rust](https://www.rust-lang.org/)

## Install 

After install minimap2 at the location included by `$PATH`, 

```
git clone 
cd 
cargo build --release 
```

Then, move `./target/release/jtk` at the location you want.


## Usage

1. Modify `example.toml` as you want.
2. Run `jtk pipeline -p example.toml`

Then, several JSON files and assmbly graphs would be created.

## *Caution*

Please do not input reads comming from a region more than 10M bp long. It has not been tested.

## Info 

Contact: Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>

Cite: 

Jittoku: [a person in the zen literature](https://en.wikipedia.org/wiki/Hanshan_and_Shide).
