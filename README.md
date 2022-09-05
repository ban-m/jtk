# JKT -- regional diploid genome assembler

JTK is a regional diploid genome assembler.

## Requirements 

- [minimap2](https://github.com/lh3/minimap2) with version >= 2.23
- [Rust](https://www.rust-lang.org/) with version > 1.59

## Installation

After installing minimap2 to the location included by `$PATH`, do

```
git clone https://github.com/ban-m/jtk.git
cd jtk
cargo build --release 
```

Then, move `./target/release/jtk` at the location you want.

## Usage

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
