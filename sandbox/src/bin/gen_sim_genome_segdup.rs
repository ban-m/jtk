// Synopsis:
// `cargo run --release --bin gen_sim_genome_segdup --${haps} ${seed}
// The output sequences have header like `>hapA/hapB`
// The hapA and hapB contains two large (1Mbp) segdup with 5% divergence rate, and 0.1% haplotype diff.
use kiley::gen_seq;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro128StarStar;
const PADDING: usize = 100_000;
const SEGDUP_LEN: usize = 1_000_000;
// Variant rate = 0.1%
const PROFILE: kiley::gen_seq::Profile = kiley::gen_seq::Profile {
    sub: 0.001 / 3f64,
    del: 0.001 / 3f64,
    ins: 0.001 / 3f64,
};

use std::io::{BufWriter, Write};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let mut hap_file = std::fs::File::create(&args[1]).map(BufWriter::new)?;
    let seed: u64 = args[2].parse().unwrap();
    let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(seed);
    let (hap_a, hap_b) = gen_haploids(&mut rng);
    writeln!(hap_file, ">hapA\n{}", String::from_utf8_lossy(&hap_a))?;
    writeln!(hap_file, ">hapB\n{}", String::from_utf8_lossy(&hap_b))?;
    Ok(())
}

fn gen_haploids<R: rand::Rng>(rng: &mut R) -> (Vec<u8>, Vec<u8>) {
    let segdup: Vec<_> = gen_seq::generate_seq(rng, SEGDUP_LEN);
    let profile = gen_seq::Profile {
        sub: 0.05 / 3f64,
        ins: 0.05 / 3f64,
        del: 0.05 / 3f64,
    };
    let segdup2 = gen_seq::introduce_randomness(&segdup, rng, &profile);
    let leading = gen_seq::generate_seq(rng, PADDING);
    let pad = gen_seq::generate_seq(rng, PADDING);
    let trail = gen_seq::generate_seq(rng, PADDING);
    let hap_a: Vec<_> = [leading, segdup, pad, segdup2, trail].concat();
    let hap_b = gen_seq::introduce_randomness(&hap_a, rng, &PROFILE);
    (hap_a, hap_b)
}
