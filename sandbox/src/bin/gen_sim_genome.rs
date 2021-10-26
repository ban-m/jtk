// Synopsis:
// `cargo run --release --bin gen_sim_genome -- ${hap1} ${hap2} ${seed}
// The output sequences have header like `>hapA/B`.
// HapA is exactly 1Mbp sequence, and hapB is constructed by introducing 0.1% mutation (indel,subs in total with equal probability) after concatinating the following sequence one after another:
// 1. hapA[..200Kbp]
// 2. hapA[220Kbp..420Kbp] (introducing L-bp deletion)
// 3. random 30Kbp sequence (introducing L-bp insertion)
// 4. hapA[420Kbp..620Kbp]
// 5. revcmp(hapA[620Kbp..630Kbp]) (introducing L-bp inversion)
// 6. hapA[630Kbp..830Kbp]
// 7. hapA[820Kbp..830Kbp]
// 8. hapA[820Kbp..830Kbp] (introducing L/3bp seg-dup x 3)
// 9. hapA[830Kbp..1Mbp].
// It would be very hard case for usual assembler or VariantCall+Phasing procedures.
use kiley::gen_seq;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro128StarStar;
const HAPA_LEN: usize = 1_000_000;
const MOD_LEN: usize = 50_000;
// Variant rate = 0.1%
const PROFILE: kiley::gen_seq::Profile = kiley::gen_seq::Profile {
    sub: 0.0005 / 3f64,
    del: 0.0005 / 3f64,
    ins: 0.0005 / 3f64,
};

use std::io::{BufWriter, Write};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let seed: u64 = args[3].parse().unwrap();
    let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(seed);
    let hap_a = gen_seq::generate_seq(&mut rng, HAPA_LEN);
    let hap_b = gen_hapb(&mut rng, &hap_a);
    let mut hapa_file = std::fs::File::create(&args[1]).map(BufWriter::new)?;
    let mut hapb_file = std::fs::File::create(&args[2]).map(BufWriter::new)?;
    writeln!(&mut hapa_file, ">hapA")?;
    for chunk in hap_a.chunks(120) {
        writeln!(&mut hapa_file, "{}", String::from_utf8_lossy(chunk))?;
    }
    writeln!(&mut hapb_file, ">hapB")?;
    for chunk in hap_b.chunks(120) {
        writeln!(&mut hapb_file, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

fn gen_hapb<R: rand::Rng>(rng: &mut R, hap_a: &[u8]) -> Vec<u8> {
    let inserted_seq = gen_seq::generate_seq(rng, MOD_LEN);
    let inverted_seq = bio_utils::revcmp(&hap_a[600_000 + MOD_LEN..620_000 + 2 * MOD_LEN]);
    let dup_region = &hap_a[800_000 + 2 * MOD_LEN..800_000 + 2 * MOD_LEN + MOD_LEN / 3];
    let hap_b: Vec<_> = hap_a[..200_000]
        .iter()
        .chain(&hap_a[200_000 + MOD_LEN..400_000 + MOD_LEN])
        .chain(inserted_seq.iter())
        .chain(&hap_a[400_000 + MOD_LEN..600_000 + MOD_LEN])
        .chain(inverted_seq.iter())
        .chain(&hap_a[600_000 + 2 * MOD_LEN..800_000 + 2 * MOD_LEN + MOD_LEN / 3])
        .chain(dup_region.clone())
        .chain(dup_region.clone())
        .chain(&hap_a[800_000 + 2 * MOD_LEN + MOD_LEN / 3..])
        .copied()
        .collect();
    gen_seq::introduce_randomness(&hap_b, rng, &PROFILE)
}
