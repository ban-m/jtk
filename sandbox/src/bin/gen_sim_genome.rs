// Synopsis:
// `gen_sim_genome --ref-out-path ${reference}\
//    --hap1-out-path ${hap1}\
//    --hap2-out-path ${hap2}\
//    --seed ${seed}
// The output sequences have header like `>reference/hapA/hapB`
// `reference` is exactly 1Mbp sequence, and hapA and hapB are constructed by introducing 0.05% mutation (indel,subs in total with equal probability) after concatinating the following sequence one after another:
// For hapA
// 1. ref[0..200Kbp] with two deletion at 50K..80K and 100K..120Kbp
// 2. ref[200Kbp..400Kbp] with inversion of 120K..170Kbp.
// 3. ref[400K..600Kbp] with 130K..180Kbp deletion
// 4. ref[600Kbp..800Kbp] with 50Kbp insertion at 700Kbp.
// 5. ref[800Kbp..].
// For hapB.
// 1. ref[0..200Kbp] with two deletion at 50K..80K and 100K..120Kbp
// 2. ref[200Kbp..400Kbp] with inversion of 120K..170Kbp.
// 3. ref[400K..600Kbp]
// 4. ref[600Kbp..800Kbp] with 50Kbp insertion at 700Kbp.
// 5. ref[800Kbp..] with 50Kbp insertion at 900Kbp..
// It would be very hard case for usual assembler or VariantCall+Phasing procedures.
use clap::Parser;
use kiley::gen_seq;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro128StarStar;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};
const HAPA_LEN: usize = 1_000_000;
const MOD_LEN: usize = 50_000;
// Variant rate = 0.1%
const PROFILE: kiley::gen_seq::Profile = kiley::gen_seq::Profile {
    sub: 0.0005 / 3f64,
    del: 0.0005 / 3f64,
    ins: 0.0005 / 3f64,
};

#[derive(Parser)]
#[command(name = "gen_sim_genome")]
#[command(version = "0.1")]
#[command(about = "Generate a random genome")]
struct Cli {
    #[arg(long)]
    ref_out_path: PathBuf,
    #[arg(long)]
    hap1_out_path: PathBuf,
    #[arg(long)]
    hap2_out_path: PathBuf,
    /// Seed for the pseudo random number generator
    #[arg(long)]
    seed: u64,
}

fn main() -> std::io::Result<()> {
    let args = Cli::parse();
    let seed = args.seed;
    let mut rng: Xoshiro128StarStar = SeedableRng::seed_from_u64(seed);
    let reference = gen_seq::generate_seq(&mut rng, HAPA_LEN);
    let (hap_a, hap_b) = gen_haploids(&mut rng, &reference);
    let mut refr_file = std::fs::File::create(&args.ref_out_path).map(BufWriter::new)?;
    let mut hapa_file = std::fs::File::create(&args.hap1_out_path).map(BufWriter::new)?;
    let mut hapb_file = std::fs::File::create(&args.hap2_out_path).map(BufWriter::new)?;
    let reference = String::from_utf8_lossy(&reference);
    writeln!(refr_file, ">reference\n{}", reference)?;
    writeln!(hapa_file, ">hapA\n{}", String::from_utf8_lossy(&hap_a))?;
    writeln!(hapb_file, ">hapB\n{}", String::from_utf8_lossy(&hap_b))?;
    Ok(())
}

fn gen_haploids<R: rand::Rng>(rng: &mut R, refr: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let seq1 = refr
        .iter()
        .take(200_000)
        .enumerate()
        .filter_map(|(pos, b)| {
            (!(50_000..80_000).contains(&pos) && !(100_000..120_000).contains(&pos)).then_some(b)
        });
    let inv_seq = bio_utils::revcmp(&refr[320_000..370_000]);
    let seq2 = refr[200_000..320_000]
        .iter()
        .chain(inv_seq.iter())
        .chain(refr[370_000..400_000].iter());
    let seq3_a = refr[400_000..600_000]
        .iter()
        .enumerate()
        .filter_map(|(pos, b)| (!(130_000..180_000).contains(&pos)).then_some(b));
    let seq3_b = refr[400_000..600_000].iter();
    let inserted_seq = gen_seq::generate_seq(rng, MOD_LEN);
    let seq4 = refr[600_000..700_000]
        .iter()
        .chain(inserted_seq.iter())
        .chain(refr[700_000..800_000].iter());
    let seq5_a = refr[800_000..].iter();
    let inserted_seq_b = gen_seq::generate_seq(rng, MOD_LEN);
    let seq5_b = refr[800_000..900_000]
        .iter()
        .chain(inserted_seq_b.iter())
        .chain(refr[900_000..].iter());
    let hap_a: Vec<_> = seq1
        .clone()
        .chain(seq2.clone())
        .chain(seq3_a)
        .chain(seq4.clone())
        .chain(seq5_a)
        .copied()
        .collect();
    let hap_b: Vec<_> = seq1
        .chain(seq2)
        .chain(seq3_b)
        .chain(seq4)
        .chain(seq5_b)
        .copied()
        .collect();
    let hap_a = gen_seq::introduce_randomness(&hap_a, rng, &PROFILE);
    let hap_b = gen_seq::introduce_randomness(&hap_b, rng, &PROFILE);
    (hap_a, hap_b)
}
