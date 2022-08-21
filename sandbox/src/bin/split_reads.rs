/// Synopsis: [fasta file] [into N chunks] [with prefix]
use bio_utils::fastq::*;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    // First argument would be the name of the binary...
    let args: Vec<_> = std::env::args().collect();
    let mut records = parse_into_vec(&args[1]).unwrap();
    let prefix = &args[2];
    let amount = parse_si(&args[3]).unwrap();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(34204);
    records.shuffle(&mut rng);
    use std::io::Write;
    let mut picked = std::fs::File::create(&format!("{prefix}.0.fq")).unwrap();
    let mut other = std::fs::File::create(&format!("{prefix}.1.fq")).unwrap();
    let mut picked_so_far = 0;
    for record in records.iter() {
        if picked_so_far < amount {
            writeln!(picked, "{record}")?;
        } else {
            writeln!(other, "{record}")?;
        }
        picked_so_far += record.seq().len();
    }
    Ok(())
}

fn parse_si(input: &str) -> Result<usize, std::num::ParseIntError> {
    assert!(!input.is_empty(), "Please specify genome size.");
    let mut input = input.to_string();
    let last = input.chars().last().unwrap();
    let mult = match last {
        'k' | 'K' => 1_000,
        'm' | 'M' => 1_000_000,
        'g' | 'G' => 1_000_000_000,
        '0'..='9' => 1,
        _ => panic!("si prefix {} is not supported yet.", last),
    };
    if last.is_ascii_alphabetic() {
        input.pop();
    }
    input.parse::<usize>().map(|digit| digit * mult)
}
