/// Synopsis: [fasta file] [into N chunks] [with prefix]
use bio_utils::fastq::*;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    // First argument would be the name of the binary...
    let args: Vec<_> = std::env::args().collect();
    let mut records = parse_into_vec(&args[1]).unwrap();
    let split_into: usize = args[2].parse().unwrap();
    let prefix = &args[3];
    let mut files: Vec<_> = (0..split_into)
        .map(|i| std::fs::File::create(&format!("{prefix}.{}.fq", i)).unwrap())
        .collect();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(34204);
    records.shuffle(&mut rng);
    use std::io::Write;
    for (i, record) in records.iter().enumerate() {
        writeln!(files[i % split_into], "{record}")?;
    }
    Ok(())
}
