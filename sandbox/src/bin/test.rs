#![allow(unused_imports)]
use definitions::*;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let prev: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    use haplotyper::multiplicity_estimation::*;
    let config = MultiplicityEstimationConfig::new(24, 2234, None);
    ds.estimate_multiplicity(&config);
    let after: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    let mut counts: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    println!("Unit\tPrev\tAfter\tCoverage");
    for (unit, prev) in prev.iter() {
        let after = after[unit];
        let cov = counts[unit];
        println!("{}\t{}\t{}\t{}", unit, prev, after, cov);
    }
    Ok(())
}
