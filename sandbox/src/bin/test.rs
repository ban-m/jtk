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
    let mut counts: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .for_each(|n| n.cluster = 0);
    use haplotyper::multiplicity_estimation::*;
    let config = MultiplicityEstimationConfig::new(1, 4390, Some("test.gfa"));
    ds.estimate_multiplicity(&config);
    // use haplotyper::copy_number_estimation_mrf::*;
    // let config = Config::new(ds.coverage.unwrap(), 9482390);
    // ds.update_copy_numbers(&config);
    for c in ds.selected_chunks.iter() {
        let prev = prev[&c.id];
        let count = counts[&c.id];
        println!("{}\t{}\t{}\t{}", c.id, c.cluster_num, prev, count);
    }
    Ok(())
}
