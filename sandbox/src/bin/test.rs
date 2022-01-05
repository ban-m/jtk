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
    use haplotyper::copy_number_estimation::*;
    let config = Config::default();
    let prev: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    ds.update_copy_numbers(&config);
    for c in ds.selected_chunks.iter() {
        let prev = prev[&c.id];
        println!("{}\t{}\t{}", c.id, c.cluster_num, prev);
    }
    Ok(())
}
