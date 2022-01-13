#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
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
    let selection: HashSet<_> = vec![1002].into_iter().collect();
    ds.selected_chunks
        .iter_mut()
        .find(|c| c.id == 1002)
        .unwrap()
        .cluster_num = 10;
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &selection);
    // use haplotyper::copy_number_estimation_mrf::*;
    // let config = Config::new(ds.coverage.unwrap(), 9482390);
    // ds.update_copy_numbers(&config);
    // for c in ds.selected_chunks.iter() {
    //     let prev = prev[&c.id];
    //     let count = counts[&c.id];
    //     println!("{}\t{}\t{}\t{}", c.id, c.cluster_num, prev, count);
    // }
    Ok(())
}
