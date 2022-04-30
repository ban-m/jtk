#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
use haplotyper::local_clustering::{local_clustering_selected, LocalClustering};
use kiley::recover;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use sandbox::IS_MOCK;
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let cov: f64 = args[2].parse().unwrap();
    ds.coverage = Some(cov);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // let selections: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    // use haplotyper::local_clustering::*;
    // local_clustering_selected(&mut ds, &selections);
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|r| r.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    // use haplotyper::multiplicity_estimation::*;
    // let threads = 56;
    // let multip_config =
    //     MultiplicityEstimationConfig::new(threads, 230493, ds.coverage, Some("multip.gfa"));
    // ds.estimate_multiplicity(&multip_config);
    Ok(())
}
