#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::local_clustering::{local_clustering_selected, LocalClustering};
use haplotyper::{assemble::*, local_clustering};
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
    let config = haplotyper::polish_units::PolishUnitConfig::new(ds.read_type, 1, 100);
    use haplotyper::polish_units::PolishUnit;
    ds.polish_unit(&config);
    for c in ds.selected_chunks.iter() {
        println!(">{}\n{}", c.id, std::str::from_utf8(c.seq()).unwrap());
    }
    // let selection: HashSet<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(1)
    //     .build()
    //     .unwrap();
    // local_clustering_selected(&mut ds, &selection);
    // println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}
