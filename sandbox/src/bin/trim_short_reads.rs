#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
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
    ds.encoded_reads.retain(|read| read.original_length > 4_000);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}
