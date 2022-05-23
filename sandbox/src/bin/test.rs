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
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let reads: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
    let read = ds.encoded_reads.iter().find(|r| r.id == 4434).unwrap();
    use haplotyper::encode::deletion_fill::{get_pileup, ReadSkelton};
    let seq: Vec<_> = read
        .nodes
        .iter()
        .map(|n| format!("{}-{}", n.unit, n.cluster))
        .collect();
    println!("REF\t{}", seq.join("\t"));
    let pileup = get_pileup(read, &reads);
    let nodes = &read.nodes;
    for (i, (n, p)) in read.nodes.iter().zip(pileup.iter()).enumerate() {
        println!("{i}\t{}\t{}\t{p:?}", n.unit, n.cluster);
        let head_cand = p.check_insertion_head(nodes, 2, i);
        println!("{head_cand:?}");
    }
    Ok(())
}
