#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
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
    use haplotyper::dense_encoding::*;
    let config = DenseEncodingConfig::new(10);
    let prev: usize = ds.encoded_reads.iter().map(|r| r.nodes.len()).sum();
    ds.dense_encoding_dev(&config);
    let now: usize = ds.encoded_reads.iter().map(|r| r.nodes.len()).sum();
    log::debug!("{}->{}", prev, now);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    // let id2desc: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.name.clone()))
    //     .collect();
    // let mut nodes: Vec<_> = vec![];
    // let mut answer = vec![];
    // // let target = (210, 0);
    // for read in ds.encoded_reads.iter() {
    //     let is_hap1 = id2desc[&read.id].contains("000251v2") as usize;
    //     // for node in read.nodes.iter().filter(|n| (n.unit, n.cluster) == target) {
    //     for node in read.nodes.iter().filter(|n| n.unit == 210) {
    //         nodes.push(node.seq());
    //         answer.push(is_hap1);
    //     }
    // }
    // let mut rng: rand_xoshiro::Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(4329890);
    // use haplotyper::local_clustering::kmeans::ClusteringConfig;
    // let config = ClusteringConfig::new(100, 4, 25f64, definitions::ReadType::CLR);
    // let (preds, gains, _, _) =
    //     haplotyper::local_clustering::kmeans::clustering(&nodes, &mut rng, &config).unwrap();
    // for ((p, a), g) in preds.iter().zip(answer.iter()).zip(gains.iter()) {
    //     println!("{}\t{}\t{}", a, p, vec2str(g));
    // }
    Ok(())
}

// fn vec2str(xs: &[f64]) -> String {
//     let xs: Vec<_> = xs.iter().map(|&x| format!("{:6.1}", x)).collect();
//     xs.join(",")
// }
