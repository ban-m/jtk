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
    //let selection: HashSet<_> = vec![1139].into_iter().collect();
    let selection: HashSet<_> = vec![1623].into_iter().collect();
    use haplotyper::local_clustering::LocalClustering;
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &selection);
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        let ans = match IS_MOCK {
            true => id2desc[&read.id].contains("hapA") as usize,
            false => id2desc[&read.id].contains("000251v2") as usize,
        };
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    let score: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u.score)).collect();
    println!("UNIT\tunit\tcluster\thap1\thap2\tpurity\tscore");
    for ((unit, cluster), counts) in counts.iter().filter(|x| selection.contains(&x.0 .0)) {
        let score = score[unit];
        let total = counts[0] + counts[1];
        let pur = counts[0].max(counts[1]) as f64 / total as f64;
        println!(
            "UNIT\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            unit, cluster, counts[0], counts[1], pur, score,
        );
    }
    Ok(())
}
