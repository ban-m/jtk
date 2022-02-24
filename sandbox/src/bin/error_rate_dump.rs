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
use std::collections::{HashMap, HashSet};
use std::io::*;
fn error(node: &definitions::Node, ref_unit: &Unit) -> (f64, f64, f64) {
    let (query, aln, refr) = node.recover(ref_unit);
    let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
    let del = query.iter().filter(|&&x| x == b' ').count() as f64;
    let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
    let aln_len = aln.len() as f64;
    (mismat / aln_len, del / aln_len, ins / aln_len)
}
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    use rayon::prelude::*;
    let error_rates: Vec<_> = ds
        .encoded_reads
        .par_iter()
        .flat_map(|r| {
            r.nodes
                .par_iter()
                .map(|n| (r.id, n.unit, n.cluster, error(n, ref_units[&n.unit])))
                .collect::<Vec<_>>()
        })
        .collect();
    println!("readid\tunit\tcluster\tmism\tins\tdel");
    for (rid, unit, cluster, (mism, ins, del)) in error_rates {
        println!("{}\t{}\t{}\t{}\t{}\t{}", rid, unit, cluster, mism, ins, del);
    }
    Ok(())
}
