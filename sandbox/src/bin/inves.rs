#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::ditch_graph::Focus;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use nalgebra::RealField;
use sandbox::IS_MOCK;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let ids: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    for read in ds.encoded_reads.iter() {
        if ids.contains(&read.id) {
            let seq = read.recover_raw_read();
            println!(">{}\n{}", read.id, std::str::from_utf8(&seq).unwrap());
        }
    }
    // let units: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    // let num_cluster: HashMap<_, _> = ds
    //     .selected_chunks
    //     .iter()
    //     .map(|c| (c.id, c.cluster_num))
    //     .collect();
    // for (i, &unit1) in units.iter().enumerate() {
    //     for &unit2 in units.iter().skip(i + 1) {
    //         let mut occs = vec![vec![0; num_cluster[&unit2]]; num_cluster[&unit1]];
    //         for read in ds.encoded_reads.iter() {
    //             for node1 in read.nodes.iter().filter(|n| n.unit == unit1) {
    //                 for node2 in read.nodes.iter().filter(|n| n.unit == unit2) {
    //                     occs[node1.cluster as usize][node2.cluster as usize] += 1;
    //                 }
    //             }
    //         }
    //         if occs.iter().flatten().sum::<u32>() == 0 {
    //             continue;
    //         }
    //         println!("{unit1}\t{unit2}");
    //         for row in occs.iter() {
    //             let row: Vec<_> = row.iter().map(|x| format!("{x}")).collect();
    //             println!("\t{}", row.join("\t"));
    //         }
    //     }
    // }
    Ok(())
}
