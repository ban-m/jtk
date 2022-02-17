#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
use log::debug;
use rand::SeedableRng;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let config = AssembleConfig::new(1, 100, false, false, 4);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    graph.remove_lightweight_edges(2, true);
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    let cov = ds.coverage.unwrap();
    println!("Type\tUnit\tCopyNumber\tCoverage");
    let counts = {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *counts.entry(node.unit).or_default() += 1;
        }
        counts
    };
    // let (copy, _) = graph.copy_number_estimation(cov, &lens);
    // for ((unit, _), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("Grad\t{}\t{}\t{}", unit, cp, cov);
    // }
    // let (copy, _) = graph.copy_number_estimation_gbs(cov, &lens);
    // for ((unit, _), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("Gibbs\t{}\t{}\t{}", unit, cp, cov);
    // }
    debug!("START");
    let (copy, _) = graph.copy_number_estimation_mcmc(cov, &lens);
    for ((unit, _), cp) in copy {
        let cov = counts[&unit];
        println!("MCMC\t{}\t{}\t{}", unit, cp, cov);
    }
    Ok(())
}
