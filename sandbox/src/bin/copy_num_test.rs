#![allow(unused_imports)]
use definitions::*;
use haplotyper::determine_units::DetermineUnit;
use haplotyper::encode::Encode;
// use haplotyper::DetermineUnit;
use haplotyper::assemble::*;
use log::debug;
use rand::SeedableRng;
use rand_xoshiro::{Seed512, Xoroshiro128PlusPlus, Xoshiro256Plus};
use std::collections::{HashMap, HashSet};
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    // use haplotyper::assemble::copy_number::estimate_copy_number_mcmc;
    // let nodes = vec![
    //     (40f64, 10),
    //     (20f64, 10),
    //     (20f64, 10),
    //     (20f64, 10),
    //     (20f64, 30),
    //     (40f64, 5),
    // ];
    // let edges = vec![
    //     (0, true, 1, false, 14f64),
    //     (0, true, 4, false, 13f64),
    //     (1, true, 2, false, 20f64),
    //     (2, true, 2, false, 3f64),
    //     (2, true, 3, false, 13f64),
    //     (3, true, 5, false, 22f64),
    //     (4, true, 5, false, 17f64),
    // ];
    // let coverage = 20f64;
    // let (nodes_cp, edges_cp) = estimate_copy_number_mcmc(&nodes, &edges, coverage);
    // eprintln!("NODE\t{:?}", nodes_cp);
    // eprintln!("EDGES\t{:?}", edges_cp);
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    // ds.encoded_reads
    //     .iter_mut()
    //     .flat_map(|r| r.nodes.iter_mut())
    //     .for_each(|n| n.cluster = 0);
    let config = AssembleConfig::new(1, 100, false, false, 4, 4f64);
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use haplotyper::assemble::ditch_graph::DitchGraph;
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), ds.read_type, &config);
    graph.remove_lightweight_edges(2, true);
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    let cov = ds.coverage.unwrap();
    // println!("Type\tUnit\tCluster\tCopyNumber\tCoverage");
    // let counts = {
    //     let mut counts: HashMap<_, u32> = HashMap::new();
    //     for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
    //         *counts.entry(node.unit).or_default() += 1;
    //     }
    //     counts
    // };
    // let (copy, _) = graph.copy_number_estimation(cov, &lens);
    // for ((unit, cl), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("Grad\t{unit}\t{cl}\t{cp}\t{cov}");
    // }
    // let (copy, _) = graph.copy_number_estimation_gbs(cov, &lens);
    // for ((unit, cl), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("Gibbs\t{unit}\t{cl}\t{cp}\t{cov}");
    // }
    // let (copy, _) = graph.copy_number_estimation_mcmc(cov, &lens);
    // for ((unit, cl), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("MCMC\t{unit}\t{cl}\t{cp}\t{cov}");
    // }
    // let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(3429);
    // let (copy, _) = graph.copy_number_estimation_mst(cov, &mut rng);
    // for ((unit, cl), cp) in copy {
    //     let cov = counts[&unit];
    //     println!("MST\t{unit}\t{cl}\t{cp}\t{cov}");
    // }
    use std::time::Instant;
    let start = Instant::now();
    let (_, copy_mcmc) = graph.copy_number_estimation_mcmc(cov, &lens);
    let mcmc_end = Instant::now();
    let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(3429);
    let (_, copy_mst) = graph.copy_number_estimation_mst(cov, &mut rng);
    let mst_end = Instant::now();
    for edge in graph.edges() {
        let key = edge.key();
        let occ = edge.occ();
        let mcmc = copy_mcmc[&key];
        let mst = copy_mst[&key];
        let (((fc, fu), _), ((tc, tu), _)) = key;
        println!("{fc}\t{fu}\t{tc}\t{tu}\t{mcmc}\t{mst}\t{occ}");
    }
    eprintln!(
        "{}\t{}",
        (mcmc_end - start).as_millis(),
        (mst_end - mcmc_end).as_millis()
    );
    Ok(())
}
