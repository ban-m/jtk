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
    // haplotyper::encode::deletion_fill::remove_weak_edges(&mut ds);
    // let mut rng: Xoshiro256Plus = SeedableRng::seed_from_u64(203);
    // use haplotyper::local_clustering::kmeans::ClusteringConfig;
    // let config = ClusteringConfig::new(100, 2, ds.coverage.unwrap(), 2.5, ds.read_type);
    // for unit in [1203, 705] {
    //     for cluster in 0..2 {
    //         let nodes: Vec<_> = ds
    //             .encoded_reads
    //             .iter()
    //             .flat_map(|r| r.nodes.iter())
    //             .filter(|n| (n.unit, n.cluster) == (unit, cluster))
    //             .map(|n| n.seq())
    //             .collect();
    //         haplotyper::local_clustering::kmeans::clustering(&nodes, &mut rng, &config);
    //     }
    // }
    let selection: HashSet<u64> = args[2..].iter().filter_map(|r| r.parse().ok()).collect();
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &selection);
    // use haplotyper::phmm_likelihood_correction::*;
    // let config = CorrectionConfig::default();
    // ds.correct_clustering_selected(&selection, &config);

    Ok(())
}
