#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::ditch_graph::Focus;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use sandbox::IS_MOCK;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::dirichlet_mixture::{
        correct_unit, ClusteringConfig, DirichletMixtureCorrection,
    };
    let config = ClusteringConfig::default();
    for unit in args[2..].iter() {
        let unit: u64 = unit.parse().unwrap();
        let k = ds
            .selected_chunks
            .iter()
            .find(|u| u.id == unit)
            .unwrap()
            .cluster_num;
        correct_unit(&ds, unit, k, &config);
    }
    Ok(())
}
