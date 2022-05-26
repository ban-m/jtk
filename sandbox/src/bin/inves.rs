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
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let selection: HashSet<_> = vec![188].into_iter().collect();
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &selection);
    Ok(())
}
