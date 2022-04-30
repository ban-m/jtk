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
    let mut nodes = vec![];
    let mut naive_cov = 0f64;
    let mut edges = vec![];
    for (i, line) in std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|l| l.ok())
        .enumerate()
    {
        let fields: Vec<_> = line.split('\t').skip(1).collect();
        if i == 0 {
            naive_cov = fields[0].parse().unwrap();
        } else if fields[0] == "NODE" {
            let cov = fields[1].parse::<f64>().unwrap();
            let len = fields[2].parse::<usize>().unwrap();
            nodes.push((cov, len));
        } else if fields[0] == "EDGE" {
            let from = fields[1].parse::<usize>().unwrap();
            let fp = fields[2] == "true";
            let to = fields[3].parse::<usize>().unwrap();
            let tp = fields[4] == "true";
            let cov = fields[5].parse::<f64>().unwrap();
            edges.push((from, fp, to, tp, cov));
        }
    }
    use haplotyper::assemble::copy_number;
    let (node_cp, edge_cp) = copy_number::estimate_copy_number_mcmc(&nodes, &edges, naive_cov);
    eprintln!("======================");
    let pot = copy_number::get_potential(&nodes, &node_cp, &edges, &edge_cp, naive_cov);
    println!("COVCP\t{naive_cov}");
    for ((cov, len), cp) in nodes.iter().zip(node_cp.iter()) {
        println!("COVCP\tNODE\t{}\t{}\t{}", cov, len, cp);
    }
    for ((f, fp, t, tp, cov), cp) in edges.iter().zip(edge_cp.iter()) {
        println!("COVCP\tEDGE\t{f}\t{fp}\t{t}\t{tp}\t{cov}\t{cp}");
    }
    println!("GET\t{pot}");
    let node_cp = vec![2, 1, 2, 1, 2, 1, 1, 1, 1];
    assert_eq!(node_cp.len(), nodes.len());
    let edge_cp = vec![1; edges.len()];
    let pot = copy_number::get_potential(&nodes, &node_cp, &edges, &edge_cp, naive_cov);
    println!("OPTIM\t{pot}");
    // let ds: DataSet =
    //     serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
    //         .unwrap();
    // let read_ids: Vec<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    // for read in ds.encoded_reads.iter().filter(|r| read_ids.contains(&r.id)) {
    //     println!("{read}");
    // }
    // let min_span_read = ds.read_type.min_span_reads();
    // let llr = ds.read_type.min_llr_value();
    // let config = AssembleConfig::new(1, 1000, false, true, min_span_read, llr);
    // check_foci(&ds, &config);
    Ok(())
}
