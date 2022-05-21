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
    let nodes: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|r| r.nodes.iter().filter(|n| n.unit == 659))
        .collect();
    let cps = &[2f64, 1f64, 1f64];
    for (i, u) in nodes.iter().enumerate() {
        for v in nodes.iter().skip(i + 1) {
            let u_post = &u.posterior;
            let v_post = &v.posterior;
            let center = haplotyper::phmm_likelihood_correction::sim(u_post, v_post, cps);
            let line1: Vec<_> = u_post.iter().map(|x| format!("{x:.2}")).collect();
            let line2: Vec<_> = v_post.iter().map(|x| format!("{x:.2}")).collect();
            println!("[{}]\t[{}]\t{center:.3}", line1.join(","), line2.join(","));
        }
    }

    Ok(())
}
