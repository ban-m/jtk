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
    for r in ds.raw_reads.iter() {
        println!("{}\t{}\t{}", r.id, r.name, r.seq().len());
    }
    // let ids: Vec<_> = args[2..]
    //     .iter()
    //     .map(|x| x.parse::<u64>().unwrap())
    //     .collect();
    // for read in ds.encoded_reads.iter().filter(|r| ids.contains(&r.id)) {
    //     let mut seq = read.recover_raw_read();
    //     seq.iter_mut().for_each(u8::make_ascii_uppercase);
    //     let id = read.id;
    //     let seq = std::str::from_utf8(&seq).unwrap();
    //     println!(">{}\n{}", id, seq);
    // }
    Ok(())
}
