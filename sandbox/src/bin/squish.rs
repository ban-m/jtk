#![allow(unused_imports)]
use definitions::*;
use haplotyper::encode::Encode;
use std::{collections::HashSet, io::*};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let ids: HashSet<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    ds.selected_chunks
        .iter_mut()
        .filter(|c| ids.contains(&c.id))
        .for_each(|c| c.cluster_num = 1);
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .filter(|n| ids.contains(&n.unit))
        .for_each(|n| {
            n.cluster = 0;
            n.posterior = vec![0f64];
        });
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}
