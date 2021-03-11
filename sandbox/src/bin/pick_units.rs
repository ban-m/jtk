// use bio_utils::fasta::Record;
use definitions::*;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufReader;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let unit: HashSet<u64> = vec![487].into_iter().collect();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.desc.clone()))
        .collect();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter().filter(|n| unit.contains(&n.unit)) {
            println!("{}\t{}\t{}", id2desc[&read.id], node.unit, node.cluster);
        }
    }
}
