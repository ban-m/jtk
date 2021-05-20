// use bio_utils::fasta::Record;
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let mut unit_count: HashMap<_, u32> = ds.selected_chunks.iter().map(|u| (u.id, 0)).collect();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *unit_count.entry(node.unit).or_default() += 1;
        }
    }
    for unit in ds.selected_chunks.iter() {
        println!(
            "{}\t{}\t{}",
            unit.id, unit_count[&unit.id], unit.cluster_num
        );
    }
}
