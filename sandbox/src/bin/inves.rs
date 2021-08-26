use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.desc.contains("252v2")))
        .collect();
    let mut counts = vec![[0; 2], [0; 2]];
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter().filter(|n| n.unit == 296) {
            counts[node.cluster as usize][id2desc[&read.id] as usize] += 1;
        }
    }
    for i in 0..2 {
        println!("{:?}", counts[i]);
    }
    Ok(())
}
