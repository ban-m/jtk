use definitions::*;
use sandbox::IS_MOCK;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    let ids: HashSet<u64> = args[2..].iter().filter_map(|r| r.parse().ok()).collect();
    for read in ds.encoded_reads.iter().filter(|r| ids.contains(&r.id)) {
        let is_hap1 = match IS_MOCK {
            false => id2desc[&read.id].contains("251v2"),
            true => id2desc[&read.id].contains("hapA"),
        };
        println!("{is_hap1}\t{read}");
    }
    Ok(())
}
