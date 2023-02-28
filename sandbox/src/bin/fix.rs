#[allow(unused_macros)]
macro_rules! elapsed {
    ($a:expr) => {{
        let start = std::time::Instant::now();
        let return_value = $a;
        let end = std::time::Instant::now();
        (return_value, (end - start))
    }};
}

use definitions::*;

use std::{collections::HashMap, io::BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let mut edges: HashMap<_, (i64, usize)> = HashMap::new();
    for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
        let key = (edge.from.min(edge.to), edge.from.max(edge.to));
        let entry = edges.entry(key).or_default();
        entry.0 += edge.offset;
        entry.1 += 1;
    }
    for ((from, to), &(total, coverage)) in edges.iter() {
        let average = total / coverage as i64;
        println!("{from}\t{to}\t{average}\t{coverage}");
    }
    Ok(())
}
