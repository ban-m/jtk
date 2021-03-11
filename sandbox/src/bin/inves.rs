use definitions::*;
use serde_json;
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let units: HashSet<u64> = vec![1634].into_iter().collect();
    let id2name: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.desc.to_string()))
        .collect();
    for read in ds.encoded_reads.iter() {
        if read.nodes.iter().any(|n| units.contains(&n.unit)) {
            let line: Vec<_> = read.nodes.iter().map(|n| format!("{}", n.unit)).collect();
            println!("{}\t{}", id2name[&read.id], line.join(":"),);
        }
        for node in read.nodes.iter().filter(|n| units.contains(&n.unit)) {
            let ref_unit = ds
                .selected_chunks
                .iter()
                .find(|unit| unit.id == node.unit)
                .unwrap();
            let (r, op, q) = node.recover(&ref_unit);
            for ((r, op), q) in r.chunks(200).zip(op.chunks(200)).zip(q.chunks(200)) {
                println!("{}", String::from_utf8_lossy(r));
                println!("{}", String::from_utf8_lossy(op));
                println!("{}", String::from_utf8_lossy(q));
                println!();
            }
            println!();
            println!();
        }
    }
    Ok(())
}
