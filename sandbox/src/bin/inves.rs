#![allow(unused_imports)]
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for r in ds.encoded_reads.iter() {
        println!("{}\t{}", r.id, r.original_length);
    }
    Ok(())
}
