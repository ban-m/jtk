use definitions::*;
use serde_json;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    for read in ds.encoded_reads.iter() {
        println!("============={}==========", read.id);
        for (n, e) in read.nodes.iter().zip(read.edges.iter()) {
            println!("{}\t{}", n, e);
        }
    }
    Ok(())
}
