use definitions::*;
use serde_json;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let remove_unit = vec![376, 1648];
    for read in ds.encoded_reads.iter_mut() {
        while let Some(idx) = read
            .nodes
            .iter()
            .position(|n| remove_unit.contains(&n.unit))
        {
            read.remove(idx);
        }
    }
    Ok(())
}
