use definitions::*;
use serde_json;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    eprintln!("{:?}", std::time::Instant::now() - start);
    for u in ds.selected_chunks.iter() {
        println!("{}\t{}", u.id, u.cluster_num);
    }
    Ok(())
}
