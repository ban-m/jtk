use definitions::DataSet;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let _ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    let _unit: u64 = args[2].parse().unwrap();
    Ok(())
}
