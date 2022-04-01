use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let ids: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    for read in ds.encoded_reads.iter().filter(|r| ids.contains(&r.id)) {
        println!("{}", read);
    }
    Ok(())
}
