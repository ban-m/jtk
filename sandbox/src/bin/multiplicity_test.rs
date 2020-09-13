use definitions::*;
const UNIT: u64 = 87;
fn main() {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1]).unwrap());
    let ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    for read in ds.encoded_reads.iter() {
        if read.nodes.iter().any(|n| n.unit == UNIT) {
            let line: Vec<_> = read.nodes.iter().map(|n| format!("{}", n.unit)).collect();
            println!("{}", line.join("-"));
        }
    }
}
