use definitions::*;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    let unit: u64 = args[2].parse().unwrap();
    let cluster: Option<u64> = args[3].parse().ok();
    let is_hap_a: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, if read.name.contains("hapA") { 1 } else { 0 }))
        .collect();
    for read in ds.encoded_reads.iter() {
        if read
            .nodes
            .iter()
            .any(|n| n.unit == unit && (Some(n.cluster) == cluster || cluster.is_none()))
        {
            let line: Vec<_> = read
                .nodes
                .iter()
                .map(|n| format!("{}-{}", n.unit, n.cluster))
                .collect();
            println!("{}\t{}\t{}", read.id, is_hap_a[&read.id], line.join(":"));
        }
    }
    Ok(())
}
