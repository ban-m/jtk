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
    use haplotyper::encode::deletion_fill::correct_unit_deletion;
    let ds = correct_unit_deletion(ds);
    use std::collections::HashMap;
    let mut edges: HashMap<_, u32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for w in read.nodes.windows(2) {
            let tuple = (w[0].unit.min(w[1].unit), w[0].unit.max(w[1].unit));
            *edges.entry(tuple).or_default() += 1;
        }
    }
    for (key, val) in edges.iter() {
        if val < &5 {
            println!("{:?}\t{}", key, val);
        }
    }
    for read in ds.encoded_reads.iter() {
        let has_light_weight = read.nodes.windows(2).any(|w| {
            let tuple = (w[0].unit.min(w[1].unit), w[0].unit.max(w[1].unit));
            // tuple == (35, 393)
            edges[&tuple] < 5
        });
        if has_light_weight {
            println!("{}", read);
        }
    }
    Ok(())
}
