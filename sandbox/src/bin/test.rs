use definitions::*;
// use log::*;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let mut ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    if args.len() > 2 {
        let unit: u64 = args[2].parse().unwrap();
        let is_hap_a: HashMap<_, _> = ds
            .raw_reads
            .iter()
            .map(|read| (read.id, if read.name.contains("hapA") { 1 } else { 0 }))
            .collect();
        let reads: Vec<_> = ds
            .encoded_reads
            .iter()
            .filter(|r| r.nodes.iter().any(|n| n.unit == unit))
            .collect();
        let config = haplotyper::em_correction::Config::with_default(unit, 2);
        let preds = haplotyper::em_correction::clustering(&ds.selected_chunks, &reads, &config);
        for cl in 0..2 {
            for (pred, read) in preds.iter().zip(reads.iter()).filter(|&(&p, _)| p == cl) {
                let line: Vec<_> = read
                    .nodes
                    .iter()
                    .map(|n| format!("{:2}-{:2}", n.unit, n.cluster))
                    .collect();
                println!("{}\t{}\t{}", is_hap_a[&read.id], pred, line.join(" "));
            }
        }
    } else {
        let units: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, true)).collect();
        haplotyper::em_correction::impute_clustering(&mut ds, &units);
        println!("{}", serde_json::ser::to_string(&ds).unwrap());
    }
    Ok(())
}
