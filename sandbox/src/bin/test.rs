use definitions::*;
#[macro_use]
extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let unit: u64 = args[2].parse().unwrap();
    debug!("Started");
    let mut dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    // let mut c = haplotyper::ClusteringConfig::clr(&dataset, 2, 100, 1000);
    let mut c = haplotyper::ClusteringConfig::ccs(&dataset, 2, 100, 1000);
    c.limit = 3000;
    c.seed = 100;
    let ref_unit = dataset
        .selected_chunks
        .iter()
        .find(|u| u.id == unit)
        .unwrap()
        .clone();
    use std::collections::HashMap;
    let id2name: HashMap<_, _> = dataset
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    let mut units: Vec<_> = dataset
        .encoded_reads
        .iter_mut()
        .filter_map(|read| {
            let id = read.id;
            read.nodes
                .iter_mut()
                .enumerate()
                .find(|(_, n)| n.unit == unit)
                .map(|(idx, n)| (id, idx, n))
        })
        .enumerate()
        .map(|(_, x)| x)
        .collect();
    Ok(())
}
