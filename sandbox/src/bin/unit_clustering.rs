use definitions::*;
#[macro_use]
extern crate log;
const UNIT: u64 = 500;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("trace")).init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let mut dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    debug!("Configuring...");
    let config = haplotyper::ClusteringConfig::with_default(&dataset, 1, 3, 100, 1000);
    debug!("Configured.");
    let ref_unit = dataset
        .selected_chunks
        .iter()
        .find(|u| u.id == UNIT)
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
                .find(|(_, n)| n.unit == UNIT)
                .map(|(idx, n)| (id, idx, n))
        })
        .collect();
    debug!("Clustering started.");
    haplotyper::unit_clustering(&mut units, &config, &ref_unit);
    debug!("Clustering ended.");
    for (id, _, unit) in units {
        eprintln!("{}\t{}", id2name[&id], unit.cluster);
    }
    Ok(())
}
