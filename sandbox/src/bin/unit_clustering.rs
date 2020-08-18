use definitions::*;
#[macro_use]
extern crate log;
const UNIT: u64 = 30;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let mut dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let c = haplotyper::ClusteringConfig::clr(&dataset, 2, 100, 1000);
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
    haplotyper::unit_clustering(&mut units, &c, &ref_unit);
    debug!("Clustering ended.");
    let mut result: Vec<_> = units
        .iter()
        .map(|(id, _, unit)| (&id2name[&id], unit.cluster))
        .collect();
    result.sort_by_key(|x| x.1);
    for (i, x) in result {
        eprintln!("{}\t{}", x, i);
    }
    Ok(())
}
