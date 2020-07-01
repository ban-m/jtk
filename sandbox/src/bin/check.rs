use definitions::*;
use std::collections::HashMap;
#[macro_use]
extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("trace")).init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let id2name: HashMap<_, _> = dataset
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    let reads: Vec<_> = dataset
        .encoded_reads
        .into_iter()
        .filter(|r| id2name[&r.id].starts_with("hapA"))
        .map(|mut r| {
            r.nodes.iter_mut().for_each(|n| n.cluster = 0);
            r
        })
        .collect();
    let graph = haplotyper::DeBruijnGraph::from_encoded_reads(&reads, 3);
    let components = graph.clustering(0);
    eprintln!("{}", components.len());
    Ok(())
}
