#[macro_use]
extern crate log;
const UNIT_ID: u64 = 12;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    use std::io::BufReader;
    let rdr = BufReader::new(std::fs::File::open(&args[1])?);
    let dataset = serde_json::de::from_reader(rdr).unwrap();
    let config = haplotyper::ClusteringConfig::with_default(&dataset, 12, 4, 100, 5000);
    let ref_unit = dataset
        .selected_chunks
        .iter()
        .find(|e| e.id == UNIT_ID)
        .unwrap();
    let units: Vec<_> = dataset
        .encoded_reads
        .iter()
        .filter_map(|r| {
            let (idx, node) = r
                .nodes
                .iter()
                .enumerate()
                .find(|&(_, n)| n.unit == UNIT_ID)?;
            Some((r.id, idx, node))
        })
        .collect();
    // Clustering
    let result = haplotyper::unit_clustering(&units, &config, ref_unit);
    for (asn, (readid, _, _)) in result.iter().zip(units) {
        debug!("{}\t{}", asn, readid);
    }
    Ok(())
}
