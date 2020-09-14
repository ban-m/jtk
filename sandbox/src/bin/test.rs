use definitions::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    use haplotyper::GlobalClustering;
    let config = haplotyper::GlobalClusteringConfig::new(3, 10, 1, -1, -2);
    ds.global_clustering(&config);
    Ok(())
}
