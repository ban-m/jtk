extern crate log;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    use std::io::BufReader;
    let rdr = BufReader::new(std::fs::File::open(&args[1])?);
    let dataset: definitions::DataSet = serde_json::de::from_reader(rdr).unwrap();
    use haplotyper::GlobalClustering;
    let config = haplotyper::GlobalClusteringConfig::new(2, 3, 50);
    let dataset = dataset.global_clustering(&config);
    let stdout = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(stdout.lock());
    match serde_json::ser::to_writer_pretty(&mut wtr, &dataset) {
        Err(why) => {
            eprintln!("{:?}", why);
            std::process::exit(1);
        }
        Ok(()) => Ok(()),
    }
}
