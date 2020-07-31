extern crate log;
use haplotyper::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    let args: Vec<_> = std::env::args().collect();
    let rdr = BufReader::new(std::fs::File::open(&args[1])?);
    let mut dataset: definitions::DataSet = serde_json::de::from_reader(rdr).unwrap();
    let chunk_len: usize = 2000;
    let margin: usize = 500;
    let skip_len: usize = 2000;
    rayon::ThreadPoolBuilder::new()
        .num_threads(24)
        .build_global()
        .unwrap();
    let config = UnitConfig::new_clr(chunk_len, 2000, skip_len, margin);
    dataset = dataset.select_chunks(&config);
    Ok(())
}
