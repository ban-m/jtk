use definitions::*;
use std::io::*;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::phmm_likelihood_correction::*;
    let selection: std::collections::HashSet<u64> =
        args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    let config = CorrectionConfig::default();
    ds.correct_clustering_selected(&selection, &config);
    Ok(())
}
