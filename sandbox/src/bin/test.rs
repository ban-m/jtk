use definitions::*;
use std::collections::HashSet;
use std::io::*;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let selection: HashSet<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    use haplotyper::phmm_likelihood_correction::*;
    let config = CorrectionConfig::default();
    // rayon::ThreadPoolBuilder::new()
    //     .num_threads(1)
    //     .build()
    //     .unwrap();
    ds.correct_clustering_selected(&selection, &config);
    Ok(())
}
