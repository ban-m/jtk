use definitions::*;
use std::io::*;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let hmm = haplotyper::model_tune::get_model(&ds).unwrap();
    let gain = haplotyper::likelihood_gains::estimate_gain(&hmm, 238019, 100, 20, 5);
    println!("{gain}");
    use std::collections::HashSet;
    let selection: HashSet<u64> = args[2..].iter().filter_map(|x| x.parse().ok()).collect();
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &selection);
    Ok(())
}
