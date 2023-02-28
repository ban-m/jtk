use definitions::DataSet;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::local_clustering::*;
    use std::collections::HashSet;
    let selection: HashSet<_> = args[2..]
        .iter()
        .map(|x| x.parse::<u64>().unwrap())
        .collect();
    ds.local_clustering_selected(&selection);
    Ok(())
}
