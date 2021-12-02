use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    use haplotyper::dirichlet_correction::*;
    let config = Config::default();
    ds.correct_clustering(&config);
    println!("{}", serde_json::ser::to_string_pretty(&ds).unwrap());
    Ok(())
}
