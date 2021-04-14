use definitions::*;
use haplotyper::em_correction::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    ds = ds.correct_clustering_em(10, 3);
    use haplotyper::assemble::*;
    ds.assignments = ds
        .raw_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    let config = AssembleConfig::new(24, 100, false);
    println!("{}", ds.assemble_as_gfa(&config));
    Ok(())
}
