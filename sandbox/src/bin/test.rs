use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use haplotyper::assemble::AssembleConfig;
    use haplotyper::Assemble;
    let config = AssembleConfig::new(2, 100, false, true);
    // ds.squish_small_contig(&config);
    let gfa = ds.assemble(&config);
    println!("{}", gfa);
    Ok(())
}
