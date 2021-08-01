use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    use haplotyper::Assemble;
    ds.squish_small_contig(&config, 25);
    let gfa = ds.assemble(&config);
    println!("{}", gfa);
    Ok(())
}
