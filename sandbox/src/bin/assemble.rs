use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let config = haplotyper::AssembleConfig::new(24, 2000, false, false);
    use haplotyper::Assemble;
    let gfa = ds.assemble(&config);
    println!("{}", gfa);
    Ok(())
}
