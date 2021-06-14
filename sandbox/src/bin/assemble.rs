use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    ds.assignments = ds
        .encoded_reads
        .iter()
        .map(|r| Assignment::new(r.id, 0))
        .collect();
    let config = haplotyper::AssembleConfig::new(24, 2000, false);
    use haplotyper::Assemble;
    let gfa = ds.assemble_as_gfa(&config);
    // let coverage = ds.coverage.unwrap();
    // let lens: Vec<_> = ds.raw_reads.iter().map(|read| read.seq().len()).collect();
    // haplotyper::assemble::copy_number::estimate_copy_number_on_gfa(&mut gfa, coverage, &lens, 2000);
    println!("{}", gfa);
    Ok(())
}
