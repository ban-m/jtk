use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    use haplotyper::encode::deletion_fill::estimate_error_rate_dev;
    let (reads, units, _) = estimate_error_rate_dev(&ds, 0.35);
    for (i, read) in reads.iter().enumerate() {
        println!("Read\t{i}\t{read}");
    }
    for (i, unit) in units.iter().enumerate() {
        for (cl, x) in unit.iter().enumerate() {
            println!("Unit\t{i}\t{cl}\t{x}");
        }
    }
    Ok(())
}
