use definitions::*;
use std::io::*;
// Fill deletion test.
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
    ds = haplotyper::encode::deletion_fill::correct_unit_deletion(ds, &mut failed_trials, 0.35);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}
