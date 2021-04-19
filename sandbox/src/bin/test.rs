use definitions::*;
use haplotyper::em_correction::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::collections::HashMap;
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let unit = 231;
    let seqs: Vec<_> = ds
        .encoded_reads
        .iter()
        .flat_map(|x| {
            x.nodes
                .iter()
                .filter_map(|n| (n.unit == unit).then(|| n.seq()))
                .collect::<Vec<_>>()
        })
        .collect();
    let start = std::time::Instant::now();
    let _consensus = kiley::consensus(&seqs, 394, 3, 100);
    let end = std::time::Instant::now();
    eprintln!("{}", (end - start).as_secs());
    Ok(())
}
