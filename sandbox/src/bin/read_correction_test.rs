use definitions::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    let reads_summary: Vec<Vec<(u64, u64)>> = ds
        .encoded_reads
        .iter()
        .map(|read| read.nodes.iter().map(|n| (n.unit, n.cluster)).collect())
        .collect();
    let rev_for_reads: Vec<_> = {
        let rev = reads_summary
            .iter()
            .map(|read| read.iter().copied().rev().collect::<Vec<_>>());
        reads_summary.iter().cloned().zip(rev).collect()
    };
    let c = haplotyper::PolishClusteringConfig::new(1, -1, -2, 0);
    let (read, _) = reads_summary
        .iter()
        .zip(ds.encoded_reads.iter())
        .find(|(_, r)| r.id == 2549)
        .unwrap();
    let before: Vec<_> = read.iter().map(|(u, c)| format!("{}:{}", u, c)).collect();
    let result = haplotyper::correct_read(read, &rev_for_reads, &c);
    let after: Vec<_> = result.iter().map(|(u, c)| format!("{}:{}", u, c)).collect();
    println!("{}\n{}", before.join("-"), after.join("-"));
    Ok(())
}
