use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let end = std::time::Instant::now();
    eprintln!("Open\t{}", (end - start).as_secs());
    use std::collections::HashMap;
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();

    // use haplotyper::em_correction::ClusteringCorrection;
    // let config = haplotyper::AssembleConfig::new(24, 2000, false, true);
    // let ds = ds.correct_clustering_em(5, 5, 3);

    let mut counts: HashMap<(u64, u64), [u32; 2]> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        // let ans = id2desc[&read.id].contains("hapA") as usize;
        let ans = id2desc[&read.id].contains("000252v2") as usize;
        for node in read.nodes.iter() {
            counts.entry((node.unit, node.cluster)).or_default()[ans] += 1;
        }
    }
    for ((unit, cluster), counts) in counts.iter() {
        println!("{}\t{}\t{}\t{}", unit, cluster, counts[0], counts[1]);
    }
    Ok(())
}
