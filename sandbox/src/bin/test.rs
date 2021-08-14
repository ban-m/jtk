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
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    let (nodes, ids): (Vec<_>, Vec<_>) = ds
        .encoded_reads
        .iter()
        .flat_map(|r| {
            r.nodes
                .iter()
                .filter(|n| n.unit == 153)
                .map(|n| (n.seq(), r.id))
                .collect::<Vec<_>>()
        })
        .unzip();
    use haplotyper::local_clustering;
    let mut config = local_clustering::kmeans::ClusteringConfig::new(100, 2, ds.coverage.unwrap());
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128Plus;
    let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(293);
    let (asn, _) = local_clustering::kmeans::clustering(&nodes, &mut rng, &mut config).unwrap();
    for ((asn, id), seq) in asn.iter().zip(ids.iter()).zip(nodes.iter()) {
        let seq = String::from_utf8_lossy(seq);
        println!(
            "{}\t{}\t{}\t{}",
            id2desc[id],
            id2desc[id].contains("252v"),
            asn,
            seq
        );
    }
    Ok(())
}
