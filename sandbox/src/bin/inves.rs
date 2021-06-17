use definitions::*;
use serde_json;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let unit: u64 = args[2].parse().unwrap();
    let cluster: u64 = args[3].parse().unwrap();
    let (ids, units): (Vec<_>, Vec<_>) = ds
        .encoded_reads
        .iter()
        .flat_map(|r| {
            r.nodes
                .iter()
                .filter(|n| n.unit == unit) //  && n.cluster == cluster)
                .map(|n| (r.id, n.seq()))
                .collect::<Vec<_>>()
        })
        .unzip();
    use rand_xoshiro::Xoroshiro128PlusPlus;
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(42);
    // let cov = ds.coverage.unwrap();
    let cov = 27f64;
    let mut config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, 4, cov);
    let (asn, _) =
        haplotyper::local_clustering::kmeans::clustering(&units, &mut rng, &mut config).unwrap();
    let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    for (id, x) in ids.iter().zip(asn.iter()) {
        println!("{}\t{}", id2desc[id], x);
    }
    Ok(())
}
