use definitions::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.name.clone()))
        .collect();
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(3214);
    log::debug!("Start");
    use rayon::prelude::*;
    let result: Vec<_> = ds
        .selected_chunks
        .par_iter()
        // .filter(|u| vec![287, 61].contains(&u.id))
        .map(|unit| {
            let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(unit.id * 23);
            let (seqs, answer): (Vec<_>, Vec<_>) = ds
                .encoded_reads
                .iter()
                .flat_map(|read| {
                    let answer = id2desc[&read.id].starts_with("hapA");
                    read.nodes
                        .iter()
                        .filter(|n| n.unit == unit.id)
                        .map(|n| (n.seq(), answer as u8))
                        .collect::<Vec<_>>()
                })
                .unzip();
            let config =
                haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, 20, 5, 2, 25);
            let start = std::time::Instant::now();
            let asn =
                haplotyper::local_clustering::kmeans::clustering(&seqs, &mut rng, &config).unwrap();
            let km = haplotyper::local_clustering::rand_index(&answer, &asn);
            let end = std::time::Instant::now();
            log::debug!("FIN\t{}\t{}", unit.id, (end - start).as_millis());
            let acc = asn
                .iter()
                .zip(answer.iter())
                .filter(|(x, y)| x == y)
                .count();
            let acc = acc as f64 / asn.len() as f64;
            (unit.id, km, acc.max(1f64 - acc))
        })
        .collect();
    log::debug!("End");
    for (id, acc, score) in result.iter() {
        println!("{}\t{}\t{}", id, acc, score);
    }
    Ok(())
}
