use definitions::*;
use std::io::BufReader;
// use rand::{Rng, SeedableRng};
// use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::collections::HashMap;
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        // .map(|r| (r.id, r.desc.clone()))
        .collect();
    println!("ID\tScore\tAcc\tCov\tClusterNum");
    for unit in ds.selected_chunks.iter() {
        // (ID,answer)
        let pick_predictions = |ds: &DataSet| -> (Vec<u8>, Vec<u8>) {
            ds.encoded_reads
                .iter()
                .flat_map(|read| {
                    let answer = id2desc[&read.id].starts_with("hapA") as u8;
                    // let answer = id2desc[&read.id].starts_with("chr6_GL000252") as u8;
                    read.nodes
                        .iter()
                        .filter_map(|n| (n.unit == unit.id).then(|| (n.cluster as u8, answer)))
                        .collect::<Vec<_>>()
                })
                .unzip()
        };
        let (answer, before): (Vec<_>, Vec<_>) = pick_predictions(&ds);
        let score = haplotyper::local_clustering::rand_index(&answer, &before);
        let acc = answer
            .iter()
            .zip(before.iter())
            .filter(|(x, y)| x != y)
            .count();
        let acc = acc as f64 / answer.len() as f64;
        let acc = acc.max(1f64 - acc);
        let cl = unit.cluster_num;
        println!("{}\t{}\t{}\t{}\t{}", unit.id, score, acc, before.len(), cl);
    }
    Ok(())
}
