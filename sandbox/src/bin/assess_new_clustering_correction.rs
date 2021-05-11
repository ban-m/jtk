use definitions::*;
use haplotyper::em_correction::*;
use std::io::BufReader;
// use rand::{Rng, SeedableRng};
// use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    // let unit = 0;
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(342);
    // let forward_len = 6;
    // let skip_prob = 0.10;
    // let error_rate = 0.10;
    // let num_reads = 20;
    // let collupsed: std::collections::HashSet<_> = vec![2, 4].into_iter().collect();
    // let reads: Vec<_> = (0..num_reads)
    //     .map(|id| {
    //         let cluster = if rng.gen_bool(error_rate) {
    //             rng.gen_range(0..2)
    //         } else {
    //             2 * id / num_reads
    //         };
    //         let forward: Vec<_> = (0..forward_len)
    //             .filter(|_| rng.gen_bool(1f64 - skip_prob))
    //             .map(|x| {
    //                 if collupsed.contains(&x) {
    //                     (x, 0)
    //                 } else {
    //                     (x, cluster)
    //                 }
    //             })
    //             .collect();
    //         Context::with_attrs(id, 0, unit, cluster, forward, vec![])
    //     })
    //     .collect();
    // let (result, lk) = em_clustering_inner(&reads, 2, &mut rng);
    // println!("{}", lk);
    // for (id, index, cluster) in result {
    //     println!("{}\t{}\t{}", id, index, cluster);
    // }
    use std::collections::HashMap;
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|r| serde_json::de::from_reader(r).unwrap())?;
    let id2desc: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|r| (r.id, r.name.clone()))
        .collect();
    println!("SCORE\tID\tScore\tType");
    for unit in ds.selected_chunks.iter() {
        // (ID,answer)
        let pick_predictions = |ds: &DataSet| -> (Vec<u8>, Vec<u8>) {
            ds.encoded_reads
                .iter()
                .flat_map(|read| {
                    let answer = id2desc[&read.id].starts_with("hapA") as u8;
                    read.nodes
                        .iter()
                        .filter_map(|n| (n.unit == unit.id).then(|| (n.cluster as u8, answer)))
                        .collect::<Vec<_>>()
                })
                .unzip()
        };
        let (answer, before): (Vec<_>, Vec<_>) = pick_predictions(&ds);
        let before_score = haplotyper::local_clustering::rand_index(&answer, &before);
        let acc = answer
            .iter()
            .zip(before.iter())
            .filter(|(x, y)| x != y)
            .count();
        let acc = acc as f64 / answer.len() as f64;
        let acc = acc.max(1f64 - acc);
        println!("SCORE\t{}\t{}\t{}", unit.id, before_score, acc);
    }
    Ok(())
}
