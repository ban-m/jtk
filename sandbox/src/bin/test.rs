use definitions::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::io::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    // let mut ds: DataSet = std::fs::File::open(&args[1])
    //     .map(BufReader::new)
    //     .map(|x| serde_json::de::from_reader(x).unwrap())?;
    // let targets = vec![280, 545];
    // let unit = 1746;
    // let (ids, seqs): (Vec<_>, Vec<_>) = ds
    //     .encoded_reads
    //     .iter()
    //     .filter_map(|r| {
    //         r.nodes
    //             .iter()
    //             .find(|n| n.unit == unit)
    //             .map(|n| (r.id, n.seq()))
    //     })
    //     .unzip();
    // let cov = ds.coverage.unwrap();
    // let mut config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, 2, cov);
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(3424);
    // let (asn, _, _) =
    //     haplotyper::local_clustering::kmeans::clustering(&seqs, &mut rng, &mut config).unwrap();
    // use std::collections::HashMap;
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    // for (id, (readid, asn)) in ids.iter().zip(asn.iter()).enumerate() {
    //     let ans = id2desc[readid].contains("000252v2") as usize;
    //     log::debug!("ANSWER\t{}\t{}\t{}", id, ans, asn);
    // }
    // use haplotyper::em_correction::Context;
    // let cluster_num = 2;
    // let error_rate = 0.13;
    // let num = 20;
    // let span_prob = 0.1;
    // let seed = 2343909;
    // println!("lk\tnull\tspan\tseed\tnum\trand_index\tclusternum");
    // for num in 15..40 {
    //     for span_prob in (0..)
    //         .map(|i| 0.05 + 0.02 * i as f64)
    //         .take_while(|&x| x < 0.5)
    //     {
    //         for seed in 0..20 {
    // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
    // let contexts: Vec<_> = (0..cluster_num)
    //     .flat_map(|i| {
    //         (0..num)
    //             .map(|id| {
    //                 let id = i * num + id;
    //                 let index = 0;
    //                 let unit = 1;
    //                 let cluster = 0; //gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                 let mut forward = vec![(4, 0)];
    //                 if rng.gen_bool(span_prob) {
    //                     let elm = gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                     forward.push((0, elm));
    //                 };
    //                 let mut backward = vec![(3, 0)];
    //                 if rng.gen_bool(span_prob) {
    //                     let elm = gen_cluster(&mut rng, error_rate, i, cluster_num);
    //                     backward.push((2, elm));
    //                 };
    //                 Context::with_attrs(id, index, unit, cluster, forward, backward)
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect();
    // let (null_asn, null_lk, null_ofs) =
    //     haplotyper::em_correction::em_clustering_inner(&contexts, 1, &mut rng);
    // let (asns, lk, ofs) = haplotyper::em_correction::em_clustering_inner(&contexts, 2, &mut rng);
    // let answer: Vec<_> = (0..cluster_num)
    //     .flat_map(|i| vec![i as u8; num as usize])
    //     .collect();
    // let (pred, k) = if null_lk - null_ofs < lk - ofs {
    //     (asns, 2)
    // } else {
    //     (null_asn, 1)
    // };
    // let pred: Vec<_> = pred.iter().map(|x| x.2 as u8).collect();
    // let randindex = haplotyper::local_clustering::rand_index(&answer, &pred);
    // println!("NULL:{},k=2:{}", null_ofs, ofs);
    // println!(
    //     "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
    //     lk, null_lk, span_prob, seed, num, randindex, k
    // );
    //         }
    //     }
    // }
    Ok(())
}

// fn gen_cluster<R: Rng>(rng: &mut R, error_rate: f64, target: u64, max: u64) -> u64 {
//     match rng.gen_bool(error_rate) {
//         true => target,
//         false => (target + rng.gen_range(1..max)) % max,
//     }
// }
