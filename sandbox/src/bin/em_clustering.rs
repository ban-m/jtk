use log::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;

use sandbox::em_clustering;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    debug!("Start");
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(123980);
    use rand::seq::SliceRandom;
    let k = 2;
    let _true_cluster = 4;
    let trials = 8;
    let error = 0.1;
    let num = 20;
    let choises: Vec<_> = (0..k).collect();
    let data: Vec<_> = vec![(0, 0), (0, 1), (1, 0), (1, 1)]
        .iter()
        .map(|&(x, y)| vec![vec![x; trials / 2], vec![y; trials / 2]].concat())
        .flat_map(|template| {
            (0..num)
                .map(|_| {
                    template
                        .iter()
                        .map(|&x| {
                            if rng.gen_bool(error) {
                                *choises.choose(&mut rng).unwrap()
                            } else {
                                x
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .collect();
    // Mock data
    // let k = 8;
    // let true_cluster = 4;
    // let trials = 8;
    // let error = 0.1;
    // let num = 20;
    // let choises: Vec<_> = (0..k).collect();
    // let data: Vec<_> = (0..true_cluster)
    //     .flat_map(|k| {
    //         (0..num)
    //             .map(|_| {
    //                 (0..trials)
    //                     .map(|_| {
    //                         if rng.gen_bool(error) {
    //                             *choises.choose(&mut rng).unwrap()
    //                         } else {
    //                             k
    //                         }
    //                     })
    //                     .collect::<Vec<_>>()
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect();
    let (result, aic) = (k / 2..=k * 2)
        .map(|i| em_clustering(&data, trials, k, i, 123901 + i as u64))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    debug!("{}", -2f64 * aic);
    for ((i, asn), xs) in result.iter().enumerate().zip(data) {
        println!("{}\t{:?}\t{}", i, xs, asn);
    }
}

// type Seq = Vec<u8>;
// use std::fs::File;
// use std::io::{BufReader, Read};
// fn main() {
//     let args: Vec<_> = std::env::args().collect();
//     let mut buf = String::new();
//     let mut rdr = BufReader::new(File::open(&args[1]).unwrap());
//     rdr.read_to_string(&mut buf).unwrap();
//     let reads: Vec<_> = buf
//         .split('>')
//         .filter_map(|record| {
//             let mut record = record.split('\n');
//             let name = record.next().unwrap();
//             let seq: Vec<_> = record.flat_map(|x| x.as_bytes()).copied().collect();
//             if seq.is_empty() {
//                 None
//             } else {
//                 Some((name, seq))
//             }
//         })
//         .collect();
//     let queries: Vec<Seq> = reads
//         .iter()
//         .map(|(_, seq)| {
//             seq.iter()
//                 .map(|&x| match x {
//                     b'a' | b'A' => 0,
//                     b'c' | b'C' => 1,
//                     b'g' | b'G' => 2,
//                     b't' | b'T' => 3,
//                     b'-' => 4,
//                     _ => 4,
//                 })
//                 .collect()
//         })
//         .collect();
//     let len = queries.iter().map(|x| x.len()).max().unwrap();
//     let seed = 120;
//     let clustering = hi_clustering(&queries, len, 5, seed);
//     for (x, (idx, r)) in clustering.iter().zip(reads.iter()) {
//         let line: String = r.iter().map(|&x| x as char).collect();
//         eprintln!("{}\t{}\t{}", x, idx, line);
//     }
// }

// fn hi_clustering(reads: &[Seq], len: usize, base: usize, seed: u64) -> Vec<usize> {
//     eprintln!("LK\tEPARAM\tAIC\tBASE_LK\tLEN\tBASE_AIC");
//     let reads: Vec<_> = reads.iter().enumerate().collect();
//     // let clusters = hi_clustering_inner(&reads, len, base, seed);
//     split_reads(&reads, len, base, seed)
//         .unwrap()
//         .into_iter()
//         .map(|x| if x { 1 } else { 0 })
//         .collect()
//     // let mut slots = vec![0; reads.len()];
//     // for (cl, members) in clusters.iter().enumerate() {
//     //     for &m in members.iter() {
//     //         slots[m] = cl;
//     //     }
//     // }
//     // slots
// }

// #[allow(dead_code)]
// fn hi_clustering_inner(
//     reads: &[(usize, &Seq)],
//     len: usize,
//     base: usize,
//     seed: u64,
// ) -> Vec<Vec<usize>> {
//     match split_reads(&reads, len, base, seed) {
//         Some(res) => {
//             let (mut cluster1, mut cluster2) = (vec![], vec![]);
//             for (read, cl) in reads.into_iter().copied().zip(res) {
//                 if cl {
//                     cluster1.push(read);
//                 } else {
//                     cluster2.push(read);
//                 }
//             }
//             let mut result = vec![];
//             result.extend(hi_clustering_inner(&cluster1, len, base, seed + 1));
//             result.extend(hi_clustering_inner(&cluster2, len, base, seed + 2));
//             result
//         }
//         None => {
//             let cl: Vec<_> = reads.iter().map(|&(i, _)| i).collect();
//             vec![cl]
//         }
//     }
// }

// fn split_reads(
//     reads_with_idx: &[(usize, &Seq)],
//     len: usize,
//     base: usize,
//     seed: u64,
// ) -> Option<Vec<bool>> {
//     let reads: Vec<_> = reads_with_idx.iter().map(|x| x.1).collect();
//     let (asn, lk, ef_param) = (0..100)
//         .map(|s| split_reads_inner(&reads, len, base, seed + s))
//         .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
//         .unwrap();
//     // K=1.
//     // let effective_parameters = (2 * len + 1) as f64;
//     let single_model = {
//         let weights = vec![vec![1.]; reads.len()];
//         EMModel::new(&reads, &weights, 1., 0, len, base, 1.)
//     };
//     let base_lk = reads.iter().map(|r| single_model.lk(r)).sum::<f64>();
//     // Compare.
//     let clustering_aic = lk - ef_param;
//     let base_aic = base_lk - len as f64;
//     eprintln!(
//         "{}\t{}\t{}\t{}\t{}\t{}",
//         lk, ef_param, clustering_aic, base_lk, len, base_aic
//     );
//     // if clustering_aic > base_aic {
//     Some(asn)
//     // } else {
//     //     None
//     // }
// }

// fn split_reads_inner(reads: &[&Seq], len: usize, base: usize, seed: u64) -> (Vec<bool>, f64, f64) {
//     // K=2
//     let mut rng: Xoshiro256PlusPlus = rand::SeedableRng::seed_from_u64(seed);
//     // Random assignment
//     let k = 2;
//     let mut weights: Vec<_> = reads
//         .iter()
//         .map(|_| {
//             let mut weights = vec![0.; k];
//             weights[rng.gen::<usize>() % k] = 1.;
//             weights
//         })
//         .collect();
//     let prior = 0.0000000001;
//     let mut models: Vec<_> = (0..k)
//         .map(|cl| {
//             let frac = weights.iter().map(|ws| ws[cl]).sum::<f64>() / reads.len() as f64;
//             EMModel::new(&reads, &weights, frac, cl, len, base, prior)
//         })
//         .collect();
//     // Log likelihood
//     let mut lk = std::f64::NEG_INFINITY;
//     loop {
//         // Update assignment
//         for (ws, r) in weights.iter_mut().zip(reads.iter()) {
//             let lks: Vec<_> = models.iter().map(|m| m.lk(r)).collect();
//             let sum = logsumexp(&lks);
//             *ws = lks.iter().map(|lk| (lk - sum).exp()).collect();
//             assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
//         }
//         // Update model
//         models = (0..k)
//             .map(|cl| {
//                 let frac = weights.iter().map(|ws| ws[cl]).sum::<f64>() / reads.len() as f64;
//                 EMModel::new(&reads, &weights, frac, cl, len, base, prior)
//             })
//             .collect();
//         // Check log-likelihood.
//         let next_lk = reads
//             .iter()
//             .map(|r| {
//                 let lks: Vec<_> = models.iter().map(|m| m.lk(r)).collect();
//                 logsumexp(&lks)
//             })
//             .sum::<f64>();
//         eprintln!("{}", next_lk);
//         assert!(next_lk >= lk, "{},{}", next_lk, lk);
//         if next_lk - lk < 0.00001 {
//             break;
//         } else {
//             lk = next_lk;
//         }
//     }
//     let effective_parameters = models[0]
//         .parameters
//         .iter()
//         .zip(models[1].parameters.iter())
//         .map(|(p1, p2)| {
//             // for (x1, x2) in p1.iter().zip(p2.iter()) {
//             //     eprint!("{:.3}:{:.3}\t", x1, x2);
//             // }
//             // eprintln!();
//             p1.iter()
//                 .zip(p2)
//                 .map(|(x1, x2)| (x1 - x2).powi(2))
//                 .sum::<f64>()
//         })
//         .sum::<f64>()
//         + 1.
//         + len as f64;
//     let asn: Vec<_> = weights.iter().map(|xs| xs[0] > xs[1]).collect();
//     eprintln!("{}\t{:?}", seed, lk);
//     (asn, lk, effective_parameters)
// }

// struct EMModel {
//     // _base: usize,
//     // len: usize,
//     frac: f64,
//     // Position -> Base
//     parameters: Vec<Vec<f64>>,
// }

// impl std::fmt::Display for EMModel {
//     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//         for p in self.parameters.iter() {
//             let line: Vec<_> = p.iter().map(|x| format!("{:.3}", x)).collect();
//             writeln!(f, "{}", line.join(","))?;
//         }
//         Ok(())
//     }
// }

// impl EMModel {
//     fn new(
//         reads: &[&Seq],
//         weights: &[Vec<f64>],
//         frac: f64,
//         cl: usize,
//         len: usize,
//         base: usize,
//         pcount: f64,
//     ) -> Self {
//         let mut parameters = vec![vec![pcount; base]; len];
//         for (read, ws) in reads.iter().zip(weights.iter()) {
//             let w = ws[cl];
//             for (idx, &b) in read.iter().enumerate() {
//                 parameters[idx][b as usize] += w;
//             }
//         }
//         for ws in parameters.iter_mut() {
//             let sum = ws.iter().sum::<f64>();
//             ws.iter_mut().for_each(|w| *w /= sum);
//         }
//         Self {
//             // base,
//             // len,
//             frac,
//             parameters,
//         }
//     }
//     // Log likelihood
//     fn lk(&self, read: &Seq) -> f64 {
//         let x1 = read
//             .iter()
//             .enumerate()
//             .map(|(idx, &b)| self.parameters[idx][b as usize].ln())
//             .sum::<f64>();
//         self.frac.ln() + x1
//     }
// }

// fn logsumexp(xs: &[f64]) -> f64 {
//     if xs.is_empty() {
//         return 0.;
//     }
//     let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
//     let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
//     assert!(sum >= 0., "{:?}->{}", xs, sum);
//     max + sum
// }
