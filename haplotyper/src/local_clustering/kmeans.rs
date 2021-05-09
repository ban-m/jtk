//! A small K-means clustering algorithm.
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    band_width: usize,
    probe_num: usize,
    retry_num: usize,
    cluster_num: u8,
    subbatch_num: usize,
}

impl ClusteringConfig {
    pub fn new(
        band_width: usize,
        probe_num: usize,
        retry_num: usize,
        cluster_num: u8,
        subbatch_num: usize,
    ) -> Self {
        Self {
            band_width,
            probe_num,
            retry_num,
            cluster_num,
            subbatch_num,
        }
    }
}

pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<Vec<u8>> {
    let ClusteringConfig {
        band_width,
        cluster_num,
        ..
    } = *config;
    let cons_template = kiley::consensus(reads, rng.gen(), 10, band_width);
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let template = kiley::padseq::PadSeq::new(cons_template.as_slice());
    let reads: Vec<_> = reads
        .iter()
        .map(|r| kiley::padseq::PadSeq::new(r.borrow()))
        .collect();
    let (hmm, _) = hmm.fit_banded_inner(&template, &reads, 100);
    let profiles: Vec<Vec<_>> = reads
        .iter()
        .map(|read| {
            let prof = banded::ProfileBanded::new(&hmm, &template, &read, 100).unwrap();
            let lk = prof.lk();
            prof.to_modification_table()
                .iter()
                .take(9 * template.len())
                .map(|x| x - lk)
                .collect()
        })
        .collect();
    let var_num = 2 * cluster_num as usize;
    let probes: Vec<_> = (0..9 * template.len())
        .map(|pos| {
            let mut xs: Vec<_> = profiles.iter().map(|xs| xs[pos]).collect();
            xs.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let sum: f64 = xs.iter().filter(|&&x| 0f64 < x).sum();
            let median = xs[xs.len() / 2];
            (pos, sum, median)
        })
        .collect();
    let mut probes: Vec<_> = probes
        .chunks_exact(9)
        .filter_map(|xs| xs.iter().max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap()))
        .copied()
        .collect();
    probes.sort_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap());
    probes.reverse();
    probes.truncate(var_num);
    let norm: f64 = probes.iter().map(|(_, s, _)| s * s).sum();
    let norm = norm.sqrt();
    let profiles: Vec<Vec<_>> = profiles
        .iter()
        .map(|xs| {
            probes
                .iter()
                .map(|&(pos, scale, _)| (xs[pos] * scale / norm))
                .collect()
        })
        .collect();
    // We can tune here!
    let (assignments, _score) = (0..10)
        .map(|_| {
            let asn = kmeans_f64_plusplus(&profiles, cluster_num, rng);
            let score = score(&profiles, &asn, cluster_num as usize);
            (asn, score)
        })
        .fold(
            (vec![], std::f64::NEG_INFINITY),
            |(argmax, max), (asn, score)| {
                if max < score {
                    (asn, score)
                } else {
                    (argmax, max)
                }
            },
        );
    Some(assignments)
}

// Return the gain of likelihood for a given dataset.
// To calculate the gain,we compute the following metrics:
// 1. Sum up the vectors for each cluster.
// 2. For each sum, the element of the positve values would be selected,
// or we acceept the edit operation at that point.
// 3. Sum up the positive values of summing-upped vector for each cluster.
fn score(data: &[Vec<f64>], asn: &[u8], k: usize) -> f64 {
    let dim = data[0].len();
    let mut sums = vec![vec![0f64; dim]; k];
    for (xs, &asn) in data.iter().zip(asn.iter()) {
        xs.iter()
            .zip(sums[asn as usize].iter_mut())
            .for_each(|(x, y)| *y += x);
    }
    sums.iter()
        .map(|xs| -> f64 { xs.iter().filter(|&&x| 0f64 < x).sum() })
        .sum()
}

// fn get_profiles<R: Rng>(
//     template: &[u8],
//     asn: &[u8],
//     reads: &[&[u8]],
//     rng: &mut R,
//     config: &ClusteringConfig,
// ) -> Option<Vec<Vec<i32>>> {
//     let probes: Vec<_> = (0..config.probe_num)
//         .map(|_| {
//             let picked = pick_reads(asn, reads, config, rng);
//             polish_until_converge_banded(&template, &picked, config.band_width)
//         })
//         .collect();
//     let mismatch = |xs: &[u8], ys: &[u8]| {
//         edlib_sys::global(xs, ys)
//             .iter()
//             .filter(|&&x| x == 3)
//             .count()
//     };
//     let profiles: Vec<Vec<_>> = reads
//         .iter()
//         .map(|r| {
//             let center = mismatch(r, &template) as i32;
//             probes
//                 .iter()
//                 .map(|p| mismatch(p, r) as i32 - center)
//                 .collect()
//         })
//         .collect();
//     let selected_position: Vec<bool> = {
//         // Transpose profiles.
//         let mut t_profiles: Vec<_> = vec![vec![]; probes.len()];
//         for prf in profiles.iter() {
//             for (i, x) in prf.iter().enumerate() {
//                 t_profiles[i].push(x);
//             }
//         }
//         t_profiles
//             .iter()
//             .map(|prf| prf.iter().any(|&x| x != prf[0]))
//             .collect()
//     };
//     selected_position.iter().any(|&x| x).then(|| {
//         profiles
//             .iter()
//             .map(|xs| {
//                 xs.iter()
//                     .zip(selected_position.iter())
//                     .filter_map(|(&x, &b)| b.then(|| x))
//                     .collect()
//             })
//             .collect()
//     })
// }

// fn pick_reads<'a, R: Rng>(
//     asn: &[u8],
//     reads: &[&'a [u8]],
//     config: &ClusteringConfig,
//     rng: &mut R,
// ) -> Vec<&'a [u8]> {
//     use rand::seq::SliceRandom;
//     let max = asn.iter().max().copied().unwrap() as usize + 1;
//     let mut bucket = vec![vec![]; max];
//     for (&asn, read) in asn.iter().zip(reads.iter()) {
//         bucket[asn as usize].push(read);
//     }
//     bucket.iter_mut().for_each(|x| x.shuffle(rng));
//     let mut counts = vec![3; max];
//     for (i, b) in bucket.iter().enumerate() {
//         if b.is_empty() {
//             counts[i] = 0;
//         }
//     }
//     let choises: Vec<_> = (0..max).collect();
//     (0..config.subbatch_num)
//         .filter_map(|_| {
//             let denom: usize = counts.iter().sum();
//             if denom == 0 {
//                 return None;
//             }
//             let chosen = *choises
//                 .choose_weighted(rng, |&a| counts[a] as f64 / denom as f64)
//                 .ok()?;
//             let picked = bucket[chosen].pop().unwrap();
//             counts[chosen] += 1;
//             if bucket[chosen].is_empty() {
//                 counts[chosen] = 0;
//             }
//             Some(*picked)
//         })
//         .collect()
// }

// fn try_clustering<R: Rng>(
//     reads: &[&[u8]],
//     rng: &mut R,
//     config: &ClusteringConfig,
// ) -> Option<Vec<u8>> {
//     let ClusteringConfig {
//         band_width,
//         cluster_num,
//         subbatch_num,
//         ..
//     } = *config;
//     let cons_template = kiley::consensus(&reads, rng.gen(), 10, band_width);
//     debug!("Polished,{}", reads.len());
//     let probes: Vec<_> = (0..4)
//         .map(|_| {
//             use kiley::*;
//             use rand::seq::SliceRandom;
//             let picked: Vec<_> = reads.choose_multiple(rng, subbatch_num).copied().collect();
//             let config = PolishConfig::new(100, 0, picked.len(), 0, 0);
//             polish_chunk_by_parts(&cons_template, &picked, &config)
//         })
//         .collect();
//     use kiley::gphmm::*;
//     let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
//     let center: Vec<_> = reads
//         .iter()
//         .map(|read| hmm.likelihood_banded(&cons_template, read, 100).unwrap())
//         .collect();
//     let profiles: Vec<Vec<_>> = reads
//         .iter()
//         .zip(center.iter())
//         .map(|(read, c)| {
//             probes
//                 .iter()
//                 .map(|p| hmm.likelihood_banded(p, read, 100).unwrap() - c)
//                 .collect()
//         })
//         .collect();
//     for prof in profiles.iter() {
//         let prof: Vec<_> = prof.iter().map(|x| format!("{:.0}", x)).collect();
//         debug!("{}", prof.join("\t"));
//     }
//     Some(kmeans_f64(&profiles, cluster_num, rng))
// }

// TODO:Maybe we need more sophisticated clustering algorithm here.
// DBSCAN?
// Kmeans with random seeds?
// Or, maybe we can tune the number of clustering?
// fn kmeans<R: Rng>(data: &[Vec<i32>], k: u8, rng: &mut R) -> Vec<u8> {
//     use rand::seq::SliceRandom;
//     let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
//     loop {
//         let centers: Vec<Vec<f64>> = (0..k)
//             .filter_map(|cl| {
//                 let (mut count, mut slots) = (0, vec![0; data[0].len()]);
//                 let filtered = data
//                     .iter()
//                     .zip(assignments.iter())
//                     .filter_map(|(d, &a)| (a == cl).then(|| d));
//                 for datum in filtered {
//                     assert_eq!(slots.len(), datum.len());
//                     slots.iter_mut().zip(datum).for_each(|(acc, x)| *acc += x);
//                     count += 1;
//                 }
//                 let center: Vec<_> = slots.iter().map(|&x| x as f64 / count as f64).collect();
//                 (count != 0).then(|| center)
//             })
//             .collect();
//         let new_assignments: Vec<_> = data
//             .iter()
//             .filter_map(|x| {
//                 let mut dists: Vec<_> = centers
//                     .iter()
//                     .enumerate()
//                     .map(|(i, center)| (i as u8, euclid_norm(center, x)))
//                     .collect();
//                 let min = dists
//                     .iter()
//                     .map(|x| x.1)
//                     .min_by(|x, y| x.partial_cmp(y).unwrap())
//                     .unwrap();
//                 dists.retain(|x| x.1 - min < 0.000001);
//                 assert!(!dists.is_empty());
//                 dists.choose(rng).map(|x| x.0)
//             })
//             .collect();
//         assert_eq!(new_assignments.len(), assignments.len());
//         if new_assignments == assignments {
//             break assignments;
//         } else {
//             assignments = new_assignments;
//         }
//     }
// }
// fn euclid_norm(xs: &[f64], ys: &[i32]) -> f64 {
//     assert_eq!(ys.len(), xs.len());
//     xs.iter()
//         .zip(ys.iter())
//         .map(|(x, &y)| (x - y as f64).powi(2))
//         .sum()
// }

fn kmeans_f64_with_init(data: &[Vec<f64>], assignments: &mut [u8], k: u8) {
    let mut is_updated = true;
    while is_updated {
        let centers: Vec<Vec<f64>> = (0..k)
            .filter_map(|cl| {
                let (mut count, mut slots) = (0, vec![0f64; data[0].len()]);
                let filtered = data
                    .iter()
                    .zip(assignments.iter())
                    .filter_map(|(d, &a)| (a == cl).then(|| d));
                for datum in filtered {
                    assert_eq!(slots.len(), datum.len());
                    slots.iter_mut().zip(datum).for_each(|(acc, x)| *acc += x);
                    count += 1;
                }
                let center: Vec<_> = slots.iter().map(|&x| x as f64 / count as f64).collect();
                (count != 0).then(|| center)
            })
            .collect();
        is_updated = false;
        for (x, asn) in data.iter().zip(assignments.iter_mut()) {
            let (new_asn, _) = centers
                .iter()
                .enumerate()
                .map(|(i, center)| (i as u8, euclid_norm_f64(center, x)))
                .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
                .unwrap();
            if new_asn != *asn {
                is_updated = true;
                *asn = new_asn;
            }
        }
    }
}

#[allow(dead_code)]
fn kmeans_f64_plusplus<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
    let mut centers: Vec<&[f64]> = vec![];
    let indices: Vec<_> = (0..data.len()).collect();
    // Choosing centers.
    use rand::seq::SliceRandom;
    while centers.len() < k as usize {
        // calculate distance to the most nearest centers.
        let mut dists: Vec<_> = data
            .iter()
            .map(|xs| {
                centers
                    .iter()
                    .map(|c| euclid_norm_f64(xs, c))
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap_or(1f64)
                    .powi(2)
            })
            .collect();
        let total: f64 = dists.iter().sum();
        dists.iter_mut().for_each(|x| *x /= total);
        let idx = *indices.choose_weighted(rng, |&idx| dists[idx]).unwrap();
        centers.push(&data[idx]);
    }
    let mut assignments: Vec<_> = data
        .iter()
        .map(|xs| {
            centers
                .iter()
                .enumerate()
                .map(|(i, center)| (i as u8, euclid_norm_f64(center, xs)))
                .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
                .unwrap()
                .0
        })
        .collect();
    kmeans_f64_with_init(data, &mut assignments, k);
    assignments
}

#[allow(dead_code)]
fn kmeans_f64<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
    kmeans_f64_with_init(data, &mut assignments, k);
    assignments
}

// Return the distance between xs and ys.
fn euclid_norm_f64(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(ys.len(), xs.len());
    xs.iter()
        .zip(ys.iter())
        .map(|(x, &y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}
