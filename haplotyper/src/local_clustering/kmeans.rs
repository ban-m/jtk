//! A small K-means clustering algorithm.
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    band_width: usize,
    pub cluster_num: u8,
}

impl ClusteringConfig {
    pub fn new(band_width: usize, cluster_num: u8) -> Self {
        Self {
            band_width,
            cluster_num,
        }
    }
}

/// Return the assignments and the consensus sequence.
/// The number of the cluster would be modified as the number of optimal clustering.
pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &mut ClusteringConfig,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let ClusteringConfig {
        band_width,
        cluster_num,
    } = config.clone();
    let cons_template = kiley::consensus(reads, rng.gen(), 10, band_width);
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let template = kiley::padseq::PadSeq::new(cons_template.as_slice());
    let reads: Vec<_> = reads
        .iter()
        .map(|r| kiley::padseq::PadSeq::new(r.borrow()))
        .collect();
    let (hmm, _) = hmm.fit_banded_inner(&template, &reads, band_width);
    let profiles: Vec<Vec<_>> = reads
        .iter()
        .map(|read| {
            let prof =
                banded::ProfileBanded::new(&hmm, &template, &read, band_width as isize).unwrap();
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
            let (sum, count): (f64, usize) = profiles.iter().fold((0f64, 0), |(sum, c), xs| {
                if 0.01 < xs[pos] {
                    (sum + xs[pos], c + 1)
                } else {
                    (sum, c)
                }
            });
            let max_lk = (1..cluster_num)
                .map(|k| poisson_lk(count, reads.len() as f64 / k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap();
            (pos, sum + max_lk)
        })
        .collect();
    let mut probes: Vec<_> = probes
        .chunks_exact(9)
        .filter_map(|xs| xs.iter().max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap()))
        .copied()
        .collect();
    probes.sort_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap());
    probes.reverse();
    probes = {
        let abs = |x: usize, y: usize| x.saturating_sub(y) + y.saturating_sub(x);
        let (mut buffer, mut pointer) = (vec![], 0);
        while buffer.len() < var_num {
            let (pos, _) = probes[pointer];
            if buffer.iter().all(|(p, _)| abs(p / 9, pos / 9) > 5) {
                buffer.push(probes[pointer]);
            }
            pointer += 1;
            if probes.len() <= pointer {
                break;
            }
        }
        buffer
    };
    debug!("{:?}", probes);
    let total_lk = {
        let lk: Vec<_> = probes.iter().map(|x| x.1).collect();
        logsumexp(&lk)
    };
    // let norm: f64 = probes.iter().map(|(_, s)| s * s).sum();
    // let norm = norm.sqrt();
    // probes.iter_mut().for_each(|x| x.1 /= norm);
    probes.iter_mut().for_each(|x| x.1 = (x.1 - total_lk).exp());
    let profiles: Vec<Vec<_>> = profiles
        .iter()
        .map(|xs| {
            probes
                .iter()
                .map(|&(pos, scale)| (xs[pos] * scale))
                .collect()
        })
        .collect();
    // let max_gain: f64 = profiles
    //     .iter()
    //     .map(|xs| -> f64 { xs.iter().map(|x| x.max(0f64)).sum() })
    //     .sum();
    // debug!("Per cluster:{:.1}", max_gain / var_num as f64);
    // let no_cluster = (vec![0; reads.len()], -max_gain / var_num as f64, 1);
    let no_cluster = (vec![0; reads.len()], std::f64::NEG_INFINITY, 1);
    let (assignments, score, cluster_num) = std::iter::once(cluster_num)
        //(cluster_num.max(4) - 2..(cluster_num + 2))
        .flat_map(|k| std::iter::repeat(k).take(1))
        .map(|cluster_num| {
            let (asn, score) = mcmc_clustering(&profiles, cluster_num, rng);
            (asn, score, cluster_num)
        })
        .fold(no_cluster, |(argmax, max, cl), (asn, score, k)| {
            if max < score {
                (asn, score, k)
            } else {
                (argmax, max, cl)
            }
        });
    debug!("SCORE:{}", score);
    let max: Vec<_> = profiles
        .iter()
        .fold(vec![0f64; profiles[0].len()], |mut acc, xs| {
            acc.iter_mut().zip(xs).for_each(|(a, x)| *a += x.max(0f64));
            acc
        });
    let max: Vec<_> = max.iter().map(|x| format!("{:.2}", x)).collect();
    debug!("MAX:{}", max.join("\t"));
    for (asn, prf) in assignments.iter().zip(profiles) {
        let prf: Vec<_> = prf.iter().map(|x| format!("{:.1}", x)).collect();
        debug!("DUMP\t{}\t{}", asn, prf.join("\t"));
    }
    config.cluster_num = cluster_num;
    Some((assignments, cons_template))
}

fn poisson_lk(count: usize, mean: f64) -> f64 {
    count as f64 * mean.ln() - mean - (1..count + 1).map(|x| (x as f64).ln()).sum::<f64>()
}

// Likelihood of each size of the clsuters.
#[allow(dead_code)]
fn poisson_likelihood(asn: &[u8], cluster_num: u8) -> f64 {
    let mut count = vec![0; cluster_num as usize];
    for &a in asn.iter() {
        count[a as usize] += 1;
    }
    let mean = asn.len() as f64 / cluster_num as f64;
    count.into_iter().map(|c| poisson_lk(c, mean)).sum()
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
        .map(|xs| -> f64 { xs.iter().map(|&x| x.max(0f64)).sum() })
        .sum()
}

// #[allow(dead_code)]
// fn clustering_features(data: &[Vec<f64>], assignments: &mut [u8], k: u8) {
//     let mut is_updated = true;
//     let mut models: Vec<Vec<_>> = {
//         let mut sum = vec![vec![0f64; data[0].len()]; k as usize];
//         for (&asn, d) in assignments.iter().zip(data.iter()) {
//             for (s, x) in sum[asn as usize].iter_mut().zip(d) {
//                 *s += x;
//             }
//         }
//         sum.iter()
//             .map(|xs| xs.iter().map(|x| x.is_sign_positive()).collect())
//             .collect()
//     };
//     let mut weights: Vec<Vec<_>> = assignments
//         .iter()
//         .map(|&asn| {
//             let mut weight = vec![0f64; k as usize];
//             weight[asn as usize] = 1f64;
//             weight
//         })
//         .collect();
//     while is_updated {
//         is_updated = false;
//         // Update weights.
//         for (weight, xs) in weights.iter_mut().zip(data.iter()) {
//             for (w, model) in weight.iter_mut().zip(models.iter()) {
//                 *w = model
//                     .iter()
//                     .zip(xs)
//                     .filter_map(|(&m, x)| m.then(|| x))
//                     .sum();
//             }
//             let tot = logsumexp(&weight);
//             weight.iter_mut().for_each(|x| *x = (*x - tot).exp());
//         }
//         for (cl, model) in models.iter_mut().enumerate() {
//             let mut sums: Vec<_> = vec![0f64; data[0].len()];
//             for (weight, xs) in weights.iter().zip(data.iter()) {
//                 for (s, x) in sums.iter_mut().zip(xs.iter()) {
//                     *s += x * weight[cl];
//                 }
//             }
//             for (m, s) in model.iter_mut().zip(sums.iter()) {
//                 let new_model = s.is_sign_positive();
//                 is_updated |= new_model != *m;
//                 *m = new_model;
//             }
//         }
//     }
//     for (asn, weight) in assignments.iter_mut().zip(weights) {
//         let (new_asn, _) = weight
//             .iter()
//             .enumerate()
//             .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
//             .unwrap();
//         *asn = new_asn as u8;
//     }
// }

#[allow(dead_code)]
fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

#[allow(dead_code)]
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
    // clustering_features(data, &mut assignments, k);
    assignments
}

fn mcmc_clustering<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> (Vec<u8>, f64) {
    // 1. Construct the first assignments.
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
    let total_lk = mcmc(data, &mut assignments, k, rng);
    (assignments, total_lk)
}

fn mcmc<R: Rng>(data: &[Vec<f64>], assign: &mut [u8], k: u8, rng: &mut R) -> f64 {
    // how many instance in a cluster.
    let mut clusters = vec![0; k as usize];
    for &asn in assign.iter() {
        clusters[asn as usize] += 1;
    }
    // Current (un-modified) likelihoods.
    let mut lks = vec![vec![0f64; data[0].len()]; k as usize];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        for (lk, x) in lks[asn as usize].iter_mut().zip(xs) {
            *lk += x;
        }
    }
    // Likelihood difference at a column.
    let diff_new = |(x, y): (&f64, &f64)| (y + x).max(0f64) - y.max(0f64);
    let diff_old = |(x, y): (&f64, &f64)| (y - x).max(0f64) - y.max(0f64);
    let total_lk = |lks: &[Vec<f64>]| -> f64 {
        lks.iter()
            .flat_map(|lk| lk.iter())
            .map(|x| x.max(0f64))
            .sum()
    };
    // TODO: Multiplicity estimation.
    // Burn in.
    let burn_in = 10_000;
    for _ in 0..burn_in {
        let (idx, new) = (rng.gen_range(0..data.len()), rng.gen_range(0..k));
        let (old, xs) = (assign[idx], &data[idx]);
        let (old, new) = (old as usize, new as usize);
        let pois_diff = (clusters[old] as f64 / (clusters[new] as f64 + 1f64)).ln();
        let new_diff: f64 = xs.iter().zip(lks[new].iter()).map(diff_new).sum();
        let old_diff: f64 = xs.iter().zip(lks[old].iter()).map(diff_old).sum();
        let total_diff = (pois_diff + old_diff + new_diff).exp().min(1f64);
        if (old != new) && rng.gen_bool(total_diff) {
            assign[idx] = new as u8;
            clusters[old] -= 1;
            clusters[new] += 1;
            for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
                *lk -= x;
            }
            for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
                *lk += x;
            }
        }
    }
    let total = 100_000;
    // MAP estimation.
    let (mut max, mut argmax) = (std::f64::NEG_INFINITY, vec![]);
    let mean = data.len() as f64 / k as f64;
    for _t in 0..total {
        let (idx, new) = (rng.gen_range(0..data.len()), rng.gen_range(0..k));
        let (old, xs) = (assign[idx], &data[idx]);
        if old == new {
            continue;
        }
        let (old, new) = (old as usize, new as usize);
        let pois_diff = (clusters[old] as f64 / (clusters[new] as f64 + 1f64)).ln();
        let old_diff: f64 = xs.iter().zip(lks[old].iter()).map(diff_old).sum();
        let new_diff: f64 = xs.iter().zip(lks[new].iter()).map(diff_new).sum();
        let total_diff = pois_diff + old_diff + new_diff;
        if rng.gen_bool(total_diff.exp().min(1f64)) {
            assign[idx] = new as u8;
            clusters[old] -= 1;
            clusters[new] += 1;
            for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
                *lk -= x;
            }
            for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
                *lk += x;
            }
        }
        let lk = total_lk(&lks) + clusters.iter().map(|&x| poisson_lk(x, mean)).sum::<f64>();
        // debug!("MCMC\t{}\t{}\t{}\tMCMC", t, lk, total_diff);
        if max < lk {
            max = lk;
            argmax = assign.to_vec();
        }
    }
    assign.iter_mut().zip(argmax).for_each(|(x, y)| *x = y);
    max
}

// fn kmeans_f64<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
//     let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
//     kmeans_f64_with_init(data, &mut assignments, k);
//     assignments
// }

// Return the distance between xs and ys.
fn euclid_norm_f64(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(ys.len(), xs.len());
    xs.iter()
        .zip(ys.iter())
        .map(|(x, &y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}
