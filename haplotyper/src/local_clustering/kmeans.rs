//! A small K-means clustering algorithm.
use rand::Rng;
const DEFAULT_GAIN: f64 = 1f64;
#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    band_width: usize,
    // Coverage for haploid.
    coverage: f64,
    pub cluster_num: u8,
}

impl ClusteringConfig {
    pub fn new(band_width: usize, cluster_num: u8, coverage: f64) -> Self {
        Self {
            band_width,
            coverage,
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
        coverage,
    } = config.clone();
    let cons_template = kiley::consensus(reads, rng.gen(), 10, band_width);
    if cluster_num == 0 {
        return Some((vec![0; reads.len()], cons_template));
    }
    let cluster_choices =
        (cluster_num.max(3) - 2..=cluster_num + 2).flat_map(|k| std::iter::repeat(k).take(10));
    let profiles = get_profiles(&cons_template, reads, band_width as isize);
    // Filter profiles by doubling.
    let var_num = 2 * (cluster_num - 1).max(1) as usize;
    let profiles: Vec<Vec<_>> = {
        let probes = filter_profiles(&profiles, cons_template.len(), cluster_num, var_num);
        let profiles: Vec<Vec<_>> = profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect();
        let (assignments, _score, cluster_num) = std::iter::repeat(cluster_num)
            .take(10)
            // cluster_choices
            // .clone()
            .map(|k| {
                let (asn, score) = mcmc_clustering(&profiles, k, coverage, rng);
                let mod_score = score - DEFAULT_GAIN * coverage * (k - 1) as f64;
                (asn, mod_score, k)
            })
            .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
            .unwrap();
        // TODO: How to determine the upper bound of the number of variants?
        let to_used = retrieve_used_positions(&assignments, &profiles, cluster_num as usize);
        if to_used.iter().all(|&x| !x) {
            return Some((vec![0; reads.len()], cons_template));
        }
        //  else if to_used.iter().all(|&x| x) && var_num < 20 {
        //     var_num *= 2;
        // } else {
        //     // Filter out unused variants. 2nd clustering.
        //     break profiles
        //         .iter()
        //         .map(|xs| {
        //             xs.iter()
        //                 .zip(to_used.iter())
        //                 .filter_map(|(&x, &y)| y.then(|| x))
        //                 .collect()
        //         })
        //         .collect();
        // }
        profiles
            .iter()
            .map(|xs| {
                xs.iter()
                    .zip(to_used.iter())
                    .filter_map(|(&x, &y)| y.then(|| x))
                    .collect()
            })
            .collect()
    };
    let (assignments, _score, cluster_num) = cluster_choices
        .map(|k| {
            let (asn, score) = mcmc_clustering(&profiles, k, coverage, rng);
            let mod_score = score - DEFAULT_GAIN * coverage * (k - 1) as f64;
            (asn, mod_score, k)
        })
        .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
        .unwrap();
    // let (assignments, _score) = std::iter::repeat(cluster_num)
    //     .take(10)
    //     .map(|k| mcmc_clustering(&profiles, k, coverage, rng))
    //     .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
    //     .unwrap();
    // for (asn, prf) in assignments.iter().zip(profiles) {
    //     let prf: Vec<_> = prf.iter().map(|x| format!("{:.1}", x)).collect();
    //     debug!("DUMP\t{}\t{}", asn, prf.join("\t"));
    // }
    config.cluster_num = cluster_num;
    Some((assignments, cons_template))
}

// i->true if there's a cluster using the i-th position to improve the total likelihood.
fn retrieve_used_positions(
    assignments: &[u8],
    profiles: &[Vec<f64>],
    cluster_num: usize,
) -> Vec<bool> {
    let len = profiles[0].len();
    let mut lks = vec![vec![0f64; len]; cluster_num];
    for (&asn, xs) in assignments.iter().zip(profiles.iter()) {
        for (l, &x) in lks[asn as usize].iter_mut().zip(xs) {
            *l += x;
        }
    }
    lks.iter().fold(vec![false; len], |mut acc, xs| {
        for (is_used, total_lk) in acc.iter_mut().zip(xs) {
            *is_used |= total_lk.is_sign_positive();
        }
        acc
    })
}

fn get_profiles<T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    band_width: isize,
) -> Vec<Vec<f64>> {
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let template = kiley::padseq::PadSeq::new(template);
    let reads: Vec<_> = reads
        .iter()
        .map(|r| kiley::padseq::PadSeq::new(r.borrow()))
        .collect();
    let (hmm, _) = hmm.fit_banded_inner(&template, &reads, band_width as usize);
    reads
        .iter()
        .map(|read| {
            let prof = banded::ProfileBanded::new(&hmm, &template, &read, band_width).unwrap();
            let lk = prof.lk();
            prof.to_modification_table()
                .iter()
                .take(9 * template.len())
                .map(|x| x - lk)
                .collect()
        })
        .collect()
}

fn filter_profiles(
    profiles: &[Vec<f64>],
    template_len: usize,
    cluster_num: u8,
    var_num: usize,
) -> Vec<(usize, f64)> {
    let total_improvement =
        profiles
            .iter()
            .fold(vec![(0f64, 0); template_len * 9], |mut acc, prof| {
                for ((sum, count), p) in acc.iter_mut().zip(prof) {
                    *sum += p.max(0f64);
                    *count += (0.01 < *p) as usize;
                }
                acc
            });
    let mut probes: Vec<(usize, f64)> = total_improvement
        .into_iter()
        .enumerate()
        .map(|(pos, (sum, count))| {
            let max_lk = (1..cluster_num + 1)
                .map(|k| poisson_lk(count, profiles.len() as f64 / k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap();
            (pos, sum + max_lk)
        })
        .collect();
    probes.sort_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap());
    probes.reverse();
    let abs = |x: usize, y: usize| x.saturating_sub(y) + y.saturating_sub(x);
    let mut buffer = vec![];
    for &(pos, lk) in probes.iter() {
        if buffer.iter().all(|(p, _)| abs(p / 9, pos / 9) > 5) {
            buffer.push((pos, lk));
        }
        if buffer.len() == var_num {
            break;
        }
    }
    buffer
}

// Return log Poiss{x|lambda}
fn poisson_lk(x: usize, lambda: f64) -> f64 {
    x as f64 * lambda.ln() - lambda - (1..x + 1).map(|c| (c as f64).ln()).sum::<f64>()
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
#[allow(dead_code)]
fn score(data: &[Vec<f64>], asn: &[u8], k: u8) -> f64 {
    let dim = data[0].len();
    let mut sums = vec![vec![0f64; dim]; k as usize];
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
    // clustering_features(data, &mut assignments, k);
    assignments
}

// Take dataset, the number of the cluster, the haploid coverage and seed generator,
// return the total likelihood.
fn mcmc_clustering<R: Rng>(data: &[Vec<f64>], k: u8, cov: f64, rng: &mut R) -> (Vec<u8>, f64) {
    if k == 1 || data.iter().all(|xs| xs.is_empty()) {
        return (vec![0; data.len()], 0f64);
    }
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
    let total_lk = mcmc(data, &mut assignments, k, cov, rng);
    (assignments, total_lk)
}

fn mcmc<R: Rng>(data: &[Vec<f64>], assign: &mut [u8], k: u8, cov: f64, rng: &mut R) -> f64 {
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
    let mut lk_sums: Vec<f64> = lks
        .iter()
        .map(|xs| xs.iter().map(|x| x.max(0f64)).sum())
        .collect();
    // Likelihood difference at a column.
    let diff_new = |(x, y): (&f64, &f64)| (y + x).max(0f64) - y.max(0f64);
    let diff_old = |(x, y): (&f64, &f64)| (y - x).max(0f64) - y.max(0f64);
    let total = 50_000 * k as usize;
    // MAP estimation.
    let (mut max, mut argmax) = (std::f64::NEG_INFINITY, vec![]);
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
            lk_sums[old] = lks[old].iter().map(|x| x.max(0f64)).sum();
            for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
                *lk += x;
            }
            lk_sums[new] = lks[new].iter().map(|x| x.max(0f64)).sum();
        }
        let lk =
            lk_sums.iter().sum::<f64>() + clusters.iter().map(|&x| poisson_lk(x, cov)).sum::<f64>();
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
