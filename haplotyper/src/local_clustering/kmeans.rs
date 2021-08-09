//! A small K-means clustering algorithm.
use rand::Rng;
#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    // TODO: Remove this?
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

/// Usual "flat(non-recursive)" clustering.
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
    let profiles = get_profiles(&cons_template, reads, band_width as isize);
    // for (i, prof) in profiles.iter().enumerate() {
    //     for (j, x) in prof.iter().enumerate() {
    //         debug!("DUMP\t{}\t{}\t{}", i, j, x);
    //     }
    // }
    let cluster_num = cluster_num as usize;
    let selected_variants: Vec<Vec<_>> = {
        let probes = filter_profiles(&profiles, cluster_num, 3, coverage);
        // for (pos, lk) in probes.iter() {
        //     let (idx, t) = (pos / 9, pos % 9);
        //     debug!("POS\t{}\t{}\t{:.3}", idx, t, lk);
        // }
        profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect()
    };
    // let selected_variants = calibrate_likelihoods(&selected_variants, coverage);
    let k_range = cluster_num.max(3) - 2..=cluster_num;
    let (assignments, _score, _cluster_num) = k_range
        .filter_map(|k| {
            let num = 10;
            let (asn, score) = (0..num)
                .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
                .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
            let penalty = (k - 1) as f64 * coverage;
            let copy_num_lk = poisson_lk(reads.len(), coverage * k as f64);
            let mod_score = copy_num_lk + score - penalty;
            // debug!("{}\t{:.2}\t{:.2}", k, score, mod_score);
            Some((asn, mod_score, k))
        })
        .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())?;
    // config.cluster_num = cluster_num as u8;
    // for (i, prf) in assignments.iter().zip(selected_variants.iter()) {
    //     let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
    //     debug!("ASN\t{}\t{}\t{}", cluster_num, i, prf.join("\t"));
    // }
    Some((assignments, cons_template))
}

#[allow(dead_code)]
fn calibrate_likelihoods(profiles: &[Vec<f64>], coverage: f64) -> Vec<Vec<f64>> {
    let mut counts = vec![0; profiles[0].len()];
    for profile in profiles.iter() {
        for (i, _) in profile.iter().enumerate().filter(|&(_, &x)| 0f64 < x) {
            counts[i] += 1;
        }
    }
    let max_copy_num = (profiles.len() as f64 / coverage).ceil() as usize * 2;
    let likelihoods: Vec<_> = counts
        .iter()
        .map(|&c| max_poisson_lk(c, coverage, 1, max_copy_num))
        .collect();
    let max_lk = likelihoods.iter().fold(f64::NEG_INFINITY, |x, y| x.max(*y));
    // TODO: Can we justify this code? I do not think so.
    let weights: Vec<_> = likelihoods.iter().map(|x| (x - max_lk) / 2f64).collect();
    // let dump: Vec<_> = weights
    //     .iter()
    //     .enumerate()
    //     .map(|(i, x)| format!("{}:{}", i, x))
    //     .collect();
    // debug!("{}", dump.join("\t"));
    profiles
        .iter()
        .map(|xs| xs.iter().zip(weights.iter()).map(|(x, y)| x - y).collect())
        .collect()
}

#[allow(dead_code)]
// Recursively clustering profiles into cluster_num size.
fn recursive_clustering<T: std::borrow::Borrow<[f64]>, R: Rng>(
    profiles: &[T],
    _cluster_num: usize,
    coverage: f64,
    rng: &mut R,
) -> (usize, Vec<u8>) {
    // Filter profiles by doubling.
    // let var_num = 2 * (cluster_num - 1).max(1) as usize;
    let cluster_num = 2;
    let selected_variants: Vec<Vec<_>> = {
        // let probes = filter_profiles(&profiles, cluster_num, 3, coverage);
        let probes = filter_profiles(&profiles, cluster_num, 3, coverage);
        let profiles: Vec<Vec<_>> = profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs.borrow()[pos]).collect())
            .collect();
        profiles
        // let (assignments, _score, cluster_num) = std::iter::repeat(cluster_num)
        //     .take(10)
        //     .map(|k| {
        //         let (asn, score) = mcmc_clustering(&profiles, k, coverage, rng);
        //         let mod_score = score - coverage * (k - 1) as f64;
        //         (asn, mod_score, k)
        //     })
        //     .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
        //     .unwrap();
        // let to_used = retrieve_used_positions(&assignments, &profiles, cluster_num as usize);
        // if to_used.iter().all(|&x| !x) {
        //     return (1, vec![0; profiles.len()]);
        // }
        // profiles
        //     .iter()
        //     .map(|xs| {
        //         xs.iter()
        //             .zip(to_used.iter())
        //             .filter_map(|(&x, &y)| y.then(|| x))
        //             .collect()
        //     })
        //     .collect()
    };
    let (assignments, score, cluster_num) = std::iter::repeat(cluster_num)
        .take(10)
        .filter_map(|k| {
            let (asn, score) = (0..10)
                .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
                .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())?;
            let expt_gain = 1.5 * coverage;
            // debug!("{}\t{}\t{:.2}\t{:.2}", profiles.len(), k, score, expt_gain);
            Some((asn, score - expt_gain, k))
        })
        .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
        .unwrap();
    // for (i, prf) in assignments.iter().zip(selected_variants.iter()) {
    //     let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
    //     debug!("{}\t{}", i, prf.join("\t"));
    // }
    if cluster_num == 1 || score < 0f64 {
        (1, vec![0; profiles.len()])
    } else {
        // Recursively call clustering method on each cluster.
        let mut current_cluster_size = 0;
        let mut final_clustering = vec![0; profiles.len()];
        for cl in 0..cluster_num {
            let (indices, profiles): (Vec<_>, Vec<_>) = profiles
                .iter()
                .enumerate()
                .filter(|&(i, _)| assignments[i] == cl as u8)
                .map(|(i, x)| (i, x.borrow()))
                .unzip();
            if !profiles.is_empty() {
                let (sub_cluster_size, asns) = recursive_clustering(&profiles, 2, coverage, rng);
                for (cl, idx) in asns.iter().zip(indices) {
                    final_clustering[idx] = cl + current_cluster_size as u8;
                }
                current_cluster_size += sub_cluster_size;
            }
        }
        (current_cluster_size, final_clustering)
    }
}

#[allow(dead_code)]
// i->true if there's a cluster using the i-th position to improve the total likelihood.
fn retrieve_used_positions<T: std::borrow::Borrow<[f64]>>(
    assignments: &[u8],
    profiles: &[T],
    cluster_num: usize,
) -> Vec<bool> {
    let len = profiles[0].borrow().len();
    let mut lks = vec![vec![0f64; len]; cluster_num];
    for (&asn, xs) in assignments.iter().zip(profiles.iter()) {
        for (l, &x) in lks[asn as usize].iter_mut().zip(xs.borrow()) {
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
            let mut modif_table = prof.to_modification_table();
            modif_table.truncate(9 * template.len());
            modif_table.extend(prof.to_deletion_table(6));
            modif_table.extend(prof.to_copy_table(6));
            modif_table.iter_mut().for_each(|x| *x -= lk);
            modif_table
        })
        .collect()
}

// Select round * (cluster_num-1) variants
fn filter_profiles<T: std::borrow::Borrow<[f64]>>(
    profiles: &[T],
    cluster_num: usize,
    round: usize,
    coverage: f64,
) -> Vec<(usize, f64)> {
    let template_len = profiles[0].borrow().len() / 9 - 1;
    // (sum, maximum gain, number of positive element)
    let mut total_improvement = vec![(0f64, 0f64, 0); template_len * 9];
    for prof in profiles.iter().map(|x| x.borrow()) {
        for ((sum, maxgain, count), p) in total_improvement.iter_mut().zip(prof) {
            *sum += p;
            *maxgain += p.max(0f64);
            *count += (0.00001 < *p) as usize;
        }
    }
    // TODO: This filtering might be implemented from recursive clustering,
    // and now we do not use recursive clustering, we do not filter any vriants.
    // Filtering out the position where the sum > 0, because it is apparently flippable. -> Or not?
    let mut probes: Vec<(usize, f64)> = total_improvement
        .into_iter()
        .enumerate()
        // .filter(|&(_, (sum, _, _))| sum < 0.01 * profiles.len() as f64)
        .map(|(pos, (_, maxgain, count))| {
            let max_lk = (1..cluster_num + 1)
                .map(|k| poisson_lk(count, profiles.len() as f64 / k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap();
            (pos, maxgain + max_lk)
        })
        .collect();
    probes.sort_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap());
    probes.reverse();
    // for (i, &(p, _)) in probes.iter().take(10).enumerate() {
    //     for (j, prof) in profiles.iter().enumerate() {
    //         debug!("PROFILE\t{}\t{}\t{}\t{}", i, j, p, prof.borrow()[p]);
    //     }
    // }
    let mut selected_variants = vec![];
    'outer: for _ in 0..round {
        let mut buffer = vec![];
        for _ in 0..(cluster_num).max(2) - 1 {
            let next_var_idx = probes.iter().position(|&(pos, _)| {
                for &(p, _) in buffer.iter() {
                    if is_the_same_bipart(profiles, p, pos, coverage) {
                        return false;
                    }
                    let cosine = cosine_similarity(profiles, p, pos);
                    if (1f64 - cosine.abs()).abs() < 0.001 {
                        return false;
                    }
                }
                for &(p, _) in selected_variants.iter() {
                    let cosine = cosine_similarity(profiles, p, pos);
                    if (1f64 - cosine.abs()).abs() < 0.001 {
                        return false;
                    }
                }
                true
            });
            match next_var_idx {
                Some(next_var_idx) => buffer.push(probes.remove(next_var_idx)),
                None => break 'outer,
            }
        }
        selected_variants.extend(buffer)
    }
    selected_variants
}

fn cosine_similarity<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> f64 {
    let (inner_prod, isumsq, jsumsq) = profiles
        .iter()
        .map(|profile| profile.borrow())
        .map(|profile| (profile[i], profile[j]))
        .fold((0f64, 0f64, 0f64), |(ip, isumsq, jsumsq), (x, y)| {
            (ip + x * y, isumsq + x * x, jsumsq + y * y)
        });
    inner_prod / isumsq.sqrt() / jsumsq.sqrt()
}

// Return true if it is the same bipartition.
fn is_the_same_bipart<T: std::borrow::Borrow<[f64]>>(
    profiles: &[T],
    i: usize,
    j: usize,
    cov: f64,
) -> bool {
    let (error_num, total) = profiles
        .iter()
        .map(|prf| (prf.borrow()[i], prf.borrow()[j]))
        .filter(|(ith, jth)| 0.5 < ith.abs() && 0.5 < jth.abs())
        .fold((0, 0), |(err, tot), (ith, jth)| {
            let is_err = ith.is_sign_positive() != jth.is_sign_positive();
            (err + is_err as usize, tot + 1)
        });
    let error_num = error_num.min(total - error_num);
    let poisson_prob: f64 = (0..2 * total / cov.round() as usize)
        .map(|cp| poisson_lk(error_num, cp as f64 * cov))
        .fold(f64::NEG_INFINITY, |x, y| x.max(y));
    let error_prob: f64 = 0.05;
    let binom_prob =
        error_prob.ln() * error_num as f64 + (1f64 - error_prob).ln() * (total - error_num) as f64;
    let pattern_lk: f64 = (0..error_num)
        .map(|i| ((total - i) as f64).ln() - ((error_num - i) as f64).ln())
        .sum();
    let binom_prob = binom_prob + pattern_lk;
    // let ipos = i / 9;
    // let jpos = j / 9;
    // debug!(
    //     "BIPAT\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
    //     ipos, i, jpos, j, error_num, total, poisson_prob, binom_prob,
    // );
    poisson_prob < binom_prob
}

// Return true if i correlate with j.
#[allow(dead_code)]
fn is_correlate<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> bool {
    fn xlogx(x: f64) -> f64 {
        if x.abs() < 0.00001 {
            0f64
        } else {
            x * x.ln()
        }
    }
    fn partlk(p: f64) -> f64 {
        xlogx(p) + xlogx(1f64 - p)
    }
    let mut counts = [0; 4];
    for prof in profiles.iter().map(|x| x.borrow()) {
        let (ith, jth) = (prof[i], prof[j]);
        if ith.abs() < 0.5 || jth.abs() < 0.5 {
            continue;
        }
        match (ith.is_sign_positive(), jth.is_sign_positive()) {
            (true, true) => counts[0] += 1,
            (true, false) => counts[1] += 1,
            (false, true) => counts[2] += 1,
            (false, false) => counts[3] += 1,
        }
    }
    let total: usize = counts.iter().sum();
    let frac_true = (counts[0] + counts[1]) as f64 / total as f64;
    let error_rate = (counts[1] + counts[2]) as f64 / total as f64;
    let corel_model_lk = (partlk(frac_true) + partlk(error_rate)) * total as f64;
    let full_model_lk: f64 = counts.iter().map(|&x| xlogx(x as f64 / total as f64)).sum();
    let full_model_lk = full_model_lk * total as f64;
    let corel_aic = -2f64 * corel_model_lk + 2f64 * 2f64;
    let full_aic = -2f64 * full_model_lk + 2f64 * 3f64;
    // let (i, ied) = (i / 9, i % 9);
    // let (j, jed) = (j / 9, j % 9);
    // debug!(
    //     "CHISQ\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
    //     i, ied, j, jed, total, corel_aic, full_aic
    // );
    corel_aic < full_aic
}

// Return log Poiss{x|lambda}
fn poisson_lk(x: usize, lambda: f64) -> f64 {
    x as f64 * lambda.ln() - lambda - (1..x + 1).map(|c| (c as f64).ln()).sum::<f64>()
}

#[allow(dead_code)]
// Return max_c log Poiss{x|c * lambda}. 1<=c.
fn max_poisson_lk(x: usize, lambda: f64, c_start: usize, c_end: usize) -> f64 {
    (c_start.max(1)..=c_end)
        .map(|c| poisson_lk(x, lambda * c as f64))
        .fold(f64::NEG_INFINITY, |x, y| x.max(y))
}

// Return log Binom(x|size=size,prob=prob)
// Prob is the probability to success.
fn binom_lk(x: usize, size: usize, prob: f64) -> f64 {
    assert!(prob <= 1.0000001);
    let combin: f64 = (1..size - x + 1)
        .map(|i| ((x as f64 / i as f64) + 1.0).ln())
        .sum();
    combin + x as f64 * prob.ln() + (size - x) as f64 * (1.0 - prob).ln()
}

//
fn max_binom_lk(x: usize, size: usize, prob: f64) -> f64 {
    (1..)
        .take_while(|&c| (prob * c as f64) < 1f64)
        .map(|c| binom_lk(x, size, prob * c as f64))
        .fold(f64::NEG_INFINITY, |x, y| x.max(y))
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
fn mcmc_clustering<R: Rng>(data: &[Vec<f64>], k: usize, cov: f64, rng: &mut R) -> (Vec<u8>, f64) {
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
    let lk = mcmc(data, &mut assignments, k, cov, rng);
    (assignments, lk)
}

// Return the maximum likelihood.
fn mcmc<R: Rng>(data: &[Vec<f64>], assign: &mut [u8], k: usize, cov: f64, rng: &mut R) -> f64 {
    // Look up table to compute the optimal copy number and its log likelihood.
    let binom_lks: Vec<_> = (0..=data.len())
        .map(|x| max_binom_lk(x, data.len(), cov / data.len() as f64))
        .collect();
    // how many instance in a cluster.
    let mut clusters = vec![0; k];
    for &asn in assign.iter() {
        clusters[asn as usize] += 1;
    }
    // Current (un-modified) likelihoods.
    let mut lks = vec![vec![0f64; data[0].len()]; k];
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
    let total = 100_000 * k as usize;
    // MAP estimation.
    let (mut max, mut argmax) = (f64::NEG_INFINITY, vec![]);
    for _t in 0..total {
        let (idx, new) = (rng.gen_range(0..data.len()), rng.gen_range(0..k));
        let (old, xs) = (assign[idx], &data[idx]);
        if old == new as u8 {
            continue;
        }
        let old = old as usize;
        let binom_diff = (binom_lks[clusters[old] - 1] + binom_lks[clusters[new] + 1])
            - (binom_lks[clusters[old]] + binom_lks[clusters[new]]);
        let old_diff: f64 = xs.iter().zip(lks[old].iter()).map(diff_old).sum();
        let new_diff: f64 = xs.iter().zip(lks[new].iter()).map(diff_new).sum();
        let total_diff = binom_diff + old_diff + new_diff;
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
        let lk = lk_sums.iter().sum::<f64>() + clusters.iter().map(|&x| binom_lks[x]).sum::<f64>();
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

// #[derive(Debug, Clone)]
// struct Model {
//     length: usize,
//     cluster: usize,
//     clusters: Vec<Vec<bool>>,
//     fractions: Vec<f64>,
// }

// impl Model {
//     fn new<T, U>(profiles: &[T], weights: &[Vec<f64>], cluster: usize) -> Self
//     where
//         T: std::borrow::Borrow<[f64]>,
//     {
//         let length = profiles[0].borrow().len();
//         assert!(profiles.iter().all(|xs| xs.borrow().len() == length));
//         assert!(weights.iter().all(|xs| xs.len() == cluster));
//         let mut fractions = vec![0f64; cluster];
//         let mut cluster_sum = vec![vec![0f64; length]; cluster];
//         for (profile, weight) in profiles.iter().zip(weights.iter()) {
//             assert!((1f64 - weights.iter().sum::<f64>()).abs() < 0.0001);
//             for (k, w) in weight.iter().enumerate() {
//                 fractions[k] += w;
//                 for (i, x) in profile.borrow().iter().enumerate() {
//                     cluster_sum[k][i] += x;
//                 }
//             }
//         }
//         let clusters: Vec<Vec<_>> = cluster_sum
//             .iter()
//             .map(|xs| xs.iter().map(|x| x.is_sign_positive()).collect())
//             .collect();
//         let total: f64 = fractions.iter().sum();
//         fractions.iter_mut().for_each(|x| *x /= total);
//         Self {
//             length,
//             cluster,
//             clusters,
//             fractions,
//         }
//     }
//     fn update_weights<T, U>(&self, weights: &[&mut Vec<f64>], profiles: &[T])
//     where
//         T: std::borrow::Borrow<[f64]>,
//     {
//         for (ws, profile) in weights.iter_mut().zip(profiles.iter()) {
//             let profile = profile.borrow();
//             for (k, w) in ws.iter_mut().enumerate() {
//                 *w = self.fractions[k].ln()
//                     + profile
//                         .iter()
//                         .zip(self.clusters[k].iter())
//                         .fold(0f64, |acc, (&x, &t)| if t { x + acc } else { acc });
//             }
//             let total = logsumexp(ws);
//             ws.iter_mut().for_each(|w| *w = (*w - total).exp());
//         }
//     }
//     fn update_clusters<T>(&mut self, weights: &[Vec<f64>], profiles: &[T])
//     where
//         T: std::borrow::Borrow<[f64]>,
//     {
//         self.fractions.for_each(|x| *x = 0f64);
//         let mut cluster_sum = vec![vec![0f64; self.length]; self.cluster];
//         for (ws, profile) in weights.iter().zip(profiles.iter()) {
//             let profile = profile.borrow();
//             for (k, w) in ws.iter().enumerate() {
//                 self.fractions[k] += w;
//             }
//         }
//     }
//     fn lk(&self) {}
// }

// fn em_clustering<R: Rng>(
//     data: &[Vec<f64>],
//     assign: &mut [u8],
//     k: usize,
//     _cov: f64,
//     rng: &mut R,
// ) -> f64 {
//     let mut weights = vec![vec![0f64; k]; data.len()];
//     for (i, &asn) in assign.iter().enumerate() {
//         weights[i][asn as usize] = 1f64;
//     }
//     let mut model = Model::new(data, &weights, k);
//     let mut lk = std::f64::NEG_INFINITY;
//     for i in 0.. {
//         model.update_weights(&mut weights, data);
//         model.update_clusters(&weights, data);
//         let next_lk: f64 = data.iter().map(|xs| model.lk(xs)).sum();
//         if next_lk < lk {
//             break;
//         }
//         debug!("EM\t{}\t{}", i, next_lk);
//     }
//     for (asn, ws) in assign.iter_mut().zip(weights.iter()) {
//         let (idx, _) = ws
//             .iter()
//             .enumerate()
//             .max_by(|x, y| (x.1).partial_cmp(&(y.1).unwrap()))
//             .unwrap();
//         *asn = idx as u8;
//     }
//     lk
// }

#[cfg(test)]
mod kmeans_test {
    use super::*;
    #[test]
    fn cosine_similarity_test() {
        {
            let profiles = vec![vec![1f64, 1f64], vec![2f64, 2f64]];
            assert!((1f64 - cosine_similarity(&profiles, 0, 1)).abs() < 0.0001);
        }
        {
            let profiles = vec![vec![1f64, -3f64], vec![1f64, -3f64]];
            let cosine = cosine_similarity(&profiles, 0, 1);
            assert!((-1f64 - cosine).abs() < 0.00001, "{}", cosine);
        }
        {
            let profiles = vec![vec![1f64, 3f64], vec![1f64, -3f64]];
            assert!(cosine_similarity(&profiles, 0, 1).abs() < 0.00001);
        }
        {
            let profiles = vec![vec![0f64, 100f64], vec![1f64, 100f64]];
            assert!(
                (std::f64::consts::FRAC_1_SQRT_2 - cosine_similarity(&profiles, 0, 1).abs()
                    < 0.00001)
            );
        }
    }
}
