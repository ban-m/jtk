//! A small K-means clustering algorithm.
// First and last `MASK_LENGTH` bases would not be considered in variant calling.
// Should be greater than the maximum length of kiley::hmm::guided::COPY_SIZE or DEL_SIZE.
const MASK_LENGTH: usize = 7;
use definitions::ReadType;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct ClusteringConfig {
    read_type: ReadType,
    band_width: usize,
    // Minimum required Likelihood gain.
    gain: f64,
    // Coverage for haploid.
    coverage: f64,
    pub copy_num: u8,
}

impl ClusteringConfig {
    pub fn new(
        band_width: usize,
        copy_num: u8,
        coverage: f64,
        gain: f64,
        read_type: ReadType,
    ) -> Self {
        Self {
            read_type,
            band_width,
            coverage,
            gain,
            copy_num,
        }
    }
}

// Assignments, posterior, likelihood, consensus.
type ClusteringResult = (Vec<usize>, Vec<Vec<f64>>, f64, Vec<u8>);
/// Usual "flat(non-recursive)" clustering. Return assignments, template, and LK.
pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<ClusteringResult> {
    let band = config.band_width;
    let mut template = kiley::ternary_consensus_by_chunk(reads, band);
    let mut hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let mut ops: Vec<_> = reads
        .iter()
        .map(|x| hmm.align(&template, x.borrow(), band).1)
        .collect();
    for t in 1..3 {
        hmm.fit_naive_with(&template, reads, &ops, band / t);
        template = hmm.polish_until_converge_with(&template, reads, &mut ops, band / t);
    }
    let result = clustering_dev(&template, reads, &mut ops, rng, &hmm, config);
    result.map(|(asn, gains, lk, _)| (asn, gains, lk, template))
}

type ClusteringDevResult = (Vec<usize>, Vec<Vec<f64>>, f64, usize);
// pub fn clustering_inner<R: Rng, T: std::borrow::Borrow<[u8]>>(
//     template: &[u8],
//     reads: &[T],
//     ops: &mut [Vec<kiley::Op>],
//     rng: &mut R,
//     hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
//     config: &ClusteringConfig,
// ) -> Option<ClusteringDevResult> {
//     let ClusteringConfig {
//         band_width,
//         copy_num,
//         coverage,
//         gain: average_lk,
//         ..
//     } = *config;
//     let profiles: Vec<Vec<_>> = reads
//         .iter()
//         .zip(ops.iter_mut())
//         .map(|(seq, op)| {
//             let (mut table, lk) = hmm.modification_table(template, seq.borrow(), band_width, op);
//             assert!(table.iter().all(|x| x.is_finite()));
//             assert!(table.iter().all(|x| !x.is_nan()));
//             table.iter_mut().for_each(|x| *x -= lk);
//             table
//         })
//         .collect();
//     let copy_num = copy_num as usize;
//     let probes = filter_profiles(&profiles, copy_num, 3, coverage, template.len());
//     let selected_variants: Vec<_> = profiles
//         .iter()
//         .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
//         .collect();
//     let num = 3;
//     let init_copy_num = copy_num.max(4) - 3;
//     let (asn, _, k) = (init_copy_num..=copy_num)
//         .flat_map(|k| std::iter::repeat(k).take(num))
//         .map(|k| {
//             let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
//             let expected_gain_per_read = calc_expected_gain_per_read(k, average_lk);
//             let expected_gain = expected_gain_per_read * reads.len() as f64;
//             (asn, score - expected_gain, k)
//         })
//         .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
//     let sequencepack = (template, reads, ops.as_ref());
//     let modelpack = (hmm, band_width);
//     let (assignments, score, posterior) =
//         re_eval_clustering(sequencepack, modelpack, &asn, coverage, k);
//     let initial_lk: f64 = reads
//         .iter()
//         .zip(ops.iter())
//         .map(|(seq, ops)| hmm.likelihood_guided(&template, seq.borrow(), ops, band_width))
//         .sum();
//     let score = score - initial_lk;
//     if log_enabled!(log::Level::Trace) {
//         for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
//             let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
//             trace!("ASN\t{}\t{}\t{}\t{}", copy_num, id, i, prf.join("\t"));
//         }
//     }
//     Some((assignments, posterior, score, k))
// }

// fn calc_expected_gain_per_read(k: usize, lk: f64) -> f64 {
//     let kf = k as f64;
//     (1f64 - kf.recip()) * lk
//     // let cons_cluster = (1f64 + (kf - 1f64) * (-lk).exp()).ln();
//     // let other_cluster = (1f64 + (kf - 1f64) * (-lk).exp()).ln() + lk;
//     // let reg_term = kf.ln();
//     // (cons_cluster + (kf - 1.0) * other_cluster) / kf - reg_term
// }

// type HMM = kiley::hmm::guided::PairHiddenMarkovModel;
// fn re_eval_clustering<T: std::borrow::Borrow<[u8]>>(
//     (template, reads, ops): (&[u8], &[T], &[Vec<kiley::Op>]),
//     (hmm, band): (&HMM, usize),
//     assignments: &[usize],
//     _coverage: f64,
//     k: usize,
// ) -> (Vec<usize>, f64, Vec<Vec<f64>>) {
//     let fracs: Vec<f64> = {
//         let mut counts = vec![0; k];
//         for &asn in assignments.iter() {
//             counts[asn as usize] += 1;
//         }
//         counts
//             .iter()
//             .map(|&x| (x as f64).max(0.00001).ln() - (reads.len() as f64).ln())
//             .collect()
//     };
//     // let mut consi = vec![];
//     let likelihood_on_clusters: Vec<Vec<f64>> = fracs
//         .iter()
//         .enumerate()
//         .map(|(cl, fr)| {
//             let mut packed_data: PackedData<T> = reads
//                 .iter()
//                 .zip(ops.to_vec())
//                 .zip(assignments.iter())
//                 .enumerate()
//                 .map(|(i, ((read, ops), asn))| (i, read, ops, *asn))
//                 .collect();
//             packed_data.sort_by_key(|x| (x.3 != cl));
//             let consensus = get_consensus_of(&template, &mut packed_data, hmm, band, cl);
//             packed_data.sort_unstable_by_key(|x| x.0);
//             packed_data
//                 .iter()
//                 .map(|(_, seq, ops, _)| hmm.likelihood_guided(&consensus, seq.borrow(), &ops, band))
//                 .map(|lk| lk + fr)
//                 .collect()
//         })
//         .collect();
//     let mut lk = 0f64;
//     let (posterior, assignments): (Vec<_>, Vec<_>) = reads
//         .iter()
//         .enumerate()
//         .map(|(i, seq)| {
//             let orig_lk = hmm.likelihood_guided(&template, seq.borrow(), &ops[i], band);
//             let mut read_lks: Vec<_> = likelihood_on_clusters.iter().map(|lks| lks[i]).collect();
//             let dump: Vec<_> = read_lks
//                 .iter()
//                 .map(|x| format!("{:.2}", x - orig_lk))
//                 .collect();
//             trace!(
//                 "{k}\t{i}\t{}\t{}\t{}",
//                 assignments[i],
//                 seq.borrow().len(),
//                 dump.join("\t")
//             );
//             let read_lk: f64 = logsumexp(&read_lks);
//             lk += read_lk;
//             read_lks.iter_mut().for_each(|x| *x -= read_lk);
//             let (asn, _) = read_lks
//                 .iter()
//                 .enumerate()
//                 .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
//                 .unwrap();
//             (read_lks, asn)
//         })
//         .unzip();
// if consi.len() > 1 {
//     let ed2k = [
//         kiley::Op::Match,
//         kiley::Op::Ins,
//         kiley::Op::Del,
//         kiley::Op::Mismatch,
//     ];
//     for i in 0..consi.len() {
//         for j in i + 1..consi.len() {
//             let ops = edlib_sys::global(&consi[i], &consi[j]);
//             let ops: Vec<_> = ops.iter().map(|&op| ed2k[op as usize]).collect();
//             let (xr, ar, yr) = kiley::recover(&consi[i], &consi[j], &ops);
//             for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
//                 eprintln!("{}", String::from_utf8_lossy(xr));
//                 eprintln!("{}", String::from_utf8_lossy(ar));
//                 eprintln!("{}\n", String::from_utf8_lossy(yr));
//             }
//         }
//     }
// }
//     (assignments, lk, posterior)
// }

// type PackedData<'a, T> = Vec<(usize, &'a T, Vec<kiley::Op>, usize)>;
// fn get_consensus_of<T: std::borrow::Borrow<[u8]>>(
//     template: &[u8],
//     packed_data: &mut PackedData<T>,
//     hmm: &HMM,
//     band: usize,
//     cl: usize,
// ) -> Vec<u8> {
//     let take = packed_data.iter().filter(|x| x.3 == cl).count();
//     let (sorted, mut ops): (Vec<_>, Vec<_>) = packed_data
//         .iter_mut()
//         .map(|(_, read, ops, _)| (read.borrow(), ops))
//         .unzip();
//     let config = kiley::hmm::guided::HMMConfig::new(band, take, MASK_LENGTH);
//     hmm.polish_until_converge_with_conf(&template, &sorted, &mut ops, &config)
// }

pub fn clustering_dev<R: Rng, T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &mut [Vec<kiley::Op>],
    rng: &mut R,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
    config: &ClusteringConfig,
) -> Option<ClusteringDevResult> {
    trace!("{}", String::from_utf8_lossy(template));
    let ClusteringConfig {
        band_width,
        copy_num,
        coverage,
        gain: average_lk,
        ..
    } = *config;
    let (profiles, _): (Vec<_>, Vec<_>) = reads
        .iter()
        .zip(ops.iter_mut())
        .map(|(seq, op)| {
            let (mut table, lk) = hmm.modification_table(template, seq.borrow(), band_width, op);
            assert!(table.iter().all(|x| x.is_finite()));
            assert!(table.iter().all(|x| !x.is_nan()));
            table.iter_mut().for_each(|x| *x -= lk);
            (table, lk)
        })
        .unzip();
    let copy_num = copy_num as usize;
    let selected_variants: Vec<_> = {
        let probes = filter_profiles(&profiles, copy_num, 3, coverage, template.len());
        for (pos, lk) in probes.iter() {
            let (pos, var) = (
                pos / kiley::hmm::guided::NUM_ROW,
                pos % kiley::hmm::guided::NUM_ROW,
            );
            trace!("DUMP\t{}\t{}\t{}", pos, var, lk);
        }
        profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect()
    };
    let init_copy_num = copy_num.max(4) - 3;
    let (assignments, _score, _k) = (init_copy_num..=copy_num)
        .map(|k| {
            let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
            let expected_gain =
                average_lk * number_of_improved_read(&selected_variants, &asn, k) as f64;
            trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t1");
            (asn, score - expected_gain, k)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    if log_enabled!(log::Level::Trace) {
        trace!("VARS\tCOPYNUM\tID\tASN\tPOS\tLKDiff\tFilter");
        for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
            for (pos, x) in prf.iter().enumerate() {
                trace!("VARS\t{copy_num}\t{id}\t{i}\t{pos}\t{x}\tBefore");
            }
        }
    }
    let selected_variants = filter_suspicious_variants(&selected_variants, &assignments);
    let (assignments, score, k) = (init_copy_num..=copy_num)
        .map(|k| {
            let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
            let expected_gain =
                average_lk * number_of_improved_read(&selected_variants, &asn, k) as f64;
            trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t1");
            (asn, score - expected_gain, k)
        })
        .fold((vec![0; reads.len()], 0f64, 1), |x, y| {
            match x.1.partial_cmp(&y.1).unwrap() {
                std::cmp::Ordering::Less => y,
                _ => x,
            }
        });
    let mut likelihood_gains = get_likelihood_gain(&selected_variants, &assignments);
    to_posterior_probability(&mut likelihood_gains);
    if log_enabled!(log::Level::Trace) {
        for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
            for (pos, x) in prf.iter().enumerate() {
                trace!("VARS\t{copy_num}\t{id}\t{i}\t{pos}\t{x}\tAfter");
            }
        }
        // for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
        //     let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
        //     trace!("ASN\t{}\t{}\t{}\t{}", copy_num, id, i, prf.join("\t"));
        // }
    }
    trace!("DUMP\t{score}\t{k}");
    likelihood_gains.iter().for_each(|x| assert_eq!(x.len(), k));
    Some((assignments, likelihood_gains, score, k))
}

// LK->LK-logsumexp(LK).
fn to_posterior_probability(lks: &mut [Vec<f64>]) {
    for xs in lks.iter_mut() {
        let total = logsumexp(xs);
        xs.iter_mut().for_each(|x| *x -= total);
    }
}

fn number_of_improved_read(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> usize {
    let dim = variants[0].len();
    let mut counts = vec![0; k];
    let mut lks: Vec<_> = vec![vec![0f64; dim]; k];
    for (&asn, xs) in assignments.iter().zip(variants.iter()) {
        counts[asn] += 1;
        lks[asn]
            .iter_mut()
            .zip(xs.iter())
            .for_each(|(x, y)| *x += y);
    }
    counts
        .iter()
        .zip(lks)
        .filter_map(|(count, sums)| sums.iter().any(|x| x.is_sign_positive()).then(|| count))
        .sum()
}

// i->k->the likelihood gain of the i-th read when clustered in the k-th cluster.
// `copy_num` is the copy number of this unit, not the *cluster number*.
// Usually, they are the same but sometimes there are exact repeats, and the number of the cluster
// would be smaller than the copy number.
fn get_likelihood_gain(variants: &[Vec<f64>], assignments: &[usize]) -> Vec<Vec<f64>> {
    let cluster_size = *assignments.iter().max().unwrap_or(&0) + 1;
    let mut total_lk_gain = vec![vec![0f64; variants[0].len()]; cluster_size];
    let mut count = vec![0.000000000001; cluster_size];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        count[asn] += 1.0;
        for (total, var) in total_lk_gain[asn].iter_mut().zip(vars.iter()) {
            *total += var;
        }
    }
    let is_used_position: Vec<Vec<_>> = total_lk_gain
        .iter()
        .map(|totals| totals.iter().map(|x| x.is_sign_positive()).collect())
        .collect();
    let len: f64 = count.iter().sum();
    let fractions: Vec<_> = count.iter().map(|&x| (x as f64 / len)).collect();
    // Padding -infinity to make sure that there are exactly k-slots.
    // This is mainly because we do not want to these `mock` slots.
    // Do NOT use std:;f64::NEG_INFINITY as it must cause under flows.
    // let pad = std::iter::repeat(SMALL_LK).take(copy_num.saturating_sub(cluster_size));
    fn pick_vars(ps: &[bool], vars: &[f64]) -> f64 {
        vars.iter()
            .zip(ps)
            .filter_map(|(lk, &p)| p.then(|| lk))
            .sum()
    }
    variants
        .iter()
        .map(|vars| {
            is_used_position
                .iter()
                .zip(fractions.iter())
                .map(|(ps, f)| f.ln() + pick_vars(ps, vars))
                .collect()
        })
        .collect()
}

// Filter column with small confidence.
// HowTo: If a column is indeed a variant site, in at least one cluster,
// the column is almost filled with positive values.
// More precisely, given the error rate e, the # of reads with the negative value
// would follow the poisson distribution with parameter lambda = |r| * e, where
// |r| is the # of the read in that cluster.
// If this value is too small, we can remove that column from `variants.`
// Also, we check "exhaution criteria", in other words,
// If we use a variant, we use that variant in almost all the reads.
const FRAC: f64 = 0.70;
const IN_POS_RATIO: f64 = 3f64;
fn filter_suspicious_variants(variants: &[Vec<f64>], assignments: &[usize]) -> Vec<Vec<f64>> {
    let max_asn = *assignments.iter().max().unwrap() + 1;
    let dim = variants[0].len();
    // Assignment -> variant position -> (total #, # of positive element, likelihood sum).
    let mut var_summary = vec![vec![(0, 0, 0f64); dim]; max_asn];
    for (&asn, xs) in assignments.iter().zip(variants.iter()) {
        assert_eq!(var_summary[asn].len(), xs.len());
        for (slot, x) in var_summary[asn].iter_mut().zip(xs.iter()) {
            slot.0 += 1;
            slot.1 += x.is_sign_positive() as u32;
            slot.2 += x;
        }
    }
    // We cound the # of the positive elements used and # of positive elm not used.
    let scattered_vars: Vec<_> = (0..dim)
        .map(|d| {
            let (in_pos, in_neg) = var_summary.iter().map(|vars| &vars[d]).fold(
                (0, 0),
                |(pos, neg), &(_, numpos, sum)| match sum.is_sign_positive() {
                    true => (pos + numpos, neg),
                    false => (pos, neg + numpos),
                },
            );
            IN_POS_RATIO < (in_pos as f64 / in_neg as f64)
        })
        .collect();
    let used_pos = var_summary
        .iter()
        .fold(vec![false; dim], |mut used_pos, cs| {
            for (pos, &(total, num_pos, lksum)) in used_pos.iter_mut().zip(cs) {
                *pos |= (0f64 < lksum) & (total as f64 * FRAC < num_pos as f64);
            }
            used_pos
        });
    let now = used_pos
        .iter()
        .zip(scattered_vars.iter())
        .filter(|&(&u, &s)| u & s)
        .count();
    trace!("FILTER\t{}\t{}", dim, now);
    variants
        .iter()
        .map(|xs| {
            xs.iter()
                .zip(used_pos.iter())
                .zip(scattered_vars.iter())
                .filter_map(|((&x, &u), &s)| (s & u).then(|| x))
                .collect()
        })
        .collect()
}

// #[allow(dead_code)]
// fn calibrate_likelihoods(profiles: &[Vec<f64>], coverage: f64) -> Vec<Vec<f64>> {
//     let mut counts = vec![0; profiles[0].len()];
//     for profile in profiles.iter() {
//         for (i, _) in profile.iter().enumerate().filter(|&(_, &x)| 0f64 < x) {
//             counts[i] += 1;
//         }
//     }
//     let max_copy_num = (profiles.len() as f64 / coverage).ceil() as usize * 2;
//     let likelihoods: Vec<_> = counts
//         .iter()
//         .map(|&c| max_poisson_lk(c, coverage, 1, max_copy_num))
//         .collect();
//     let max_lk = likelihoods.iter().fold(f64::NEG_INFINITY, |x, y| x.max(*y));
// let weights: Vec<_> = likelihoods.iter().map(|x| (x - max_lk) / 2f64).collect();
// let dump: Vec<_> = weights
//     .iter()
//     .enumerate()
//     .map(|(i, x)| format!("{}:{}", i, x))
//     .collect();
// debug!("{}", dump.join("\t"));
//     profiles
//         .iter()
//         .map(|xs| xs.iter().zip(weights.iter()).map(|(x, y)| x - y).collect())
//         .collect()
// }

// #[allow(dead_code)]
// Recursively clustering profiles into cluster_num size.
// fn recursive_clustering<T: std::borrow::Borrow<[f64]>, R: Rng>(
//     template_len: usize,
//     profiles: &[T],
//     _cluster_num: usize,
//     coverage: f64,
//     rng: &mut R,
// ) -> (usize, Vec<u8>) {
//     // Filter profiles by doubling.
//     // let var_num = 2 * (cluster_num - 1).max(1) as usize;
//     let cluster_num = 2;
//     let selected_variants: Vec<Vec<_>> = {
//         let probes = filter_profiles(profiles, cluster_num, 3, coverage, template_len);
//         let profiles: Vec<Vec<_>> = profiles
//             .iter()
//             .map(|xs| probes.iter().map(|&(pos, _)| xs.borrow()[pos]).collect())
//             .collect();
//         profiles
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
//     };
//     let (assignments, score, cluster_num) = std::iter::repeat(cluster_num)
//         .take(10)
//         .filter_map(|k| {
//             let (asn, score) = (0..10)
//                 .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
//                 .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())?;
//             let expt_gain = 1.5 * coverage;
//             // debug!("{}\t{}\t{:.2}\t{:.2}", profiles.len(), k, score, expt_gain);
//             Some((asn, score - expt_gain, k))
//         })
//         .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
//         .unwrap();
//     // for (i, prf) in assignments.iter().zip(selected_variants.iter()) {
//     //     let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
//     //     debug!("{}\t{}", i, prf.join("\t"));
//     // }
//     if cluster_num == 1 || score < 0f64 {
//         (1, vec![0; profiles.len()])
//     } else {
//         // Recursively call clustering method on each cluster.
//         let mut current_cluster_size = 0;
//         let mut final_clustering = vec![0; profiles.len()];
//         for cl in 0..cluster_num {
//             let (indices, profiles): (Vec<_>, Vec<_>) = profiles
//                 .iter()
//                 .enumerate()
//                 .filter(|&(i, _)| assignments[i] == cl as u8)
//                 .map(|(i, x)| (i, x.borrow()))
//                 .unzip();
//             if !profiles.is_empty() {
//                 let (sub_cluster_size, asns) =
//                     recursive_clustering(template_len, &profiles, 2, coverage, rng);
//                 for (cl, idx) in asns.iter().zip(indices) {
//                     final_clustering[idx] = cl + current_cluster_size as u8;
//                 }
//                 current_cluster_size += sub_cluster_size;
//             }
//         }
//         (current_cluster_size, final_clustering)
//     }
// }
// i->true if there's a cluster using the i-th position to improve the total likelihood.
// fn retrieve_used_positions<T: std::borrow::Borrow<[f64]>>(
//     assignments: &[u8],
//     profiles: &[T],
//     cluster_num: usize,
// ) -> Vec<bool> {
//     let len = profiles[0].borrow().len();
//     let mut lks = vec![vec![0f64; len]; cluster_num];
//     for (&asn, xs) in assignments.iter().zip(profiles.iter()) {
//         for (l, &x) in lks[asn as usize].iter_mut().zip(xs.borrow()) {
//             *l += x;
//         }
//     }
//     lks.iter().fold(vec![false; len], |mut acc, xs| {
//         for (is_used, total_lk) in acc.iter_mut().zip(xs) {
//             *is_used |= total_lk.is_sign_positive();
//         }
//         acc
//     })
// }

// fn get_profiles<T: std::borrow::Borrow<[u8]>>(
//     template: &[u8],
//     hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
//     reads: &[T],
//     band_width: isize,
// ) -> Option<Vec<Vec<f64>>> {
//     use kiley::gphmm::*;
//     let template = kiley::padseq::PadSeq::new(template);
//     let mut profiles = vec![];
//     for read in reads.iter() {
//         let read = kiley::padseq::PadSeq::new(read.borrow());
//         let prof = match banded::ProfileBanded::new(hmm, &template, &read, band_width) {
//             Some(res) => res,
//             None => {
//                 for read in reads.iter() {
//                     let read = String::from_utf8(read.borrow().to_vec()).unwrap();
//                     error!("READ\tKMEANS\t{}", read);
//                 }
//                 let template = String::from_utf8(template.clone().into()).unwrap();
//                 error!("TEMPLATE\tKMEANS\t{}", template);
//                 error!("{}", hmm);
//                 return None;
//             }
//         };
//         let lk = prof.lk();
//         // Modif.
//         let mut modif_table = prof.to_modification_table();
//         modif_table.truncate(9 * template.len());
//         // Del table
//         let del_table = prof.to_deletion_table(DEL_SIZE);
//         assert_eq!(
//             del_table.len(),
//             (DEL_SIZE - 1) * (template.len() - DEL_SIZE)
//         );
//         modif_table.extend(del_table);
//         // Copy Table.
//         let copy_table = prof.to_copy_table(REP_SIZE);
//         assert_eq!(copy_table.len(), REP_SIZE * (template.len() - REP_SIZE));
//         modif_table.extend(copy_table);
//         // Normalize.
//         modif_table.iter_mut().for_each(|x| *x -= lk);
//         profiles.push(modif_table);
//     }
//     Some(profiles)
// }

// Select round * (cluster_num-1) variants
// fn filter_profiles<T: std::borrow::Borrow<[f64]>>(
//     profiles: &[T],
//     cluster_num: usize,
//     round: usize,
//     coverage: f64,
//     template_len: usize,
// ) -> Vec<(usize, f64)> {
//     // (sum, maximum gain, number of positive element)
//     let mut total_improvement = vec![(0f64, 0); profiles[0].borrow().len()];
//     for prof in profiles.iter().map(|x| x.borrow()) {
//         for ((maxgain, count), p) in total_improvement.iter_mut().zip(prof) {
//             *maxgain += p.max(0f64);
//             *count += (0.00001 < *p) as usize;
//         }
//     }
//     let base_position = |pos: usize| match pos / (9 * template_len) {
//         0 => pos / 9,
//         _ => {
//             let pos = pos - 9 * template_len;
//             match pos / ((DEL_SIZE - 1) * (template_len - DEL_SIZE)) {
//                 0 => pos / (DEL_SIZE - 1),
//                 _ => (pos - (DEL_SIZE - 1) * (template_len - DEL_SIZE)) / REP_SIZE,
//             }
//         }
//     };
//     let in_mask = |pos: usize| {
//         let pos = base_position(pos);
//         pos < MASK_LENGTH || (template_len - MASK_LENGTH) < pos
//     };
//     // trace!("LKSUM\tpos\tbp\tlen\ttype\tlkdelta");
//     let probes: Vec<(usize, f64)> = total_improvement
//         .into_iter()
//         .enumerate()
//         .filter(|&(pos, _)| !in_mask(pos))
//         .map(|(pos, (maxgain, count))| {
//             let max_lk = (1..cluster_num + 1)
//                 .map(|k| poisson_lk(count, coverage * k as f64))
//                 .max_by(|x, y| x.partial_cmp(y).unwrap())
//                 .unwrap_or_else(|| panic!("{}", cluster_num));
//             let total_lk = max_lk + maxgain;
//             (pos, total_lk)
//         })
//         .collect();
//     // 0-> not selected yet. 1-> selected. 2-> removed because of co-occurence.
//     // Currently, find and sweep. So, to select one variant,
//     // It takes O(L) time to find a variant, and O(L) time to sweep linked variant.
//     // So, in total, O(ML) time to select M variants. It is OK I think, because
//     // usually L is 2K, M is around 3-20.
//     let mut is_selected = vec![0; probes.len()];
//     'outer: for _ in 0..round {
//         let mut weights: Vec<_> = vec![1f64; probes.len()];
//         for _ in 0..cluster_num.max(2) - 1 {
//             let next_var_idx = probes
//                 .iter()
//                 .map(|x| x.1)
//                 .zip(weights.iter())
//                 .zip(is_selected.iter())
//                 .enumerate()
//                 .max_by(|(_, ((lk1, w1), picked1)), (_, ((lk2, w2), picked2))| {
//                     match (picked1, picked2) {
//                         (1, _) | (2, _) => std::cmp::Ordering::Less,
//                         (_, 1) | (_, 2) => std::cmp::Ordering::Greater,
//                         _ => (lk1 * *w1).partial_cmp(&(lk2 * *w2)).unwrap(),
//                     }
//                 });
//             if let Some((next_var_idx, _)) = next_var_idx {
//                 let picked_pos = probes[next_var_idx].0;
//                 let picked_pos_in_bp = base_position(picked_pos);
//                 is_selected[next_var_idx] = 1;
//                 for ((&(pos, _), weight), selected) in probes
//                     .iter()
//                     .zip(weights.iter_mut())
//                     .zip(is_selected.iter_mut())
//                     .filter(|&(_, &mut selected)| selected == 0)
//                 {
//                     let pos_in_bp = base_position(pos);
//                     let diff_in_bp =
//                         pos_in_bp.max(picked_pos_in_bp) - pos_in_bp.min(picked_pos_in_bp);
//                     let sim = sokal_michener(profiles, picked_pos, pos);
//                     let cos_sim = cosine_similarity(profiles, picked_pos, pos);
//                     // Note that, 1/40 = 0.975. So, if the 39 sign out of 40 pair is positive, then,
//                     // the sokal michener's similarity would be 0.975.
//                     if 0.95 < sim || 1f64 - cos_sim.abs() < 0.01 || diff_in_bp < MASK_LENGTH {
//                         *selected = 2;
//                     }
//                     if 0.75 < sim {
//                         *weight *= 1f64 - sim;
//                     }
//                 }
//             } else {
//                 break 'outer;
//             }
//         }
//     }
//     let selected_variants: Vec<_> = probes
//         .into_iter()
//         .zip(is_selected)
//         .filter_map(|(x, y)| (y == 1).then(|| x))
//         .collect();
//     selected_variants
// }

fn filter_profiles<T: std::borrow::Borrow<[f64]>>(
    profiles: &[T],
    cluster_num: usize,
    round: usize,
    coverage: f64,
    template_len: usize,
) -> Vec<(usize, f64)> {
    // (sum, maximum gain, number of positive element)
    let mut total_improvement = vec![(0f64, 0); profiles[0].borrow().len()];
    for prof in profiles.iter().map(|x| x.borrow()) {
        for ((maxgain, count), p) in total_improvement.iter_mut().zip(prof) {
            *maxgain += p.max(0f64);
            *count += (0.00001 < *p) as usize;
        }
    }
    let base_position = |pos: usize| pos / kiley::hmm::guided::NUM_ROW;
    let in_mask = |pos: usize| {
        let pos = base_position(pos);
        pos < MASK_LENGTH || (template_len - MASK_LENGTH) < pos
    };
    let probes: Vec<(usize, f64)> = total_improvement
        .into_iter()
        .enumerate()
        .filter(|&(pos, _)| !in_mask(pos))
        .map(|(pos, (maxgain, count))| {
            let max_lk = (1..cluster_num + 1)
                .map(|k| poisson_lk(count, coverage * k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap_or_else(|| panic!("{}", cluster_num));
            let total_lk = max_lk + maxgain;
            (pos, total_lk)
        })
        .collect();
    // 0-> not selected yet. 1-> selected. 2-> removed because of co-occurence.
    // Currently, find and sweep. So, to select one variant,
    // It takes O(L) time to find a variant, and O(L) time to sweep linked variant.
    // So, in total, O(ML) time to select M variants. It is OK I think, because
    // usually L is 2K, M is around 3-20.
    let mut is_selected = vec![0; probes.len()];
    'outer: for _ in 0..round {
        let mut weights: Vec<_> = vec![1f64; probes.len()];
        for _ in 0..cluster_num.max(2) - 1 {
            let next_var_idx = probes
                .iter()
                .map(|x| x.1)
                .zip(weights.iter())
                .zip(is_selected.iter())
                .enumerate()
                .max_by(|(_, ((lk1, w1), picked1)), (_, ((lk2, w2), picked2))| {
                    match (picked1, picked2) {
                        (1, _) | (2, _) => std::cmp::Ordering::Less,
                        (_, 1) | (_, 2) => std::cmp::Ordering::Greater,
                        _ => (lk1 * *w1).partial_cmp(&(lk2 * *w2)).unwrap(),
                    }
                });
            if let Some((next_var_idx, _)) = next_var_idx {
                let picked_pos = probes[next_var_idx].0;
                let picked_pos_in_bp = base_position(picked_pos);
                is_selected[next_var_idx] = 1;
                for ((&(pos, _), weight), selected) in probes
                    .iter()
                    .zip(weights.iter_mut())
                    .zip(is_selected.iter_mut())
                    .filter(|&(_, &mut selected)| selected == 0)
                {
                    let pos_in_bp = base_position(pos);
                    let diff_in_bp =
                        pos_in_bp.max(picked_pos_in_bp) - pos_in_bp.min(picked_pos_in_bp);
                    let sim = sokal_michener(profiles, picked_pos, pos);
                    let cos_sim = cosine_similarity(profiles, picked_pos, pos);
                    // Note that, 1/40 = 0.975. So, if the 39 sign out of 40 pair is positive, then,
                    // the sokal michener's similarity would be 0.975.
                    if 0.95 < sim || 1f64 - cos_sim.abs() < 0.01 || diff_in_bp < MASK_LENGTH {
                        *selected = 2;
                    }
                    if 0.75 < sim {
                        *weight *= 1f64 - sim;
                    }
                }
            } else {
                break 'outer;
            }
        }
    }
    let selected_variants: Vec<_> = probes
        .into_iter()
        .zip(is_selected)
        .filter_map(|(x, y)| (y == 1).then(|| x))
        .collect();
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
// fn is_the_same_bipart<T: std::borrow::Borrow<[f64]>>(
//     profiles: &[T],
//     i: usize,
//     j: usize,
//     cov: f64,
// ) -> bool {
//     let (error_num, total) = profiles
//         .iter()
//         .map(|prf| (prf.borrow()[i], prf.borrow()[j]))
//         .filter(|(ith, jth)| 0.5 < ith.abs() && 0.5 < jth.abs())
//         .fold((0, 0), |(err, tot), (ith, jth)| {
//             let is_err = ith.is_sign_positive() != jth.is_sign_positive();
//             (err + is_err as usize, tot + 1)
//         });
//     let error_num = error_num.min(total - error_num);
//     let poisson_prob: f64 = (0..2 * total / cov.round() as usize)
//         .map(|cp| poisson_lk(error_num, cp as f64 * cov))
//         .fold(f64::NEG_INFINITY, |x, y| x.max(y));
//     let error_prob: f64 = 0.05;
//     let binom_prob =
//         error_prob.ln() * error_num as f64 + (1f64 - error_prob).ln() * (total - error_num) as f64;
//     let pattern_lk: f64 = (0..error_num)
//         .map(|i| ((total - i) as f64).ln() - ((error_num - i) as f64).ln())
//         .sum();
//     let binom_prob = binom_prob + pattern_lk;
//     poisson_prob < binom_prob
// }

// Return Max-Sokal-Michener similarity of the i-th column and j-th (flipped) column.
fn sokal_michener<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> f64 {
    let (matches, mismatches) = profiles
        .iter()
        .map(|prof| {
            let prof = prof.borrow();
            (prof[i] * prof[j]).is_sign_positive()
        })
        .fold((0, 0), |(mat, mism), x| match x {
            true => (mat + 1, mism),
            false => (mat, mism + 1),
        });
    mismatches.max(matches) as f64 / profiles.len() as f64
}

// fn sokal_michener_sign<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> f64 {
//     let matches = profiles
//         .iter()
//         .filter(|prof| {
//             let prof = prof.borrow();
//             (prof[i] * prof[j]).is_sign_positive()
//         })
//         .count();
//     matches as f64 / profiles.len() as f64
// }

// Return true if i correlate with j.
// fn is_correlate<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> bool {
//     fn xlogx(x: f64) -> f64 {
//         if x.abs() < 0.00001 {
//             0f64
//         } else {
//             x * x.ln()
//         }
//     }
//     fn partlk(p: f64) -> f64 {
//         xlogx(p) + xlogx(1f64 - p)
//     }
//     let mut counts = [0; 4];
//     for prof in profiles.iter().map(|x| x.borrow()) {
//         let (ith, jth) = (prof[i], prof[j]);
//         if ith.abs() < 0.5 || jth.abs() < 0.5 {
//             continue;
//         }
//         match (ith.is_sign_positive(), jth.is_sign_positive()) {
//             (true, true) => counts[0] += 1,
//             (true, false) => counts[1] += 1,
//             (false, true) => counts[2] += 1,
//             (false, false) => counts[3] += 1,
//         }
//     }
//     let total: usize = counts.iter().sum();
//     let frac_true = (counts[0] + counts[1]) as f64 / total as f64;
//     let error_rate = (counts[1] + counts[2]) as f64 / total as f64;
//     let corel_model_lk = (partlk(frac_true) + partlk(error_rate)) * total as f64;
//     let full_model_lk: f64 = counts.iter().map(|&x| xlogx(x as f64 / total as f64)).sum();
//     let full_model_lk = full_model_lk * total as f64;
//     let corel_aic = -2f64 * corel_model_lk + 2f64 * 2f64;
//     let full_aic = -2f64 * full_model_lk + 2f64 * 3f64;
//     // let (i, ied) = (i / 9, i % 9);
//     // let (j, jed) = (j / 9, j % 9);
//     // debug!(
//     //     "CHISQ\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
//     //     i, ied, j, jed, total, corel_aic, full_aic
//     // );
//     corel_aic < full_aic
// }

// fn p2p(pos: usize, template_len: usize) -> usize {
//     let idx = pos / 9;
//     if idx < template_len {
//         idx
//     } else {
//         let idx = pos - 9 * template_len;
//         if idx < (DEL_SIZE - 1) * (template_len - DEL_SIZE) {
//             idx / (DEL_SIZE - 1)
//         } else {
//             let idx = idx - (DEL_SIZE - 1) * (template_len - DEL_SIZE);
//             idx / REP_SIZE
//         }
//     }
// }

// fn feature_extract<R: Rng>(
//     profiles: &[Vec<f64>],
//     cluster_num: usize,
//     coverage: f64,
//     template_len: usize,
//     rng: &mut R,
// ) -> Vec<Vec<f64>> {
//     let mut total_improvement = vec![(0f64, 0); profiles[0].len()];
//     for prof in profiles.iter() {
//         for ((maxgain, count), p) in total_improvement.iter_mut().zip(prof) {
//             *maxgain += p.max(0f64);
//             *count += (0.00001 < *p) as usize;
//         }
//     }
//     let mut probes: Vec<(usize, f64)> = total_improvement
//         .iter()
//         .enumerate()
//         .map(|(pos, &(maxgain, count))| {
//             let max_lk = (1..cluster_num + 1)
//                 .map(|k| poisson_lk(count, coverage * k as f64))
//                 .max_by(|x, y| x.partial_cmp(y).unwrap())
//                 .unwrap_or_else(|| panic!("{}", cluster_num));
//             (pos, max_lk + maxgain)
//         })
//         .filter(|&(_, lk)| lk.is_sign_positive())
//         .collect();
//     probes.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap());
//     // There's no feature currently.
//     let mut features: Vec<Vec<f64>> = vec![vec![]; profiles.len()];
//     let mut max_lk_gains = vec![];
//     for _ in 0..3 * (cluster_num.max(2) - 1) {
//         // Pop the current maximum variants.
//         let (pos, mut lk_gain) = match probes.pop() {
//             Some(res) => res,
//             None => break,
//         };
//         let pos_in_bp = p2p(pos, template_len);
//         trace!("SEED\t{}\t{}", pos_in_bp, lk_gain);
//         let mut feature: Vec<_> = profiles.iter().map(|fs| fs[pos]).collect();
//         let mut pushed_pos = vec![pos_in_bp];
//         probes.retain(|&(pos2, max_lk)| {
//             let pos2_in_bp = p2p(pos2, template_len);
//             let is_the_same = 0.99 < cosine_similarity(profiles, pos, pos2);
//             let sokal_mich_sim = sokal_michener_sign(profiles, pos, pos2);
//             let is_near = pushed_pos
//                 .iter()
//                 .any(|&pos_bp| pos_bp.max(pos2_in_bp) - pos_bp.min(pos2_in_bp) < MASK_LENGTH);
//             if is_the_same {
//                 false
//             } else if 0.95 < sokal_mich_sim {
//                 if !is_near {
//                     feature
//                         .iter_mut()
//                         .zip(profiles.iter().map(|fs| fs[pos2]))
//                         .for_each(|(x, y)| *x += y);
//                     trace!("FOLLOW\t{}\t+\t{}", p2p(pos2, template_len), max_lk);
//                     lk_gain += max_lk;
//                 }
//                 pushed_pos.push(pos2_in_bp);
//                 false
//             } else if sokal_mich_sim < 0.05 {
//                 if !is_near {
//                     feature
//                         .iter_mut()
//                         .zip(profiles.iter().map(|fs| fs[pos2]))
//                         .for_each(|(x, y)| *x -= y);
//                     trace!("FOLLOW\t{}\t-\t{}", p2p(pos2, template_len), max_lk);
//                     lk_gain += max_lk;
//                 }
//                 pushed_pos.push(pos2_in_bp);
//                 false
//             } else {
//                 true
//             }
//         });
//         for (fs, &x) in features.iter_mut().zip(feature.iter()) {
//             fs.push(x);
//         }
//         max_lk_gains.push(lk_gain);
//     }
//     let gains = try_clustering(&features, cluster_num, coverage, rng);
//     trace!("{:?}", max_lk_gains);
//     trace!("{:?}", gains);
//     let is_ok_position: Vec<_> = gains
//         .iter()
//         .zip(max_lk_gains.iter())
//         .map(|(&gain, &upperbound)| upperbound / 2f64 < gain)
//         .collect();
//     for fs in features.iter_mut() {
//         let mut idx = 0;
//         fs.retain(|_| {
//             idx += 1;
//             is_ok_position[idx - 1]
//         });
//     }
//     features
// }

// fn try_clustering<R: Rng>(
//     data: &[Vec<f64>],
//     cluster_num: usize,
//     coverage: f64,
//     rng: &mut R,
// ) -> Vec<f64> {
//     let (asn, _): (Vec<_>, _) = (0..2)
//         .map(|_| mcmc_clustering(data, cluster_num, coverage, rng))
//         .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
//         .unwrap();
//     let mut gains = vec![vec![0f64; data[0].len()]; cluster_num];
//     for (asn, prof) in asn.iter().zip(data.iter()) {
//         assert_eq!(prof.len(), gains[*asn as usize].len());
//         for (slot, x) in gains[*asn as usize].iter_mut().zip(prof.iter()) {
//             *slot += x;
//         }
//     }
//     (0..data[0].len())
//         .map(|pos| gains.iter().map(|lks| lks[pos].max(0f64)).sum::<f64>())
//         .collect()
// }

// Return log Poiss{x|lambda}
fn poisson_lk(x: usize, lambda: f64) -> f64 {
    x as f64 * lambda.ln() - lambda - (1..x + 1).map(|c| (c as f64).ln()).sum::<f64>()
}

// Return max_c log Poiss{x|c * lambda}. 1<=c.
fn max_poisson_lk(x: usize, lambda: f64, c_start: usize, c_end: usize) -> f64 {
    (c_start.max(1)..=c_end)
        .map(|c| poisson_lk(x, lambda * c as f64))
        .fold(f64::NEG_INFINITY, |x, y| x.max(y))
}

// Return log Binom(x|size=size,prob=prob)
// Prob is the probability to success.
// fn binom_lk(x: usize, size: usize, prob: f64) -> f64 {
//     assert!(prob <= 1.0000001);
//     let combin: f64 = (1..size - x + 1)
//         .map(|i| ((x as f64 / i as f64) + 1.0).ln())
//         .sum();
//     combin + x as f64 * prob.ln() + (size - x) as f64 * (1.0 - prob).ln()
// }

// fn max_binom_lk(x: usize, size: usize, prob: f64) -> f64 {
//     (1..)
//         .take_while(|&c| (prob * c as f64) < 1f64)
//         .map(|c| binom_lk(x, size, prob * c as f64))
//         .fold(f64::NEG_INFINITY, |x, y| x.max(y))
// }

// // Likelihood of each size of the clsuters.
// fn poisson_likelihood(asn: &[u8], cluster_num: u8) -> f64 {
//     let mut count = vec![0; cluster_num as usize];
//     for &a in asn.iter() {
//         count[a as usize] += 1;
//     }
//     let mean = asn.len() as f64 / cluster_num as f64;
//     count.into_iter().map(|c| poisson_lk(c, mean)).sum()
// }

// Return the gain of likelihood for a given dataset.
// To calculate the gain,we compute the following metrics:
// 1. Sum up the vectors for each cluster.
// 2. For each sum, the element of the positve values would be selected,
// or we acceept the edit operation at that point.
// 3. Sum up the positive values of summing-upped vector for each cluster.
// fn score(data: &[Vec<f64>], asn: &[u8], k: u8) -> f64 {
//     let dim = data[0].len();
//     let mut sums = vec![vec![0f64; dim]; k as usize];
//     for (xs, &asn) in data.iter().zip(asn.iter()) {
//         xs.iter()
//             .zip(sums[asn as usize].iter_mut())
//             .for_each(|(x, y)| *y += x);
//     }
//     sums.iter()
//         .map(|xs| -> f64 { xs.iter().map(|&x| x.max(0f64)).sum() })
//         .sum()
// }

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

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

// #[allow(dead_code)]
// fn kmeans_f64_with_init(data: &[Vec<f64>], assignments: &mut [u8], k: u8) {
//     let mut is_updated = true;
//     while is_updated {
//         let centers: Vec<Vec<f64>> = (0..k)
//             .filter_map(|cl| {
//                 let (mut count, mut slots) = (0, vec![0f64; data[0].len()]);
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
//         is_updated = false;
//         for (x, asn) in data.iter().zip(assignments.iter_mut()) {
//             let (new_asn, _) = centers
//                 .iter()
//                 .enumerate()
//                 .map(|(i, center)| (i as u8, euclid_norm_f64(center, x)))
//                 .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
//                 .unwrap();
//             if new_asn != *asn {
//                 is_updated = true;
//                 *asn = new_asn;
//             }
//         }
//     }
// }

// #[allow(dead_code)]
// fn kmeans_f64_plusplus<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
//     let mut centers: Vec<&[f64]> = vec![];
//     let indices: Vec<_> = (0..data.len()).collect();
//     // Choosing centers.
//     use rand::seq::SliceRandom;
//     while centers.len() < k as usize {
//         // calculate distance to the most nearest centers.
//         let mut dists: Vec<_> = data
//             .iter()
//             .map(|xs| {
//                 centers
//                     .iter()
//                     .map(|c| euclid_norm_f64(xs, c))
//                     .min_by(|x, y| x.partial_cmp(y).unwrap())
//                     .unwrap_or(1f64)
//                     .powi(2)
//             })
//             .collect();
//         let total: f64 = dists.iter().sum();
//         dists.iter_mut().for_each(|x| *x /= total);
//         let idx = *indices.choose_weighted(rng, |&idx| dists[idx]).unwrap();
//         centers.push(&data[idx]);
//     }
//     let mut assignments: Vec<_> = data
//         .iter()
//         .map(|xs| {
//             centers
//                 .iter()
//                 .enumerate()
//                 .map(|(i, center)| (i as u8, euclid_norm_f64(center, xs)))
//                 .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
//                 .unwrap()
//                 .0
//         })
//         .collect();
//     kmeans_f64_with_init(data, &mut assignments, k);
//     // clustering_features(data, &mut assignments, k);
//     assignments
// }

// Take dataset, the number of the cluster, the haploid coverage and seed generator,
// return the total likelihood.
fn mcmc_clustering<R: Rng>(
    data: &[Vec<f64>],
    k: usize,
    cov: f64,
    rng: &mut R,
) -> (Vec<usize>, f64) {
    if k <= 1 || data.iter().all(|xs| xs.is_empty()) || data.len() <= k {
        return (vec![0; data.len()], 0f64);
    }
    (0..20)
        .map(|i| {
            let mut assignments: Vec<_> = match i % 2 == 0 {
                true => (0..data.len()).map(|_| rng.gen_range(0..k)).collect(),
                false => kmeans(data, k, rng),
            };
            let lk = mcmc(data, &mut assignments, k, cov, rng);
            (assignments, lk)
        })
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap()
}

fn kmeans<R: Rng>(data: &[Vec<f64>], k: usize, rng: &mut R) -> Vec<usize> {
    assert!(1 < k);
    fn update_assignments(data: &[Vec<f64>], centers: &[Vec<f64>], assignments: &mut [usize]) {
        data.iter().zip(assignments.iter_mut()).for_each(|(xs, c)| {
            let (new_center, _) = centers
                .iter()
                .map(|cs| -> f64 { cs.iter().zip(xs).map(|(x, y)| (x - y).powi(2)).sum() })
                .enumerate()
                .min_by(|x, y| (x.1.partial_cmp(&(y.1)).unwrap()))
                .unwrap_or_else(|| panic!("{:?}\t{:?}", centers, xs));
            *c = new_center;
        });
    }
    fn update_centers(
        data: &[Vec<f64>],
        centers: &mut [Vec<f64>],
        counts: &mut [usize],
        assignments: &[usize],
    ) {
        centers
            .iter_mut()
            .for_each(|cs| cs.iter_mut().for_each(|x| *x = 0f64));
        counts.iter_mut().for_each(|c| *c = 0);
        for (&asn, xs) in assignments.iter().zip(data.iter()) {
            centers[asn].iter_mut().zip(xs).for_each(|(c, x)| *c += x);
            counts[asn] += 1;
        }
        centers
            .iter_mut()
            .zip(counts.iter())
            .filter(|&(_, &cou)| 0 < cou)
            .for_each(|(cen, &cou)| cen.iter_mut().for_each(|x| *x /= cou as f64));
    }
    fn get_dist(data: &[Vec<f64>], centers: &[Vec<f64>], assignments: &[usize]) -> f64 {
        data.iter()
            .zip(assignments)
            .map(|(xs, &asn)| {
                xs.iter()
                    .zip(centers[asn].iter())
                    .map(|(x, y)| (x - y).powi(2))
                    .sum::<f64>()
            })
            .sum()
    }
    let dim = data[0].len();
    assert!(0 < dim);
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
    let (mut centers, mut counts) = (vec![vec![0f64; dim]; k], vec![0; k]);
    let mut dist = get_dist(data, &centers, &assignments);
    loop {
        // Update
        update_centers(data, &mut centers, &mut counts, &assignments);
        update_assignments(data, &centers, &mut assignments);
        let new_dist = get_dist(data, &centers, &assignments);
        if new_dist - dist < 0.00001 {
            break;
        } else {
            dist = new_dist;
        }
    }
    assignments
}

// #[allow(dead_code)]
// fn em_clustering<R: Rng>(data: &[Vec<f64>], k: usize, _cov: f64, rng: &mut R) -> (Vec<u8>, f64) {
//     if k == 1 || data.iter().all(|xs| xs.is_empty()) || data.len() <= k {
//         return (vec![0; data.len()], 0f64);
//     }
//     // 1. Construct the first assignments.
//     let mut centers: Vec<&[f64]> = vec![];
//     let indices: Vec<_> = (0..data.len()).collect();
//     // Choosing centers.
//     use rand::seq::SliceRandom;
//     while centers.len() < k as usize {
//         // calculate distance to the most nearest centers.
//         let mut dists: Vec<_> = data
//             .iter()
//             .map(|xs| {
//                 centers
//                     .iter()
//                     .map(|c| euclid_norm_f64(xs, c))
//                     .min_by(|x, y| x.partial_cmp(y).unwrap())
//                     .unwrap_or(1f64)
//                     .powi(2)
//             })
//             .collect();
//         let total: f64 = dists.iter().sum();
//         dists.iter_mut().for_each(|x| *x /= total);
//         let idx = match indices.choose_weighted(rng, |&idx| dists[idx]) {
//             Ok(res) => *res,
//             Err(why) => {
//                 for d in data.iter() {
//                     error!("{:?}", d);
//                 }
//                 panic!("{:?}", why);
//             }
//         };
//         centers.push(&data[idx]);
//     }
//     let mut weights: Vec<_> = data
//         .iter()
//         .map(|xs| {
//             let mut ws: Vec<_> = centers
//                 .iter()
//                 .map(|center| (1f64 + euclid_norm_f64(center, xs)).recip())
//                 .collect();
//             let sum: f64 = ws.iter().copied().sum();
//             for w in ws.iter_mut() {
//                 *w /= sum;
//             }
//             ws
//         })
//         .collect();
//     let lk = em_cl(data, &mut weights, k);
//     let choises: Vec<_> = (0..k).collect();
//     let assignments: Vec<_> = weights
//         .iter()
//         .map(|ws| *choises.choose_weighted(rng, |&k| ws[k]).unwrap() as u8)
//         .collect();
//     //for (ws, asn) in weights.iter().zip(assignments.iter()) {
//     //     let ws: Vec<_> = ws.iter().map(|x| format!("{:.3}", x)).collect();
//     // }
//     (assignments, lk)
// }

// Return the P(n_1,...,n_k|theta), where the model behide is Chinese Restaurant Process.
// Note that the zero-sized cluster would be supressed.
// fn log_partition(clusters: &[usize], theta: f64) -> f64 {
//     // Return log theta^K * Prod((n_k-1)!) / theta / (theta+1) / ... / (theta+n-1)
//     // where K= # non-zero elements
//     let log_theta: f64 = theta.ln();
//     fn log_fact(n: usize) -> f64 {
//         match n {
//             0 => 0f64,
//             _ => (1..n + 1).map(|x| (x as f64).ln()).sum(),
//         }
//     }
//     let numerator: f64 = clusters
//         .iter()
//         .filter(|&&x| x != 0)
//         .map(|&x| log_theta + log_fact(x - 1))
//         .sum();
//     let total: usize = clusters.iter().copied().sum();
//     let denominator: f64 = (0..total).map(|i| (theta + i as f64).ln()).sum();
//     numerator - denominator
// }

// Return the maximum likelihood.
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).
fn mcmc<R: Rng>(data: &[Vec<f64>], assign: &mut [usize], k: usize, cov: f64, rng: &mut R) -> f64 {
    let partition_lks: Vec<_> = (0..=data.len())
        .map(|x| max_poisson_lk(x, cov, 1, k))
        .collect();
    assert!(
        partition_lks.iter().all(|x| !x.is_nan()),
        "{:?}\t{}",
        partition_lks,
        cov
    );
    // how many instance in a cluster.
    let mut clusters = vec![0; k];
    // Current (un-modified) likelihoods.
    let mut lks = vec![vec![0f64; data[0].len()]; k];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            *lk += x;
        }
    }
    // Current Likelihood of the data.
    let data_lk: f64 = lks.iter().flatten().fold(0f64, |x, y| x + y.max(0f64));
    let partition_lk = clusters.iter().map(|&x| partition_lks[x]).sum::<f64>();
    let mut lk = data_lk + partition_lk;
    // MAP estimation.
    let (mut max, mut argmax) = (lk, assign.to_vec());
    // Temprature from 1 -> 0, Linear scale?
    let total = 2000 * data.len();
    for (_t, temperature) in (1..total).map(|t| (t, (t as f64).recip())) {
        let params = (partition_lks.as_slice(), temperature);
        let (is_success, diff) = update(data, assign, k, &mut lks, &mut clusters, params, rng);
        if is_success {
            lk += diff;
            if max < lk {
                max = lk;
                argmax
                    .iter_mut()
                    .zip(assign.iter())
                    .for_each(|(x, &y)| *x = y);
            }
        }
    }
    // get max value.
    clusters.iter_mut().for_each(|c| *c = 0);
    lks.iter_mut().flatten().for_each(|x| *x = 0f64);
    assign.iter_mut().zip(argmax).for_each(|(x, y)| *x = y);
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            *lk += x;
        }
    }
    let data_lk = lks.iter().flatten().fold(0f64, |x, y| x + y.max(0f64));
    let partition_lk = clusters.iter().map(|&x| partition_lks[x]).sum::<f64>();
    data_lk + partition_lk
}

fn update<R: Rng>(
    data: &[Vec<f64>],
    assign: &mut [usize],
    k: usize,
    lks: &mut [Vec<f64>],
    clusters: &mut [usize],
    (partition_lks, temperature): (&[f64], f64),
    rng: &mut R,
) -> (bool, f64) {
    let (idx, new) = (rng.gen_range(0..data.len()), rng.gen_range(0..k));
    let (old, xs) = (assign[idx], &data[idx]);
    if old == new {
        return (false, 0f64);
    }
    let prev_data_lk = lks[old]
        .iter()
        .chain(lks[new].iter())
        .fold(0f64, |x, y| x + y.max(0f64));
    let prev_part_lk = partition_lks[clusters[old]] + partition_lks[clusters[new]];
    let prev_lk = prev_data_lk + prev_part_lk;
    // Change the assignment.
    assign[idx] = new;
    clusters[old] -= 1;
    clusters[new] += 1;
    for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
        *lk -= x;
    }
    for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
        *lk += x;
    }
    let prop_data_lk = lks[old]
        .iter()
        .chain(lks[new].iter())
        .fold(0f64, |x, y| x + y.max(0f64));
    let prop_part_lk = partition_lks[clusters[old]] + partition_lks[clusters[new]];
    let prop_lk = prop_data_lk + prop_part_lk;
    let diff = prop_lk - prev_lk;
    // Calc likelihoods.
    let prob = (diff / temperature).exp();
    if 0f64 < diff || rng.gen_bool(prob) {
        (true, diff)
    } else {
        assign[idx] = old;
        clusters[old] += 1;
        clusters[new] -= 1;
        for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
            *lk += x;
        }
        for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
            *lk -= x;
        }
        (false, diff)
    }
}

// Return the maximum likelihood.
// Traditional EM clustering.
// pub fn em_cl<R: Rng>(data: &[Vec<f64>], k: usize, rng: &mut R) -> (Vec<usize>, f64) {
//     fn model_fraction(
//         data: &[Vec<f64>],
//         weights: &[Vec<f64>],
//         k: usize,
//     ) -> (Vec<Vec<bool>>, Vec<f64>) {
//         let mut fractions = vec![0f64; k];
//         for ws in weights.iter() {
//             for (c, w) in fractions.iter_mut().zip(ws.iter()) {
//                 *c += w;
//             }
//         }
//         for w in fractions.iter_mut() {
//             *w /= data.len() as f64;
//         }
//         // Current (un-modified) likelihoods.
//         let mut lks = vec![vec![0f64; data[0].len()]; k];
//         for (xs, ws) in data.iter().zip(weights.iter()) {
//             for (w, lks) in ws.iter().zip(lks.iter_mut()) {
//                 for (lk, x) in lks.iter_mut().zip(xs.iter()) {
//                     *lk += x * w;
//                 }
//             }
//         }
//         let models: Vec<Vec<bool>> = lks
//             .iter()
//             .map(|lks| lks.iter().map(|lk| lk.is_sign_positive()).collect())
//             .collect();
//         (models, fractions)
//     }
//     fn get_lk(data: &[Vec<f64>], models: &[Vec<bool>], fractions: &[f64]) -> f64 {
//         data.iter()
//             .map(|xs| -> f64 {
//                 let xs: Vec<_> = models
//                     .iter()
//                     .map(|ms| -> f64 {
//                         xs.iter()
//                             .zip(ms.iter())
//                             .filter_map(|(x, &b)| b.then(|| x))
//                             .sum()
//                     })
//                     .zip(fractions.iter())
//                     .map(|(lk, w)| lk + w.ln())
//                     .collect();
//                 logsumexp(&xs)
//             })
//             .sum()
//     }
//     fn update_weight(weights: &mut [f64], xs: &[f64], models: &[Vec<bool>], fractions: &[f64]) {
//         for ((m, f), w) in models.iter().zip(fractions.iter()).zip(weights.iter_mut()) {
//             *w = xs.iter().zip(m).filter_map(|(&x, &b)| b.then(|| x)).sum();
//             *w += f.ln();
//         }
//         let total = logsumexp(weights);
//         for w in weights.iter_mut() {
//             *w = (*w - total).exp();
//         }
//     }
//     let mut weights: Vec<Vec<f64>> = data
//         .iter()
//         .map(|_| {
//             let mut bucket = vec![0; k];
//             let len = 100;
//             for _ in 0..len {
//                 bucket[rng.gen_range(0..k)] += 1;
//             }
//             bucket.iter().map(|&x| x as f64 / 100.0).collect()
//         })
//         .collect();
//     // how many instance in a cluster.
//     let (mut models, mut fractions) = model_fraction(data, &weights, k);
//     let mut lk = get_lk(data, &models, &fractions);
//     loop {
//         for (ws, xs) in weights.iter_mut().zip(data.iter()) {
//             update_weight(ws, xs, &models, &fractions);
//         }
//         let (new_models, new_fractions) = model_fraction(data, &weights, k);
//         let new_lk = get_lk(data, &new_models, &new_fractions);
//         if lk + 0.00001 < new_lk {
//             lk = new_lk;
//             models = new_models;
//             fractions = new_fractions;
//         } else {
//             break;
//         }
//     }
//     let assignments: Vec<usize> = weights
//         .iter()
//         .map(|ws| {
//             ws.iter()
//                 .enumerate()
//                 .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
//                 .map(|x| x.0)
//                 .unwrap()
//         })
//         .collect();
//     (assignments, lk)
// }

// fn kmeans_f64<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
//     let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
//     kmeans_f64_with_init(data, &mut assignments, k);
//     assignments
// }

// Return the distance between xs and ys.
// fn euclid_norm_f64(xs: &[f64], ys: &[f64]) -> f64 {
//     assert_eq!(ys.len(), xs.len());
//     xs.iter()
//         .zip(ys.iter())
//         .map(|(x, &y)| (x - y).powi(2))
//         .sum::<f64>()
//         .sqrt()
// }

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
