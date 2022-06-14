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
/// Clustering given sequences. Return assignments, template, and LK.
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
    let strands = vec![true; reads.len()];
    let result = clustering_dev(&template, reads, &mut ops, &strands, rng, &hmm, config);
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
    strands: &[bool],
    rng: &mut R,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
    config: &ClusteringConfig,
) -> Option<ClusteringDevResult> {
    //trace!("{}", String::from_utf8_lossy(template));
    // for (i, &is_forward) in strands.iter().enumerate() {
    //     trace!("VARS\t{i}\t-1\t{}", is_forward as u8);
    // }
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
    let selected_variants: Vec<Vec<_>> = {
        let probes = filter_profiles(&profiles, strands, copy_num, coverage, average_lk, rng);
        for &(pos, _) in probes.iter() {
            for (i, (prof, is_forward)) in profiles.iter().zip(strands.iter()).enumerate() {
                let lk = prof[pos];
                trace!("STRAND\t{pos}\t{i}\t{lk:.3}\t{is_forward}");
            }
            let lks = profiles.iter().map(|p| p[pos]);
            let paired = lks.zip(strands.iter().copied());
            let (diff, thr) = is_explainable_by_strandedness(paired, rng, 0.01, pos);
            trace!("STATTEST\t{pos}\t{diff}\t{thr}");
        }
        for &(pos, lk) in probes.iter() {
            let counts = profiles.iter().filter(|xs| 0f64 < xs[pos]).count();
            use kiley::hmm::guided::NUM_ROW;
            let (pos, var) = (pos / NUM_ROW, pos % NUM_ROW);
            trace!("DUMP\t{pos}\t{var}\t{lk}\t{counts}");
        }
        profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect()
    };
    clustering_selected_variants(&selected_variants, copy_num, coverage, average_lk, rng)
}

// Make this function faster?
fn is_explainable_by_strandedness<I, R: Rng>(
    profiles: I,
    rng: &mut R,
    fprate: f64,
    pos: usize,
) -> (f64, f64)
where
    I: std::iter::Iterator<Item = (f64, bool)>,
{
    let (mut fcount, mut fsum, mut bcount, mut bsum) = (0, 0f64, 0, 0f64);
    let mut lks = Vec::with_capacity(50);
    for (lkdiff, strand) in profiles {
        lks.push(lkdiff);
        if strand {
            fsum += lkdiff;
            fcount += 1;
        } else {
            bsum += lkdiff;
            bcount += 1;
        }
    }
    // If fdiv == 0 or bdiv == 0, then all the mean_diff would be the same,
    // and diff < mean_diff[thr] + 0.01 is always hold.
    if fcount == 0 || bcount == 0 {
        return (0f64, 0.01);
    }
    let diff = (fsum / fcount as f64 - bsum / bcount as f64).abs();
    let (fcount, bcount) = (fcount.min(bcount), fcount.max(bcount));
    use rand::seq::SliceRandom;
    const SAMPLE_NUM: usize = 1500;
    let total: f64 = lks.iter().sum();
    let mut mean_diffs: Vec<_> = (0..SAMPLE_NUM)
        .map(|_| {
            let (frec, _) = lks.partial_shuffle(rng, fcount);
            let fsum: f64 = frec.iter().sum();
            let bsum = total - fsum;
            (fsum / fcount as f64 - bsum / bcount as f64).abs()
        })
        .collect();
    mean_diffs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    mean_diffs.reverse();
    let thr = (SAMPLE_NUM as f64 * fprate).ceil() as usize;
    for x in mean_diffs.iter() {
        trace!("SHUFFLE\t{pos}\t{x}");
    }
    (diff, mean_diffs[thr] + 0.01)
}

fn clustering_selected_variants<R: Rng>(
    selected_variants: &[Vec<f64>],
    copy_num: usize,
    coverage: f64,
    average_lk: f64,
    rng: &mut R,
) -> Option<ClusteringDevResult> {
    let init_copy_num = copy_num.max(5) - 3;
    let datasize = selected_variants.len();
    let coverage_imp_thr = (coverage / 2f64).ceil() as usize;
    let (mut assignments, mut max) = (vec![0; datasize], 0f64);
    let mut read_lk_gains = vec![0f64; datasize];
    for k in init_copy_num..=copy_num {
        let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
        let new_read_lk_gains = get_read_lk_gains(&selected_variants, &asn, k);
        let improved_reads = read_lk_gains
            .iter()
            .zip(new_read_lk_gains.iter())
            .filter(|(&x, &y)| x + 0.1 < y)
            .count();
        let expected_gain = average_lk * improved_reads as f64;
        trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t0");
        if expected_gain < score - max && coverage_imp_thr < improved_reads {
            assignments = asn;
            max = score;
            read_lk_gains = new_read_lk_gains;
        } else {
            break;
        }
    }
    if log_enabled!(log::Level::Trace) {
        for (i, asn) in assignments.iter().enumerate() {
            trace!("VARS\t{i}\t-2\t{asn}");
        }
        for (i, prf) in selected_variants.iter().enumerate() {
            for (pos, x) in prf.iter().enumerate() {
                trace!("VARS\t{i}\t{pos}\t{x}");
            }
        }
        // trace!("VARS\tCOPYNUM\tID\tASN\tPOS\tLKDiff\tFilter");
        // let asn_iter = assignments.iter().zip(selected_variants.iter());
        // let mut id = 0;
        // for c in 0..copy_num {
        //     for (i, prf) in asn_iter.clone().filter(|x| *x.0 == c) {
        //         for (pos, x) in prf.iter().enumerate() {
        //             trace!("VARS\t{copy_num}\t{id}\t{i}\t{pos}\t{x}\tBefore");
        //         }
        //         id += 1;
        //     }
        // }
    }
    let selected_variants = filter_suspicious_variants(&selected_variants, &assignments);
    let (mut assignments, mut max, mut max_k) = (vec![0; datasize], 0f64, 1);
    let mut read_lk_gains = vec![0f64; datasize];
    for k in init_copy_num..=copy_num {
        let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
        let new_read_lk_gains = get_read_lk_gains(&selected_variants, &asn, k);
        let improved_reads = read_lk_gains
            .iter()
            .zip(new_read_lk_gains.iter())
            .filter(|(&x, &y)| x + 0.1 < y)
            .count();
        let expected_gain = average_lk * improved_reads as f64;
        trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t1");
        if expected_gain < score - max && coverage_imp_thr < improved_reads {
            assignments = asn;
            max = score;
            max_k = k;
            read_lk_gains = new_read_lk_gains;
        } else {
            break;
        }
    }
    let mut likelihood_gains = get_likelihood_gain(&selected_variants, &assignments);
    to_posterior_probability(&mut likelihood_gains);
    // if log_enabled!(log::Level::Trace) {
    //     let asn_iter = assignments.iter().zip(selected_variants.iter());
    //     let mut id = 0;
    //     for c in 0..copy_num {
    //         for (i, prf) in asn_iter.clone().filter(|x| *x.0 == c) {
    //             for (pos, x) in prf.iter().enumerate() {
    //                 trace!("VARS\t{copy_num}\t{id}\t{i}\t{pos}\t{x}\tAfter");
    //             }
    //             id += 1;
    //         }
    //     }
    // }
    likelihood_gains
        .iter()
        .for_each(|x| assert_eq!(x.len(), max_k));
    Some((assignments, likelihood_gains, max, max_k))
}

// LK->LK-logsumexp(LK).
fn to_posterior_probability(lks: &mut [Vec<f64>]) {
    for xs in lks.iter_mut() {
        let total = logsumexp(xs);
        xs.iter_mut().for_each(|x| *x -= total);
    }
}

#[allow(dead_code)]
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

fn get_read_lk_gains(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> Vec<f64> {
    let var_num = variants[0].len();
    let is_used_position: Vec<Vec<_>> = (0..k)
        .map(|cl| {
            (0..var_num)
                .map(|pos| {
                    let lk_sum: f64 = variants
                        .iter()
                        .zip(assignments.iter())
                        .filter_map(|(vars, &asn)| (asn == cl).then(|| vars[pos]))
                        .sum();
                    0f64 < lk_sum
                })
                .collect()
        })
        .collect();
    variants
        .iter()
        .zip(assignments.iter())
        .map(|(vars, &asn)| -> f64 {
            vars.iter()
                .zip(is_used_position[asn].iter())
                .filter_map(|(&v, &b)| b.then(|| v))
                .sum()
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

// ROUND * cluster num variants would be selected.
const ROUND: usize = 3;
const MIN_COV_FRACTION: usize = 10;
fn filter_profiles<T: std::borrow::Borrow<[f64]>, R: Rng>(
    profiles: &[T],
    strands: &[bool],
    cluster_num: usize,
    coverage: f64,
    average_lk: f64,
    rng: &mut R,
) -> Vec<(usize, f64)> {
    // (sum, maximum gain, number of positive element)
    let mut total_improvement = vec![(0f64, 0); profiles[0].borrow().len()];
    for prof in profiles.iter().map(|x| x.borrow()) {
        for ((maxgain, count), p) in total_improvement.iter_mut().zip(prof) {
            *maxgain += p.max(0f64);
            *count += (0.00001 < *p) as usize;
        }
    }
    let template_len = profiles[0].borrow().len() / kiley::hmm::guided::NUM_ROW;
    let base_position = |pos: usize| pos / kiley::hmm::guided::NUM_ROW;
    let in_mask = |pos: usize| {
        let pos = base_position(pos);
        pos < MASK_LENGTH || (template_len - MASK_LENGTH) < pos
    };
    let probes: Vec<(usize, f64)> = total_improvement
        .into_iter()
        .enumerate()
        .filter(|&(pos, _)| !in_mask(pos))
        .filter(|&(_, (_, count))| coverage.ceil() as usize / MIN_COV_FRACTION < count)
        .filter(|&(_, (gain, count))| count as f64 * average_lk < gain)
        .filter(|&(pos, _)| {
            let lks = profiles.iter().map(|p| p.borrow()[pos]);
            let paired = lks.zip(strands.iter().copied());
            let (diff, thr) = is_explainable_by_strandedness(paired, rng, 0.01, pos);
            diff < thr
        })
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
    'outer: for _ in 0..ROUND {
        let mut weights: Vec<_> = vec![1f64; probes.len()];
        // for _ in 0..cluster_num.max(2) - 1 {
        for _ in 0..cluster_num.max(2) {
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

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

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
    // for (_t, temperature) in (1..total).map(|t| (t, (t as f64).recip())) {
    for _ in 1..total {
        let temperature = 1f64;
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
    if 0f64 < diff || rng.gen_bool((diff / temperature).exp()) {
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
