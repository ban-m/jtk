//! A small K-means clustering algorithm.
// First and last `MASK_LENGTH` bases would not be considered in variant calling.
// Should be greater than the maximum length of kiley::hmm::guided::COPY_SIZE or DEL_SIZE.
const MASK_LENGTH: usize = 7;
use kiley::hmm::guided::NUM_ROW;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct ClusteringConfig {
    band_width: usize,
    // Minimum required Likelihood gain.
    gain: f64,
    // Coverage for haploid.
    coverage: f64,
    pub copy_num: u8,
}

impl ClusteringConfig {
    pub fn new(band_width: usize, copy_num: u8, coverage: f64, gain: f64) -> Self {
        Self {
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

// pub fn clustering_dev<R: Rng, T: std::borrow::Borrow<[u8]>>(
//     template: &[u8],
//     reads: &[T],
//     ops: &mut [Vec<kiley::Op>],
//     strands: &[bool],
//     rng: &mut R,
//     hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
//     config: &ClusteringConfig,
// ) -> Option<ClusteringDevResult> {
//     // trace!("{}", String::from_utf8_lossy(template));
//     for (i, &is_forward) in strands.iter().enumerate() {
//         trace!("VARS\t{i}\t-1\t{}", is_forward as u8);
//     }
//     let ClusteringConfig {
//         band_width,
//         copy_num,
//         coverage,
//         gain: average_lk,
//         ..
//     } = *config;
//     let (profiles, _): (Vec<_>, Vec<_>) = reads
//         .iter()
//         .zip(ops.iter_mut())
//         .map(|(seq, op)| {
//             let (mut table, lk) = hmm.modification_table(template, seq.borrow(), band_width, op);
//             assert!(table.iter().all(|x| x.is_finite()));
//             assert!(table.iter().all(|x| !x.is_nan()));
//             table.iter_mut().for_each(|x| *x -= lk);
//             (table, lk)
//         })
//         .unzip();
//     let copy_num = copy_num as usize;
//     let selected_variants: Vec<Vec<_>> = {
//         let probes = filter_profiles(&profiles, strands, copy_num, coverage, average_lk, rng);
//         // for &(pos, _) in probes.iter() {
//         //     for (i, (prof, is_forward)) in profiles.iter().zip(strands.iter()).enumerate() {
//         //         let lk = prof[pos];
//         //         trace!("STRAND\t{pos}\t{i}\t{lk:.3}\t{is_forward}");
//         //     }
//         //     let lks = profiles.iter().map(|p| p[pos]);
//         //     let paired = lks.zip(strands.iter().copied());
//         //     let (diff, thr) = is_explainable_by_strandedness(paired, rng, 0.01, pos);
//         //     trace!("STATTEST\t{pos}\t{diff}\t{thr}");
//         // }
//         for &(pos, lk) in probes.iter() {
//             let counts = profiles.iter().filter(|xs| 0f64 < xs[pos]).count();
//             use kiley::hmm::guided::NUM_ROW;
//             let (pos, var) = (pos / NUM_ROW, pos % NUM_ROW);
//             trace!("DUMP\t{pos}\t{var}\t{lk}\t{counts}");
//         }
//         profiles
//             .iter()
//             .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
//             .collect()
//     };
//     clustering_selected_variants(&selected_variants, copy_num, coverage, average_lk, rng)
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
    let ClusteringConfig {
        band_width,
        copy_num,
        coverage,
        gain: average_lk,
        ..
    } = *config;
    let profiles: Vec<_> = reads
        .iter()
        .zip(ops.iter_mut())
        .map(|(seq, op)| {
            let (mut table, lk) = hmm.modification_table(template, seq.borrow(), band_width, op);
            table.iter_mut().for_each(|x| *x -= lk);
            table
        })
        .collect();
    let copy_num = copy_num as usize;
    let coverage_imp_thr = get_read_thr(reads.len() as f64 / copy_num as f64, 3f64);
    let req_improve = average_lk / 2f64;
    let req_cov_filter = coverage_imp_thr.saturating_sub(2);
    let filter_params = (average_lk / 2f64, req_cov_filter, req_improve);
    let probes = filter_profiles(&profiles, strands, copy_num, coverage, filter_params, rng);
    for (i, (pos, lk)) in probes.iter().enumerate() {
        let sum: f64 = profiles.iter().map(|prof| prof[*pos].max(0f64)).sum();
        let (pos, ed) = (pos / NUM_ROW, pos % NUM_ROW);
        trace!("DUMP\t{i}\t{pos}\t{ed}\t{lk:.1}\t{sum:.1}");
    }
    let selected_variants: Vec<Vec<_>> = profiles
        .iter()
        .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
        .collect();
    for (i, prof) in selected_variants.iter().enumerate() {
        for (idx, x) in prof.iter().enumerate() {
            trace!("VARS\t{i}\t{idx}\t{x}");
        }
    }
    let datasize = reads.len();
    let (mut assignments, mut max, mut max_k, mut read_lk_gains) =
        (vec![0; datasize], 0f64, 1, vec![0f64; datasize]);
    // let mut likelihood_gains = vec![vec![0f64]; datasize];
    let init_copy_num = copy_num.max(3) - 1;
    for k in init_copy_num..=copy_num {
        let (asn, score, new_lk_gains) = mcmc_clustering(&selected_variants, k, coverage, rng);
        trace!("LK\t{k}\t{score:.3}\t1");
        let score = score + get_lk_of_coverage(&asn, coverage, k);
        let improved_reads = std::iter::zip(read_lk_gains.iter(), new_lk_gains.iter())
            .filter(|(&x, &y)| x + req_improve < y)
            .count();
        // let (score, improved_reads) = std::iter::zip(read_lk_gains.iter(), new_lk_gains.iter())
        //     .filter(|(&x, &y)| x + req_improve < y)
        //     .fold((0f64, 0), |(lk, count), x| (lk + x.1, count + 1));
        let expected_gain = average_lk * datasize as f64 / copy_num as f64;
        // let expected_gain = average_lk * improved_reads as f64;
        trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t{coverage_imp_thr}");
        if expected_gain < score - max && coverage_imp_thr < improved_reads {
            assignments = asn;
            max = score;
            max_k = k;
            read_lk_gains = new_lk_gains;
        } else {
            break;
        }
    }
    if log_enabled!(log::Level::Trace) {
        for (i, asn) in assignments.iter().enumerate() {
            trace!("VARS\t{i}\t-2\t{asn}");
        }
        // for (i, &is_forward) in strands.iter().enumerate() {
        //     trace!("VARS\t{i}\t-1\t{}", is_forward as u8);
        // }
    }
    let mut likelihood_gains = get_likelihood_gain(&selected_variants, &assignments, max_k);
    // Tune the assignments, is it ok .... ?
    for (lks, asn) in likelihood_gains.iter().zip(assignments.iter_mut()) {
        let iter = lks.iter().enumerate();
        let (i, max) = iter.max_by(|x, y| x.1.partial_cmp(&y.1).unwrap()).unwrap();
        if lks[*asn] + 0.001 < *max {
            *asn = i;
        }
    }
    to_posterior_probability(&mut likelihood_gains);
    Some((assignments, likelihood_gains, max, max_k as usize))
}

fn get_lk_of_coverage(asn: &[usize], cov: f64, k: usize) -> f64 {
    let mut cluster = vec![0; k];
    for &asn in asn.iter() {
        cluster[asn] += 1;
    }
    cluster.iter().map(|&n| max_poisson_lk(n, cov, 1, k)).sum()
}

#[allow(dead_code)]
fn get_read_thr(cov: f64, sigma: f64) -> usize {
    (cov - sigma * cov.sqrt()).floor() as usize
}

// fn get_homopolymer_length(seq: &[u8]) -> Vec<usize> {
//     if seq.is_empty() {
//         return vec![];
//     }
//     let mut length = Vec::with_capacity(seq.len());
//     let mut seq = seq.iter();
//     let (mut current, mut len) = (seq.next().unwrap(), 1);
//     for base in seq {
//         if base == current {
//             len += 1;
//         } else {
//             length.extend(std::iter::repeat(len).take(len));
//             current = base;
//             len = 1;
//         }
//     }
//     length.extend(std::iter::repeat(len).take(len));
//     length
// }

// fn to_table(asn: &[usize], k: usize) -> Vec<usize> {
//     let mut counts = vec![0; k];
//     for &asn in asn.iter() {
//         counts[asn] += 1;
//     }
//     counts
// }

// fn re_define_variants<T: std::borrow::Borrow<[f64]>>(
//     profiles: &[T],
//     asn: &[usize],
//     counts: &[usize],
//     _coverage: f64,
//     average_lk: f64,
//     k: usize,
// ) -> Vec<Vec<f64>> {
//     let is_ok_position = get_ok_position(profiles, asn, counts, average_lk, k);
//     let vars: usize = is_ok_position.iter().filter(|&&x| x).count();
//     trace!("FILTERED\t{}", vars);
//     profiles
//         .iter()
//         .map(|prof| {
//             prof.borrow()
//                 .iter()
//                 .zip(is_ok_position.iter())
//                 .filter_map(|(&x, y)| y.then(|| x))
//                 .collect()
//         })
//         .collect()
// }

// fn get_ok_position<T: std::borrow::Borrow<[f64]>>(
//     profiles: &[T],
//     asn: &[usize],
//     counts: &[usize],
//     average_lk: f64,
//     k: usize,
// ) -> Vec<bool> {
//     let len = profiles[0].borrow().len();
//     let mut sums = vec![vec![(0, 0f64); len]; k];
//     for (&cl, prof) in asn.iter().zip(profiles.iter()) {
//         for ((count, lk), &x) in sums[cl].iter_mut().zip(prof.borrow().iter()) {
//             *count += (0f64 < x) as usize;
//             *lk += x;
//         }
//     }
//     sums.iter()
//         .zip(counts.iter())
//         .map(|(sums, &count)| pick_positive_column(sums, count, average_lk))
//         .fold(vec![false; len], |mut is_ok, filter| {
//             is_ok
//                 .iter_mut()
//                 .zip(filter.iter())
//                 .for_each(|(x, y)| *x |= *y);
//             is_ok
//         })
// }

// fn pick_positive_column(sums: &[(usize, f64)], count: usize, average_lk: f64) -> Vec<bool> {
//     const MAX_VAR_NUM: usize = 15;
//     let lk_thr = average_lk * count as f64 / 2f64;
//     let ok_position: Vec<_> = sums
//         .iter()
//         .map(|&(positive_num, lk_sum)| lk_thr < lk_sum && 0.5 < positive_num as f64 / count as f64)
//         .collect();
//     if ok_position.iter().filter(|&&x| x).count() <= MAX_VAR_NUM {
//         ok_position
//     } else {
//         let mut filtered_lks: Vec<_> = sums
//             .iter()
//             .zip(ok_position.iter())
//             .enumerate()
//             .filter_map(|(i, ((_, lk), &is_ok))| is_ok.then(|| (i, lk)))
//             .collect();
//         filtered_lks.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap());
//         let mut ok_position_filtered = ok_position;
//         for (pos, _) in filtered_lks.into_iter().skip(MAX_VAR_NUM) {
//             ok_position_filtered[pos] = false;
//         }
//         ok_position_filtered
//     }
// }

// Make this function faster?
fn is_explainable_by_strandedness<I, R: Rng>(
    profiles: I,
    rng: &mut R,
    fprate: f64,
    _pos: usize,
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
    const SAMPLE_NUM: usize = 3000;
    let total: f64 = lks.iter().sum();
    let mut mean_diffs: Vec<_> = (0..SAMPLE_NUM)
        .map(|_| {
            lks.shuffle(rng);
            let fsum: f64 = lks.iter().take(fcount).sum();
            // let (frec, _) = lks.partial_shuffle(rng, fcount);
            // let fsum: f64 = frec.iter().sum();
            let bsum = total - fsum;
            (fsum / fcount as f64 - bsum / bcount as f64).abs()
        })
        .collect();
    mean_diffs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    mean_diffs.reverse();
    let thr = (SAMPLE_NUM as f64 * fprate).ceil() as usize;
    (diff, mean_diffs[thr] + 0.01)
}

// fn clustering_selected_variants<R: Rng>(
//     selected_variants: &[Vec<f64>],
//     copy_num: usize,
//     coverage: f64,
//     average_lk: f64,
//     rng: &mut R,
// ) -> Option<ClusteringDevResult> {
//     let init_copy_num = copy_num.max(5) - 3;
//     let datasize = selected_variants.len();
//     let coverage_imp_thr = (coverage / 2f64).ceil() as usize;
//     let (mut assignments, mut max) = (vec![0; datasize], 0f64);
//     let mut read_lk_gains = vec![0f64; datasize];
//     for k in init_copy_num..=copy_num {
//         let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
//         let new_read_lk_gains = get_read_lk_gains(&selected_variants, &asn, k);
//         let improved_reads = read_lk_gains
//             .iter()
//             .zip(new_read_lk_gains.iter())
//             .filter(|(&x, &y)| x + 0.1 < y)
//             .count();
//         let expected_gain = average_lk * improved_reads as f64;
//         trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t0");
//         if expected_gain < score - max && coverage_imp_thr < improved_reads {
//             assignments = asn;
//             max = score;
//             read_lk_gains = new_read_lk_gains;
//         } else {
//             break;
//         }
//     }
//     if log_enabled!(log::Level::Trace) {
//         for (i, asn) in assignments.iter().enumerate() {
//             trace!("VARS\t{i}\t-2\t{asn}");
//         }
//         for (i, prf) in selected_variants.iter().enumerate() {
//             for (pos, x) in prf.iter().enumerate() {
//                 trace!("VARS\t{i}\t{pos}\t{x}");
//             }
//         }
//     }
//     let selected_variants = filter_suspicious_variants(&selected_variants, &assignments);
//     let (mut assignments, mut max, mut max_k) = (vec![0; datasize], 0f64, 1);
//     let mut read_lk_gains = vec![0f64; datasize];
//     for k in init_copy_num..=copy_num {
//         let (asn, score) = mcmc_clustering(&selected_variants, k, coverage, rng);
//         let new_read_lk_gains = get_read_lk_gains(&selected_variants, &asn, k);
//         let improved_reads = read_lk_gains
//             .iter()
//             .zip(new_read_lk_gains.iter())
//             .filter(|(&x, &y)| x + 0.1 < y)
//             .count();
//         let expected_gain = average_lk * improved_reads as f64;
//         trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t1");
//         if expected_gain < score - max && coverage_imp_thr < improved_reads {
//             assignments = asn;
//             max = score;
//             max_k = k;
//             read_lk_gains = new_read_lk_gains;
//         } else {
//             break;
//         }
//     }
//     let mut likelihood_gains = get_likelihood_gain(&selected_variants, &assignments);
//     to_posterior_probability(&mut likelihood_gains);
//     likelihood_gains
//         .iter()
//         .for_each(|x| assert_eq!(x.len(), max_k));
//     Some((assignments, likelihood_gains, max, max_k))
// }

// LK->LK-logsumexp(LK).
fn to_posterior_probability(lks: &mut [Vec<f64>]) {
    for xs in lks.iter_mut() {
        let total = logsumexp(xs);
        xs.iter_mut().for_each(|x| *x -= total);
    }
}

// fn number_of_improved_read(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> usize {
//     let dim = variants[0].len();
//     let mut counts = vec![0; k];
//     let mut lks: Vec<_> = vec![vec![0f64; dim]; k];
//     for (&asn, xs) in assignments.iter().zip(variants.iter()) {
//         counts[asn] += 1;
//         lks[asn]
//             .iter_mut()
//             .zip(xs.iter())
//             .for_each(|(x, y)| *x += y);
//     }
//     counts
//         .iter()
//         .zip(lks)
//         .filter_map(|(count, sums)| sums.iter().any(|x| x.is_sign_positive()).then(|| count))
//         .sum()
// }

// i->k->the likelihood gain of the i-th read when clustered in the k-th cluster.
// `copy_num` is the copy number of this unit, not the *cluster number*.
// Usually, they are the same but sometimes there are exact repeats, and the number of the cluster
// would be smaller than the copy number.
fn get_likelihood_gain(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> Vec<Vec<f64>> {
    let mut lks = vec![vec![(0f64, 0); variants[0].len()]; k];
    let mut clusters = vec![0; k];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        clusters[asn] += 1;
        for (slot, x) in lks[asn].iter_mut().zip(vars.iter()) {
            slot.0 += x;
            slot.1 += x.is_sign_positive() as usize;
        }
    }
    let use_columns = get_used_columns(&lks, &clusters);
    trace!("FILTER\t{use_columns:?}");
    variants
        .iter()
        .map(|vars| {
            lks.iter()
                .map(|slots| -> f64 {
                    vars.iter()
                        .zip(slots.iter())
                        .zip(use_columns.iter())
                        .filter(|&((_, (sum, _)), to_use)| to_use & sum.is_positive())
                        .map(|x| (x.0).0)
                        .sum()
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn get_read_lk_gains(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> Vec<f64> {
    let mut lks = vec![vec![(0f64, 0); variants[0].len()]; k];
    let mut clusters = vec![0; k];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        clusters[asn] += 1;
        for (slot, x) in lks[asn].iter_mut().zip(vars.iter()) {
            slot.0 += x;
            slot.1 += x.is_sign_positive() as usize;
        }
    }
    let use_columns = get_used_columns(&lks, &clusters);
    variants
        .iter()
        .zip(assignments.iter())
        .map(|(vars, &asn)| {
            vars.iter()
                .zip(lks[asn].iter())
                .zip(use_columns.iter())
                .filter(|&((_, (sum, _)), to_use)| to_use & sum.is_positive())
                .map(|x| (x.0).0)
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
// TODO: Can we determine FRAC and IN_POS_RATIO from the data?
// These values are critical.
const POS_FRAC: f64 = 0.70;
const IN_POS_RATIO: f64 = 3f64;
#[allow(dead_code)]
fn filter_suspicious_variants(variants: &[Vec<f64>], assignments: &[usize]) -> Vec<Vec<f64>> {
    let max_asn = *assignments.iter().max().unwrap() + 1;
    let dim = variants[0].len();
    // Assignment -> variant position -> (total #, # of positive element, likelihood sum).
    // TODO: maybe x.is_sign_poitive() is not so good. We need to do
    // THR < x, where THR is mean gain / 2 or so.
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
            in_neg as f64 * IN_POS_RATIO < in_pos as f64
        })
        .collect();
    let used_pos = var_summary
        .iter()
        .fold(vec![false; dim], |mut used_pos, cs| {
            for (pos, &(total, num_pos, lksum)) in used_pos.iter_mut().zip(cs) {
                *pos |= (0f64 < lksum) & (total as f64 * POS_FRAC < num_pos as f64);
            }
            used_pos
        });
    let bits: Vec<_> = used_pos
        .iter()
        .zip(scattered_vars.iter())
        .map(|(&b, &c)| (b & c) as usize)
        .collect();
    let now: usize = bits.iter().sum();
    trace!("FILTER\t{}\t{}\t{bits:?}", dim, now);
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
const MAX_GAIN_FAC: f64 = 0.6;
fn filter_profiles<T: std::borrow::Borrow<[f64]>, R: Rng>(
    profiles: &[T],
    strands: &[bool],
    cluster_num: usize,
    coverage: f64,
    (average_lk, req_count, req_lk): (f64, usize, f64),
    rng: &mut R,
) -> Vec<(usize, f64)> {
    // (sum, maximum gain, number of positive element)
    let total_improvement = column_sum_of(profiles, req_lk);
    let template_len = total_improvement.len() / NUM_ROW;
    let base_position = |pos: usize| (pos / NUM_ROW, pos % NUM_ROW);
    let probes: Vec<(usize, f64)> = total_improvement
        .iter()
        .enumerate()
        .filter(|&(pos, _)| {
            let pos = base_position(pos).0;
            MASK_LENGTH <= pos && pos <= (template_len - MASK_LENGTH)
        })
        .filter(|&(_, &(gain, _, max_gain))| max_gain * MAX_GAIN_FAC < gain)
        .filter(|&(pos, _)| pos % NUM_ROW < 4 + 4 + kiley::hmm::guided::COPY_SIZE)
        .filter(|&(_, &(gain, _, _))| 0f64 < gain)
        .filter(|&(_, &(_, count, _))| req_count < count)
        .filter(|&(_, &(gain, count, _))| count as f64 * average_lk < gain)
        .map(|(pos, &(maxgain, count, _))| {
            let max_lk = (1..cluster_num + 1)
                .map(|k| poisson_lk(count, coverage * k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap_or_else(|| panic!("{}", cluster_num));
            let total_lk = max_lk + maxgain;
            (pos, total_lk)
        })
        .filter(|&(_, gain)| 0f64 < gain)
        .filter(|&(pos, _)| {
            let lks = profiles.iter().map(|p| p.borrow()[pos]);
            let paired = lks.zip(strands.iter().copied());
            let (diff, thr) = is_explainable_by_strandedness(paired, rng, 0.01, pos);
            diff < thr
        })
        .collect();
    for &(pos, lk) in probes.iter() {
        let count = total_improvement[pos].1;
        let (pos, ed) = (pos / NUM_ROW, pos % NUM_ROW);
        trace!("CAND\t{pos}\t{ed}\t{lk:.1}\t{count}");
    }
    // 0-> not selected yet. 1-> selected. 2-> removed because of co-occurence.
    // Currently, find and sweep. So, to select one variant,
    // It takes O(L) time to find a variant, and O(L) time to sweep linked variant.
    // So, in total, O(ML) time to select M variants. It is OK I think, because
    // usually L is 2K, M is around 3-20.
    let mut is_selected = vec![0; probes.len()];
    'outer: for _ in 0..ROUND {
        let mut weights: Vec<_> = vec![1f64; probes.len()];
        for _ in 0..cluster_num.max(2) {
            let next_var_idx = find_next_variants(&probes, &weights, &is_selected);
            if let Some(next_var_idx) = next_var_idx {
                let picked_pos = probes[next_var_idx].0;
                let picked_pos_in_bp = base_position(picked_pos).0;
                is_selected[next_var_idx] = 1;
                for ((&(pos, _), weight), selected) in probes
                    .iter()
                    .zip(weights.iter_mut())
                    .zip(is_selected.iter_mut())
                    .filter(|&(_, &mut selected)| selected == 0)
                {
                    let pos_in_bp = base_position(pos).0;
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

// Return the sum of the positive values of the columns.
// Each position have exactly one chance to contribute to the total sum.
// In other words, for each position, one of the largest value would be
// chosen as "varinats".
fn column_sum_of<T: std::borrow::Borrow<[f64]>>(
    profiles: &[T],
    req_lk: f64,
) -> Vec<(f64, usize, f64)> {
    let mut total_improvement = vec![(0.0, 0, 0f64); profiles[0].borrow().len()];
    for (i, prof) in profiles.iter().map(|x| x.borrow()).enumerate() {
        let total_per_row = total_improvement.chunks_exact_mut(NUM_ROW);
        for (row, total_row) in prof.chunks_exact(NUM_ROW).zip(total_per_row) {
            let (i, max) = pick_one_of_the_max(row, i);
            total_row[i].2 += max.max(0f64);
            if req_lk < max {
                total_row[i].0 += max.max(0f64);
                total_row[i].1 += 1;
            }
        }
    }
    total_improvement
}

fn pick_one_of_the_max(xs: &[f64], seed: usize) -> (usize, f64) {
    let thr = 0.000001;
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let count = xs.iter().filter(|&x| (x - max).abs() < thr).count();
    let pick = seed % count;
    xs.iter()
        .copied()
        .enumerate()
        .filter(|&(_, x)| (x - max).abs() < thr)
        .nth(pick)
        .unwrap()
}

fn find_next_variants(
    probes: &[(usize, f64)],
    weights: &[f64],
    is_selected: &[u8],
) -> Option<usize> {
    probes
        .iter()
        .map(|x| x.1)
        .zip(weights.iter())
        .zip(is_selected.iter())
        .enumerate()
        .max_by(
            |(_, ((lk1, w1), picked1)), (_, ((lk2, w2), picked2))| match (picked1, picked2) {
                (1, _) | (2, _) => std::cmp::Ordering::Less,
                (_, 1) | (_, 2) => std::cmp::Ordering::Greater,
                _ => (lk1 * *w1).partial_cmp(&(lk2 * *w2)).unwrap(),
            },
        )
        .map(|x| x.0)
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
) -> (Vec<usize>, f64, Vec<f64>) {
    if k <= 1 || data.iter().all(|xs| xs.is_empty()) || data.len() <= k {
        return (vec![0; data.len()], 0f64, vec![0f64; data.len()]);
    }
    let (assignment, score) = (0..20)
        .map(|i| {
            let mut assignments: Vec<_> = match i % 2 == 0 {
                true => (0..data.len()).map(|_| rng.gen_range(0..k)).collect(),
                false => kmeans(data, k, rng),
            };
            let lk = mcmc_with_filter(data, &mut assignments, k, cov, rng);
            (assignments, lk)
        })
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    let lk_gains = get_read_lk_gains(data, &assignment, k);
    // Checking...
    let total_gain: f64 = lk_gains.iter().sum();
    assert!((score - total_gain).abs() < 0.0001, "The reg term matters?");
    (assignment, score, lk_gains)
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
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).
#[allow(dead_code)]
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
    let total = 2000 * data.len();
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

// Return the maximum likelihood.
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).

use rand::seq::IteratorRandom;
use rand_distr::num_traits::Signed;
#[allow(dead_code)]
fn mcmc_with_filter<R: Rng>(
    data: &[Vec<f64>],
    assign: &mut [usize],
    k: usize,
    cov: f64,
    rng: &mut R,
) -> f64 {
    let size_to_lk: Vec<_> = (0..=data.len())
        .map(|x| max_poisson_lk(x, cov, 1, k))
        .collect();
    let is_valid_lk = size_to_lk.iter().all(|x| !x.is_nan());
    assert!(is_valid_lk, "{:?}\t{}", size_to_lk, cov);
    // how many instance in a cluster.
    let mut clusters = vec![0; k];
    // Cluster Id -> Column ID -> (total lk, # of positive values in the cluster)
    let mut lks = vec![vec![(0f64, 0); data[0].len()]; k];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            lk.0 += x;
            lk.1 += x.is_sign_positive() as usize;
        }
    }
    let mut lk = get_lk(&lks, &clusters, &size_to_lk);
    let (mut max, mut argmax) = (lk, assign.to_vec());
    let total = 2000 * data.len();
    for _t in 1..total {
        let idx = rng.gen_range(0..data.len());
        let old = assign[idx];
        let new = (0..k).filter(|&k| k != old).choose(rng).unwrap();
        flip(data, assign, idx, new, &mut lks, &mut clusters);
        let proposed_lk = get_lk(&lks, &clusters, &size_to_lk);
        let diff = proposed_lk - lk;
        if 0f64 < diff || rng.gen_bool(diff.exp()) {
            lk = proposed_lk;
            if max < lk {
                max = proposed_lk;
                argmax
                    .iter_mut()
                    .zip(assign.iter())
                    .for_each(|(x, &y)| *x = y);
            }
        } else {
            flip(data, assign, idx, old, &mut lks, &mut clusters);
        }
    }
    // get max value.
    assign.iter_mut().zip(argmax).for_each(|(x, y)| *x = y);
    let mut clusters = vec![0; k];
    let mut lks = vec![vec![(0f64, 0); data[0].len()]; k];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            lk.0 += x;
            lk.1 += x.is_sign_positive() as usize;
        }
    }
    let lk = get_lk(&lks, &clusters, &size_to_lk);
    assert!((max - lk).abs() < 0.0001);
    max
}

fn flip(
    data: &[Vec<f64>],
    assign: &mut [usize],
    idx: usize,
    to: usize,
    lks: &mut [Vec<(f64, usize)>],
    clusters: &mut [usize],
) {
    let from = assign[idx];
    // Remove data[idx] from lks[from] and clusters[from.]
    clusters[from] -= 1;
    for (lks, x) in lks[from].iter_mut().zip(data[idx].iter()) {
        lks.0 -= x;
        lks.1 -= x.is_sign_positive() as usize;
    }
    assign[idx] = to;
    clusters[to] += 1;
    for (lks, x) in lks[to].iter_mut().zip(data[idx].iter()) {
        lks.0 += x;
        lks.1 += x.is_sign_positive() as usize;
    }
}

fn get_lk(lks: &[Vec<(f64, usize)>], clusters: &[usize], _size_to_lk: &[f64]) -> f64 {
    // If true, then that column would not be used.
    let use_columns = get_used_columns(lks, clusters);
    let mut lk = 0f64;
    //let mut lk: f64 = clusters.iter().map(|&size| size_to_lk[size]).sum();
    for lks in lks.iter() {
        for ((total, _), _) in lks.iter().zip(use_columns.iter()).filter(|x| *x.1) {
            lk += total.max(0f64);
        }
    }
    lk
}

fn get_used_columns(lks: &[Vec<(f64, usize)>], clusters: &[usize]) -> Vec<bool> {
    let mut to_uses = vec![false; lks[0].len()];
    for (lks, &size) in lks.iter().zip(clusters.iter()) {
        assert_eq!(to_uses.len(), lks.len());
        for (to_use, &(lk, pos_count)) in to_uses.iter_mut().zip(lks.iter()) {
            *to_use |= 0f64 < lk && size as f64 * POS_FRAC < pos_count as f64;
        }
    }
    for (d, to_use) in to_uses.iter_mut().enumerate() {
        let column = lks.iter().map(|lks| lks[d]);
        let pos_in_use: usize = column
            .clone()
            .filter_map(|x| (0f64 < x.0).then(|| x.1))
            .sum();
        let pos_in_neg: usize = column
            .clone()
            .filter_map(|x| (x.0 <= 0f64).then(|| x.1))
            .sum();
        *to_use &= pos_in_neg as f64 * IN_POS_RATIO < pos_in_use as f64;
    }
    to_uses
}

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
    // #[test]
    // fn cacl_homopolymer_length() {
    //     let seq = b"AC";
    //     let homop = get_homopolymer_length(seq);
    //     assert_eq!(homop, vec![1, 1]);
    //     let seq = b"ACC";
    //     let homop = get_homopolymer_length(seq);
    //     assert_eq!(homop, vec![1, 2, 2]);
    //     let seq = b"ACCA";
    //     let homop = get_homopolymer_length(seq);
    //     assert_eq!(homop, vec![1, 2, 2, 1]);
    //     let seq = b"AAACAAAGC";
    //     let homop = get_homopolymer_length(seq);
    //     assert_eq!(homop, vec![3, 3, 3, 1, 3, 3, 3, 1, 1]);
    // }
}
