//! A small K-means clustering algorithm.
// First and last `MASK_LENGTH` bases would not be considered in variant calling.
// Should be greater than the maximum length of kiley::hmm::guided::COPY_SIZE or DEL_SIZE.
const MASK_LENGTH: usize = 7;
const MAX_HOMOP_LENGTH: usize = 2;
use crate::likelihood_gains::{Gains, Pvalues};
use kiley::hmm::guided::NUM_ROW;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig<'a> {
    band_width: usize,
    gains: &'a Gains,
    // Coverage for haploid.
    coverage: f64,
    pub copy_num: u8,
}

impl<'a> ClusteringConfig<'a> {
    pub fn new(band_width: usize, copy_num: u8, coverage: f64, gains: &'a Gains) -> Self {
        Self {
            band_width,
            coverage,
            gains,
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
    cluster_num: usize,
    band: usize,
) -> Option<ClusteringResult> {
    let mut template = kiley::ternary_consensus_by_chunk(reads, band);
    let mut hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    hmm.fit_naive(&template, reads, band);
    let gains = crate::likelihood_gains::estimate_gain(&hmm, 4283094, 100, 20, 5);
    let cov = reads.len() as f64 / cluster_num as f64;
    let config = ClusteringConfig::new(band, cluster_num as u8, cov, &gains);
    let mut ops: Vec<_> = reads
        .iter()
        .map(|x| hmm.align(&template, x.borrow(), band).1)
        .collect();
    for t in 1..3 {
        hmm.fit_naive_with(&template, reads, &ops, band / t);
        template = hmm.polish_until_converge_with(&template, reads, &mut ops, band / t);
    }
    let strands = vec![true; reads.len()];
    let result = clustering_dev(&template, reads, &mut ops, &strands, rng, &hmm, &config);
    result.map(|(asn, gains, lk, _)| (asn, gains, lk, template))
}

fn modification_table<T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &mut [Vec<kiley::Op>],
    band: usize,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
) -> Vec<Vec<f64>> {
    reads
        .iter()
        .zip(ops.iter_mut())
        .map(|(seq, op)| {
            let (mut table, lk) = hmm.modification_table(template, seq.borrow(), band, op);
            table.iter_mut().for_each(|x| *x -= lk);
            table
        })
        .collect()
}

fn filter_by(profiles: &[Vec<f64>], probes: &[(usize, f64)]) -> Vec<Vec<f64>> {
    profiles
        .iter()
        .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
        .collect()
}

type ClusteringDevResult = (Vec<usize>, Vec<Vec<f64>>, f64, usize);
pub fn clustering_dev<R: Rng, T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &mut [Vec<kiley::Op>],
    strands: &[bool],
    rng: &mut R,
    hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
    config: &ClusteringConfig,
) -> Option<ClusteringDevResult> {
    let profiles = modification_table(template, reads, ops, config.band_width, hmm);
    let probes = filter_profiles(template, &profiles, strands, config, rng);
    let op_and_homop = operation_and_homopolymer_length(template, &probes);
    let variants = filter_by(&profiles, &probes);
    let feature_vectors = (variants.as_slice(), op_and_homop.as_slice());
    let clustering_result = cluster_filtered_variants(feature_vectors, config, rng);
    let (mut assignments, mut likelihood_gains, max, max_k) = clustering_result;
    if log_enabled!(log::Level::Trace) {
        for (i, (pos, lk)) in probes.iter().enumerate() {
            let sum: f64 = profiles.iter().map(|prof| prof[*pos].max(0f64)).sum();
            let (pos, ed) = (pos / NUM_ROW, pos % NUM_ROW);
            trace!("DUMP\t{i}\t{pos}\t{ed}\t{lk:.1}\t{sum:.1}");
        }
        for (i, prof) in variants.iter().enumerate() {
            for (idx, x) in prof.iter().enumerate() {
                trace!("VARS\t{i}\t{idx}\t{x}");
            }
        }
        for (i, &strand) in strands.iter().enumerate() {
            trace!("VARS\t{i}\t-1\t{}", strand as usize)
        }
        for (i, asn) in assignments.iter().enumerate() {
            trace!("VARS\t{i}\t-2\t{asn}");
        }
    }
    // Tune the assignments, is it ok .... ?
    for (lks, asn) in likelihood_gains.iter().zip(assignments.iter_mut()) {
        let iter = lks.iter().enumerate();
        let (i, max) = iter.max_by(|x, y| x.1.partial_cmp(y.1).unwrap()).unwrap();
        if lks[*asn] + 0.001 < *max {
            *asn = i;
        }
    }
    to_posterior_probability(&mut likelihood_gains);
    Some((assignments, likelihood_gains, max, max_k as usize))
}

use crate::likelihood_gains::DiffType;
fn pos_to_bp_and_difftype(pos: usize) -> (usize, DiffType) {
    let (bp, op) = (pos / NUM_ROW, pos % NUM_ROW);
    let diff_type = if op < 4 {
        DiffType::Subst
    } else if op < 8 + kiley::hmm::guided::COPY_SIZE {
        DiffType::Ins
    } else {
        DiffType::Del
    };
    (bp, diff_type)
}

fn operation_and_homopolymer_length(
    template: &[u8],
    probes: &[(usize, f64)],
) -> Vec<(usize, crate::likelihood_gains::DiffType)> {
    let homop_length = homopolymer_length(template);
    probes
        .iter()
        .map(|&(pos, _)| {
            let (bp_pos, diff_type) = pos_to_bp_and_difftype(pos);
            let homop_len = *homop_length.get(bp_pos).unwrap_or(&0);
            (homop_len, diff_type)
        })
        .collect()
}

fn homopolymer_length(xs: &[u8]) -> Vec<usize> {
    let mut current = xs[0];
    let mut length = 0;
    let mut homop = vec![];
    for &x in xs.iter() {
        if x == current {
            length += 1;
        } else {
            homop.extend(std::iter::repeat(length).take(length));
            current = x;
            length = 1;
        }
    }
    homop.extend(std::iter::repeat(length).take(length));
    assert_eq!(homop.len(), xs.len());
    homop
}

fn cluster_filtered_variants<R: Rng>(
    (variants, variant_type): (&[Vec<f64>], &[(usize, DiffType)]),
    config: &ClusteringConfig,
    rng: &mut R,
) -> ClusteringDevResult {
    let copy_num = config.copy_num as usize;
    let coverage = config.coverage;
    let gains = config.gains;
    if copy_num <= 1 || variants.iter().all(|xs| xs.is_empty()) || variants.len() <= copy_num {
        let asn = vec![0; variants.len()];
        let lk_gains = vec![vec![0f64; 1]; variants.len()];
        return (asn, lk_gains, 0f64, 1);
    }
    let coverage_imp_thr = 0;
    //    let coverage_imp_thr = get_read_thr(variants.len() as f64 / copy_num as f64, 3f64);
    let datasize = variants.len();
    let (mut assignments, mut max, mut max_k, mut read_lk_gains) =
        (vec![0; datasize], 0f64, 1, vec![0f64; datasize]);
    let mut prev_used_columns = vec![false; variants[0].len()];
    let init_copy_num = copy_num.max(4) - 2;
    for k in init_copy_num..=copy_num {
        let (asn, score, new_lk_gains, used_columns) = match k == 2 {
            false => mcmc_clustering(variants, k, coverage, rng),
            true => {
                let mcmc_gains = mcmc_clustering(variants, k, coverage, rng);
                let highest_gain = use_highest_gain(variants);
                if mcmc_gains.1 < highest_gain.1 {
                    highest_gain
                } else {
                    mcmc_gains
                }
            }
        };
        trace!("LK\t{k}\t{score:.3}");
        let min_gain = min_gain(gains, variant_type, &used_columns);
        let improved_reads = count_improved_reads(&new_lk_gains, &read_lk_gains, min_gain);
        let expected_gain_per_read =
            expected_gains(gains, variant_type, &prev_used_columns, &used_columns);
        let expected_gain = expected_gain_per_read * datasize as f64 / copy_num as f64 + 0.1;
        trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}\t{coverage_imp_thr}");
        if expected_gain < score - max && coverage_imp_thr < improved_reads {
            let mut counts = vec![0; k];
            for x in asn.iter() {
                counts[*x] += 1;
            }
            trace!("COUNTS\t{counts:?}");
            assignments = asn;
            max = score;
            max_k = k;
            read_lk_gains = new_lk_gains;
            prev_used_columns = used_columns;
        } else {
            break;
        }
    }
    let likelihood_gains = get_likelihood_gain(variants, &assignments, max_k);
    (assignments, likelihood_gains, max, max_k as usize)
}

// fn get_lk_of_coverage(asn: &[usize], cov: f64, k: usize) -> f64 {
//     let mut cluster = vec![0; k];
//     for &asn in asn.iter() {
//         cluster[asn] += 1;
//     }
//     cluster.iter().map(|&n| max_poisson_lk(n, cov, 1, k)).sum()
// }

// fn find_substitution(
//     used_columns: &[bool],
//     prev_columns: &[bool],
//     variant_type: &[(usize, DiffType)],
// ) -> Option<usize> {
//     std::iter::zip(used_columns, prev_columns)
//         .zip(variant_type)
//         .find(|((&now, &prev), &(_, t))| !prev && now && t == DiffType::Subst)
//         .map(|(_, &(pos, _))| pos)
// }

fn min_gain(gains: &Gains, variant_type: &[(usize, DiffType)], used_columns: &[bool]) -> f64 {
    std::iter::zip(variant_type, used_columns)
        .filter_map(|(&(homop_len, diff_type), is_used)| match is_used {
            true => Some(gains.expected(homop_len, diff_type) / 3f64),
            false => None,
        })
        .min_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(1f64)
}

const EXPT_GAIN_FACTOR: f64 = 0.6;
fn expected_gains(
    gains: &Gains,
    variant_type: &[(usize, DiffType)],
    prev_columns: &[bool],
    used_columns: &[bool],
) -> f64 {
    assert_eq!(variant_type.len(), used_columns.len());
    // Previously not used, currently used.
    let newly_used = std::iter::zip(prev_columns, used_columns).map(|(&p, &c)| !p & c);
    let expt_sum: f64 = std::iter::zip(variant_type, newly_used)
        .map(|(&(homop, diff_type), is_used)| match is_used {
            true => gains.expected(homop, diff_type),
            false => 0.00001,
        })
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(0f64);
    //        .sum();
    trace!("EXPTGAIN\t{expt_sum:.4}\t{used_columns:?}");
    (EXPT_GAIN_FACTOR * expt_sum).max(0.1)
}

fn count_improved_reads(new_gains: &[f64], old_gains: &[f64], min_gain: f64) -> usize {
    std::iter::zip(new_gains.iter(), old_gains.iter())
        .filter(|(&new, &old)| old + min_gain < new)
        .count()
}

// fn get_read_thr(cov: f64, sigma: f64) -> usize {
//     (cov - sigma * cov.sqrt()).floor() as usize
// }

fn is_explainable_by_strandedness<I, R: Rng>(
    profiles: I,
    rng: &mut R,
    fprate: f64,
    _pos: usize,
) -> bool
where
    I: std::iter::Iterator<Item = (f64, bool)>,
{
    let (mut fcount, mut fsum, mut bcount, mut bsum) = (0, 0f64, 0, 0f64);
    let mut lks = Vec::with_capacity(50);
    for (lkdiff, strand) in profiles {
        lks.push(lkdiff);
        if strand {
            fsum += lkdiff.max(0f64);
            fcount += 1;
        } else {
            bsum += lkdiff.max(0f64);
            bcount += 1;
        }
    }
    if fcount == 0 || bcount == 0 {
        return true;
    }
    let diff = fsum - bsum;
    const SAMPLE_NUM: usize = 3000;
    let mut null_diffs: Vec<_> = (0..SAMPLE_NUM)
        .map(|_| {
            lks.shuffle(rng);
            let fsum: f64 = lks[..fcount].iter().map(|x| x.max(0f64)).sum();
            let bsum: f64 = lks[fcount..].iter().map(|x| x.max(0f64)).sum();
            fsum - bsum
        })
        .collect();
    null_diffs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let lower_thr_idx = (SAMPLE_NUM as f64 * fprate).ceil() as usize;
    let upper_thr_idx = (SAMPLE_NUM as f64 * (1f64 - fprate)).ceil() as usize;
    let (lower_thr, upper_thr) = (null_diffs[lower_thr_idx], null_diffs[upper_thr_idx]);
    if log_enabled!(log::Level::Trace) {
        let (pos, ed) = pos_to_bp_and_difftype(_pos);
        trace!("STRAND\t{pos}\t{ed}\t{fsum:.2}\t{fcount}\t{bsum:.2}\t{bcount}\t{lower_thr:.2}\t{upper_thr:.2}");
    }
    lower_thr < diff && diff < upper_thr
}

// LK->LK-logsumexp(LK).
fn to_posterior_probability(lks: &mut [Vec<f64>]) {
    for xs in lks.iter_mut() {
        let total = crate::misc::logsumexp(xs);
        xs.iter_mut().for_each(|x| *x -= total);
    }
}

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
                        .filter(|&((_, (sum, _)), to_use)| to_use & sum.is_sign_positive())
                        .map(|x| (x.0).0)
                        .sum()
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

fn get_read_lk_gains(
    variants: &[Vec<f64>],
    assignments: &[usize],
    k: usize,
) -> (Vec<bool>, Vec<f64>) {
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
    let gain_on_read: Vec<_> = variants
        .iter()
        .zip(assignments.iter())
        .map(|(vars, &asn)| {
            vars.iter()
                .zip(lks[asn].iter())
                .zip(use_columns.iter())
                .filter(|&((_, (sum, _)), to_use)| to_use & sum.is_sign_positive())
                .map(|x| (x.0).0)
                .sum()
        })
        .collect();
    (use_columns, gain_on_read)
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

// ROUND * cluster num variants would be selected.
const ROUND: usize = 3;
const PVALUE: f64 = 0.05;
// False positive rate to determine the strand bias. In other words,
// The probability that the variant is regarded as biased even if it is not
// is 0.05. It essentially sacrifice 5% variants under the name of the strand bias.
const FP_RATE: f64 = 0.05;
const MIN_REQ_FRACTION: f64 = 0.5;
fn filter_profiles<T: std::borrow::Borrow<[f64]>, R: Rng>(
    template: &[u8],
    profiles: &[T],
    strands: &[bool],
    config: &ClusteringConfig,
    rng: &mut R,
) -> Vec<(usize, f64)> {
    let cluster_num = config.copy_num as usize;
    let coverage = config.coverage;
    let gains = config.gains;
    trace!("\n{gains}");
    let pvalues = gains.pvalues(profiles.len());
    let homopolymer_length = homopolymer_length(template);
    let min_req: Vec<_> = (0..profiles[0].borrow().len())
        .map(|pos| {
            let (bp, diff_type) = pos_to_bp_and_difftype(pos);
            let homop_len = *homopolymer_length.get(bp).unwrap_or(&1);
            gains.expected(homop_len, diff_type) * MIN_REQ_FRACTION
        })
        .collect();
    let total_improvement = column_sum_of(profiles, &min_req);
    let temp_len = total_improvement.len() / NUM_ROW;
    let probes: Vec<(usize, f64)> = total_improvement
        .iter()
        .enumerate()
        .filter(|&(pos, _)| {
            let (pos, _) = pos_to_bp_and_difftype(pos);
            MASK_LENGTH <= pos && pos <= (temp_len - MASK_LENGTH)
        })
        .filter(|&(pos, _)| pos % NUM_ROW < 8 || pos % NUM_ROW == 8 + kiley::hmm::guided::COPY_SIZE) // 8 -> 1 length copy = same as insertion
        .filter(|&(_, &(gain, _))| 0f64 < gain)
        .filter(|&(pos, _)| is_in_short_homopolymer(pos, &homopolymer_length, template))
        .filter(|&(pos, &improve)| {
            has_small_pvalue(pos, improve, &homopolymer_length, &pvalues, gains, temp_len)
        })
        .map(|(pos, &(maxgain, count))| {
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
            is_explainable_by_strandedness(paired, rng, FP_RATE, pos)
        })
        .collect();
    trace!("TOTAL\t{}", probes.len());
    for &(pos, lk) in probes.iter() {
        let count = total_improvement[pos].1;
        let (pos, ed) = (pos / NUM_ROW, pos % NUM_ROW);
        trace!("CAND\t{pos}\t{ed}\t{lk:.1}\t{count}");
    }
    pick_filtered_profiles(&probes, profiles, cluster_num)
}

fn has_small_pvalue(
    pos: usize,
    (gain, count): (f64, usize),
    homopolymer_length: &[usize],
    pvalues: &Pvalues,
    gains: &Gains,
    template_len: usize,
) -> bool {
    let (bp_pos, diff_type) = pos_to_bp_and_difftype(pos);
    let homop_len = *homopolymer_length.get(bp_pos).unwrap_or(&0);
    let pvalue = pvalues.pvalue(homop_len, diff_type, count);
    let expt = gains.expected(homop_len, diff_type) * EXPT_GAIN_FACTOR;
    if 10 < count {
        let pvalue = template_len as f64 * pvalue;
        let (pos, ed) = pos_to_bp_and_difftype(pos);
        let homop = &homopolymer_length[pos - 1..=pos + 1];
        trace!("PVALUE\t{pos}\t{ed}\t{gain:.3}\t{count}\t{expt:.3}\t{pvalue:.3}\t{homop:?}");
    }
    (count as f64) * expt < gain && pvalue < PVALUE / template_len as f64
}

fn is_in_short_homopolymer(pos: usize, homopolymer_length: &[usize], template: &[u8]) -> bool {
    match pos_to_bp_and_difftype(pos) {
        (x, DiffType::Ins) => {
            let base = *b"ACGT".get(pos % NUM_ROW - 4).unwrap_or(&0);
            let prev_len = match 0 < x {
                true => homopolymer_length[x - 1] + (template[x - 1] == base) as usize,
                false => (template[x - 1] == base) as usize,
            };
            let next_len = match x < template.len() {
                true => homopolymer_length[x] + (template[x] == base) as usize,
                false => (template[x] == base) as usize,
            };
            prev_len <= MAX_HOMOP_LENGTH && next_len <= MAX_HOMOP_LENGTH
        }
        (x, DiffType::Del) if x < template.len() => homopolymer_length[x] <= MAX_HOMOP_LENGTH,
        _ => true,
    }
}

fn pick_filtered_profiles<T: std::borrow::Borrow<[f64]>>(
    probes: &[(usize, f64)],
    profiles: &[T],
    cluster_num: usize,
) -> Vec<(usize, f64)> {
    // 0-> not selected yet. 1-> selected. 2-> removed because of co-occurence.
    // Currently, find and sweep. So, to select one variant,
    // It takes O(L) time to find a variant, and O(L) time to sweep linked variant.
    // So, in total, O(ML) time to select M variants. It is OK I think, because
    // usually L is 2K, M is around 3-20.
    let mut is_selected = vec![0; probes.len()];
    'outer: for _ in 0..ROUND {
        let mut weights: Vec<_> = vec![1f64; probes.len()];
        for _ in 0..cluster_num.max(2) {
            let next_var_idx = find_next_variants(probes, &weights, &is_selected);
            if let Some(next_var_idx) = next_var_idx {
                let picked_pos = probes[next_var_idx].0;
                let (picked_pos_in_bp, _) = pos_to_bp_and_difftype(picked_pos);
                is_selected[next_var_idx] = 1;
                for ((&(pos, _), weight), selected) in probes
                    .iter()
                    .zip(weights.iter_mut())
                    .zip(is_selected.iter_mut())
                    .filter(|&(_, &mut selected)| selected == 0)
                {
                    let (pos_in_bp, _) = pos_to_bp_and_difftype(pos);
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
        .iter()
        .zip(is_selected)
        .filter_map(|(x, y)| (y == 1).then(|| *x))
        .collect();
    selected_variants
}

// Return the sum of the positive values of the columns.
// Each position have exactly one chance to contribute to the total sum.
// In other words, for each position, one of the largest value would be
// chosen as "varinats".
fn column_sum_of<T: std::borrow::Borrow<[f64]>>(
    profiles: &[T],
    min_req: &[f64],
) -> Vec<(f64, usize)> {
    let mut total_improvement = vec![(0.0, 0); profiles[0].borrow().len()];
    assert_eq!(total_improvement.len(), min_req.len());
    for prof in profiles.iter().map(|x| x.borrow()) {
        for ((slot, gain), req) in total_improvement.iter_mut().zip(prof).zip(min_req) {
            if req < gain {
                slot.0 += gain;
                slot.1 += 1;
            }
        }
    }
    total_improvement
}

// fn pick_one_of_the_max(xs: &[f64], min_req: &[f64], _seed: usize) -> Option<(usize, f64)> {
//     assert_eq!(xs.len(), min_req.len());
//     std::iter::zip(xs, min_req)
//         .enumerate()
//         .max_by(|(_, x), (_, y)| x.0.partial_cmp(y.0).unwrap())
//         .filter(|&(_, (x, req))| req < x)
//         .map(|(idx, (&x, _))| (idx, x))
// }

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
        .filter(|&(_, (_, &flag))| flag == 0)
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

// Take dataset, the number of the cluster, the haploid coverage and seed generator,
// return the total likelihood.
fn mcmc_clustering<R: Rng>(
    data: &[Vec<f64>],
    k: usize,
    cov: f64,
    rng: &mut R,
) -> (Vec<usize>, f64, Vec<f64>, Vec<bool>) {
    let (assignment, score) = (0..20)
        .map(|_| {
            let mut assignments = crate::misc::kmeans(data, k, rng).1;
            let lk = mcmc_with_filter(data, &mut assignments, k, cov, rng);
            // debug!("MCMC\t{lk:.2}");
            (assignments, lk)
        })
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    let (used_columns, lk_gains) = get_read_lk_gains(data, &assignment, k);
    let mut counts = vec![0; k];
    for &a in assignment.iter() {
        counts[a] += 1;
    }
    let cluster_lk: f64 = counts.iter().map(|&c| max_poisson_lk(c, cov, 1, k)).sum();
    (assignment, score - cluster_lk, lk_gains, used_columns)
}

// Take dataset and return the total likelihood when using the most powerful feature.
fn use_highest_gain(data: &[Vec<f64>]) -> (Vec<usize>, f64, Vec<f64>, Vec<bool>) {
    let dim = data[0].len();
    let mut gains = vec![0f64; dim];
    for xs in data.iter() {
        for (slot, x) in gains.iter_mut().zip(xs) {
            *slot += x.max(0f64);
        }
    }
    let (max_idx, _) = gains
        .iter()
        .enumerate()
        .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
        .unwrap();
    let assignments: Vec<_> = data
        .iter()
        .map(|xs| (0f64 < xs[max_idx]) as usize)
        .collect();
    let (used_columns, lk_gains) = get_read_lk_gains(data, &assignments, 2);
    let score: f64 = lk_gains.iter().sum();
    (assignments, score, lk_gains, used_columns)
}

// Return the maximum likelihood.
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).

use rand::prelude::SliceRandom;
use rand::seq::IteratorRandom;

// use super::HOMOP_LEN;
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
    // let init = lk;
    for _t in 0..total {
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
    let mut lk: f64 = clusters.iter().map(|&size| _size_to_lk[size]).sum();
    for lks in lks.iter() {
        for ((total, _), _) in lks.iter().zip(use_columns.iter()).filter(|x| *x.1) {
            lk += total.max(0f64);
        }
    }
    lk
}
const POS_FRAC: f64 = 0.70;
// TODO: this is art.
const IN_POS_RATIO: f64 = 2f64;
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
    #[test]
    fn homop_length_test() {
        let xs = b"ACCCCGTTTGGTT";
        let answer = vec![1, 4, 4, 4, 4, 1, 3, 3, 3, 2, 2, 2, 2];
        let homop_len = homopolymer_length(xs);
        assert_eq!(homop_len, answer);
    }
}
