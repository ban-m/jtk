//! A small K-means clustering algorithm.
// First and last `MASK_LENGTH` bases would not be considered in variant calling.
const MASK_LENGTH: usize = 7;
const MAX_HOMOP_LENGTH: usize = 2;
const POS_THR: f64 = 0.00001;
use crate::likelihood_gains::{Gains, Pvalues};
use kiley::hmm::NUM_ROW;
use log::*;
use rand::Rng;

type FeatureVector = (
    Vec<Vec<f64>>,
    Vec<(usize, crate::likelihood_gains::DiffType)>,
);
type ClusteringDevResult = (Vec<usize>, Vec<Vec<f64>>, f64, usize);

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig<'a> {
    pub band_width: usize,
    gains: &'a Gains,
    // Coverage for haploid.
    pub coverage: f64,
    pub copy_num: usize,
    pub local_coverage: f64,
}

impl<'a> ClusteringConfig<'a> {
    pub fn new(
        band_width: usize,
        copy_num: usize,
        coverage: f64,
        local_coverage: f64,
        gains: &'a Gains,
    ) -> Self {
        Self {
            band_width,
            coverage,
            gains,
            copy_num,
            local_coverage,
        }
    }
}

fn modification_table<T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &[Vec<kiley::Op>],
    strands: &[bool],
    band: usize,
    hmm: &kiley::hmm::PairHiddenMarkovModelOnStrands,
) -> Vec<Vec<f64>> {
    reads
        .iter()
        .zip(ops)
        .zip(strands)
        .map(|((seq, op), strand)| {
            let hmm = match strand {
                true => hmm.forward(),
                false => hmm.reverse(),
            };
            let (mut table, lk) =
                hmm.modification_table_antidiagonal(template, seq.borrow(), op, band);
            table.iter_mut().for_each(|x| *x -= lk);
            table
        })
        .collect()
}

fn filter_by<_T>(profiles: &[Vec<f64>], probes: &[(usize, _T)]) -> Vec<Vec<f64>> {
    profiles
        .iter()
        .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
        .collect()
}

pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &[Vec<kiley::Op>],
    strands: &[bool],
    rng: &mut R,
    hmm: &kiley::hmm::PairHiddenMarkovModelOnStrands,
    config: &ClusteringConfig,
) -> ClusteringDevResult {
    if config.copy_num < 2 {
        return (vec![0; reads.len()], vec![vec![0f64]; reads.len()], 0f64, 1);
    }
    let feature_vectors = search_variants(template, reads, ops, strands, hmm, config);
    let clustering_result = cluster_filtered_variants(&feature_vectors, config, rng);
    let (mut assignments, mut likelihood_gains, max, max_k) = clustering_result;
    if log_enabled!(log::Level::Trace) {
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
    (assignments, likelihood_gains, max, max_k)
}

pub fn search_variants<T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    ops: &[Vec<kiley::Op>],
    strands: &[bool],
    hmm: &kiley::hmm::PairHiddenMarkovModelOnStrands,
    config: &ClusteringConfig,
) -> FeatureVector {
    let profiles = modification_table(template, reads, ops, strands, config.band_width, hmm);
    let profiles = compress_small_gains(profiles, template, config.gains);
    let probes = filter_profiles(template, &profiles, strands, config);
    let op_and_homop = operation_and_homopolymer_length(template, &probes);
    let variants = filter_by(&profiles, &probes);
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
    }
    (variants, op_and_homop)
}

const MIN_REQ_FRACTION: f64 = 0.5;
fn compress_small_gains(
    mut profiles: Vec<Vec<f64>>,
    template: &[u8],
    gains: &Gains,
) -> Vec<Vec<f64>> {
    if profiles.is_empty() {
        return profiles;
    }
    let homopolymer_length = homopolymer_length(template);
    let min_req: Vec<_> = (0..profiles[0].len())
        .map(|pos| {
            let (bp, diff_type) = pos_to_bp_and_difftype(pos);
            let homop_len = *homopolymer_length.get(bp).unwrap_or(&1);
            gains.expected(homop_len, diff_type) * MIN_REQ_FRACTION
        })
        .collect();
    for prof in profiles.iter_mut() {
        for (&min_req, x) in min_req.iter().zip(prof.iter_mut()) {
            if x.abs() < min_req {
                *x = 0f64;
            }
        }
    }
    profiles
}

use crate::likelihood_gains::DiffType;
fn pos_to_bp_and_difftype(pos: usize) -> (usize, DiffType) {
    let (bp, op) = (pos / NUM_ROW, pos % NUM_ROW);
    let diff_type = if op < 4 {
        DiffType::Subst
    } else if op < 8 + kiley::hmm::COPY_SIZE {
        DiffType::Ins
    } else {
        DiffType::Del
    };
    (bp, diff_type)
}

fn operation_and_homopolymer_length<_T>(
    template: &[u8],
    probes: &[(usize, _T)],
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

pub fn cluster_filtered_variants<R: Rng>(
    (variants, variant_type): &FeatureVector,
    config: &ClusteringConfig,
    rng: &mut R,
) -> ClusteringDevResult {
    let copy_num = config.copy_num;
    let coverage = config.coverage;
    let gains = config.gains;
    if copy_num <= 1 || variants.iter().all(|xs| xs.is_empty()) || variants.len() <= copy_num {
        let asn = vec![0; variants.len()];
        let lk_gains = vec![vec![0f64; 1]; variants.len()];
        return (asn, lk_gains, 0f64, 1);
    }
    let datasize = variants.len();
    let per_cluster_cov = config.local_coverage;
    let (mut assignments, mut max, mut max_k, mut read_lk_gains) =
        (vec![0; datasize], 0f64, 1, vec![0f64; datasize]);
    let mut prev_used_columns = vec![false; variants[0].len()];
    let range = {
        let end = copy_num.min(1 + 2 * variant_type.len());
        let start = end.max(5) - 3;
        start..=end
    };
    trace!("RANGE\t{:?}", range);
    for k in range {
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
        let expected_gain = expected_gain_per_read * per_cluster_cov + 0.1;
        trace!("LK\t{k}\t{score:.3}\t{expected_gain:.3}\t{improved_reads}");
        if expected_gain < score - max {
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
    (assignments, likelihood_gains, max, max_k)
}

fn min_gain(gains: &Gains, variant_type: &[(usize, DiffType)], used_columns: &[bool]) -> f64 {
    std::iter::zip(variant_type, used_columns)
        .filter_map(|(&(homop_len, diff_type), is_used)| match is_used {
            true => Some(gains.expected(homop_len, diff_type) / 3f64),
            false => None,
        })
        .min_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(1f64)
}

const EXPT_GAIN_FACTOR: f64 = 0.8;
fn expected_gains(
    gains: &Gains,
    variant_type: &[(usize, DiffType)],
    prev_columns: &[bool],
    used_columns: &[bool],
) -> f64 {
    assert_eq!(variant_type.len(), used_columns.len());
    let no_new_variants = prev_columns == used_columns;
    // Previously not used, currently used.
    let newly_used = std::iter::zip(prev_columns, used_columns).map(|(&p, &c)| !p & c);
    let check_column = newly_used.map(|b| b | no_new_variants);
    let expt_gain = std::iter::zip(variant_type, check_column)
        .map(|(&(homop, diff_type), is_used)| match is_used {
            true => gains.expected(homop, diff_type),
            false => 0.0000001,
        })
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(0f64);
    (EXPT_GAIN_FACTOR * expt_gain).max(0.1)
}

fn count_improved_reads(new_gains: &[f64], old_gains: &[f64], min_gain: f64) -> usize {
    std::iter::zip(new_gains.iter(), old_gains.iter())
        .filter(|(&new, &old)| old + min_gain < new)
        .count()
}

fn is_explainable_by_strandedness<I: std::iter::Iterator<Item = (f64, bool)>>(profiles: I) -> bool {
    let mut strand_count = [0; 2];
    let mut sign_count = [0; 2];
    let mut obs_count = [[0; 2]; 2];
    for (lkdiff, strand) in profiles.filter(|(lk, _)| lk.abs() > 0.0001f64) {
        strand_count[strand as usize] += 1;
        sign_count[lkdiff.is_sign_positive() as usize] += 1;
        obs_count[strand as usize][lkdiff.is_sign_positive() as usize] += 1;
    }
    let sum: usize = strand_count.iter().sum();
    if sum == 0 {
        return false;
    }
    trace!("RAWCOUNT\t{obs_count:?}");
    let chisq: f64 = std::iter::zip(&obs_count, &strand_count)
        .map(|(obs_sign, strand)| {
            std::iter::zip(obs_sign, &sign_count)
                .map(|(&obs, sign)| {
                    let expected = (strand * sign) as f64 / sum as f64;
                    (obs as f64 - expected).powi(2) / expected
                })
                .sum::<f64>()
        })
        .sum();
    chisq < 10f64
}

// LK->LK-logsumexp(LK).
fn to_posterior_probability(lks: &mut [Vec<f64>]) {
    for xs in lks.iter_mut() {
        let total = crate::misc::logsumexp(xs);
        xs.iter_mut().for_each(|x| *x -= total);
    }
}

// i->k->the likelihood gain of the i-th read when clustered in the k-th cluster.
// `copy_num` is the copy number of this chunk, not the *cluster number*.
// Usually, they are the same but sometimes there are exact repeats, and the number of the cluster
// would be smaller than the copy number.
fn get_likelihood_gain(variants: &[Vec<f64>], assignments: &[usize], k: usize) -> Vec<Vec<f64>> {
    let mut lks = vec![vec![LKCount::default(); variants[0].len()]; k];
    let mut clusters = vec![0; k];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        clusters[asn] += 1;
        for (slot, x) in lks[asn].iter_mut().zip(vars.iter()) {
            slot.add(*x);
        }
    }
    let use_columns = get_used_columns(&lks);
    trace!("FILTER\t{use_columns:?}");
    variants
        .iter()
        .map(|vars| {
            lks.iter()
                .map(|slots| -> f64 {
                    vars.iter()
                        .zip(slots.iter())
                        .zip(use_columns.iter())
                        .filter(|&((_, lkc), to_use)| to_use & (POS_THR < lkc.total_gain))
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
    let mut lks = vec![vec![LKCount::default(); variants[0].len()]; k];
    let mut clusters = vec![0; k];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        clusters[asn] += 1;
        for (slot, x) in lks[asn].iter_mut().zip(vars.iter()) {
            slot.add(*x);
        }
    }
    let use_columns = get_used_columns(&lks);
    let gain_on_read: Vec<_> = variants
        .iter()
        .zip(assignments.iter())
        .map(|(vars, &asn)| {
            vars.iter()
                .zip(lks[asn].iter())
                .zip(use_columns.iter())
                .filter(|&((_, lkc), to_use)| to_use & (POS_THR < lkc.total_gain))
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

// ROUND * cluster num variants would be selected.
const ROUND: usize = 3;
const PVALUE: f64 = 0.05;
// False positive rate to determine the strand bias. In other words,
// The probability that the variant is regarded as biased even if it is not
// is 0.05. It essentially sacrifice 5% variants under the name of the strand bias.
fn filter_profiles<T: std::borrow::Borrow<[f64]>>(
    template: &[u8],
    profiles: &[T],
    strands: &[bool],
    config: &ClusteringConfig,
) -> Vec<(usize, f64)> {
    let cluster_num = config.copy_num;
    let coverage = config.coverage;
    let gains = config.gains;
    trace!("\n{gains}");
    let pvalues = gains.pvalues(profiles.len());
    let homopolymer_length = homopolymer_length(template);
    let total_improvement = column_sum(profiles);
    let temp_len = total_improvement.len() / NUM_ROW;
    let probes: Vec<(usize, f64)> = total_improvement
        .iter()
        .enumerate()
        .filter(|&(pos, _)| {
            let (pos, _) = pos_to_bp_and_difftype(pos);
            MASK_LENGTH <= pos && pos <= (temp_len - MASK_LENGTH)
        })
        .filter(|&(pos, _)| pos % NUM_ROW < 8 || pos % NUM_ROW == 8 + kiley::hmm::COPY_SIZE) // 8 -> 1 length copy = same as insertion
        .filter(|&(pos, _)| is_in_short_homopolymer(pos, &homopolymer_length, template))
        .filter(|&(pos, &improve)| {
            has_small_pvalue(pos, improve, &homopolymer_length, &pvalues, gains, temp_len)
        })
        .filter(|&(pos, _)| {
            let lks = profiles.iter().map(|p| p.borrow()[pos]);
            let paired = lks.zip(strands.iter().copied());
            is_explainable_by_strandedness(paired)
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
    let pvalue = template_len as f64 * pvalue;
    let (pos, ed) = pos_to_bp_and_difftype(pos);
    let homop = &homopolymer_length[pos - 1..=pos + 1];
    if pvalue < PVALUE / template_len as f64 {
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
    // 0-> not selected yet. 1-> selected. 2-> removed because of co-occurence, 3-> temporary supressed.
    // Currently, find and sweep. So, to select one variant,
    // It takes O(L) time to find a variant, and O(L) time to sweep linked variant.
    // So, in total, O(ML) time to select M variants. It is OK I think, because
    // usually L is 2K, M is around 3-20.
    let mut is_selected = vec![0; probes.len()];
    for _ in 0..ROUND {
        is_selected
            .iter_mut()
            .filter(|b| **b == 3)
            .for_each(|b| *b = 0);
        for _ in 0..cluster_num.max(2) {
            let next_var_idx = match find_next_variants(probes, &is_selected) {
                Some(idx) => idx,
                None => break,
            };
            let (picked_pos, lk) = probes[next_var_idx];
            let (picked_pos_in_bp, ed) = pos_to_bp_and_difftype(picked_pos);
            trace!("PICK\t{picked_pos_in_bp}\t{ed}\t{lk:.3}");
            is_selected[next_var_idx] = 1;
            for (&(pos, _), selected) in probes
                .iter()
                .zip(is_selected.iter_mut())
                .filter(|&(_, &mut selected)| selected == 0 || selected == 3)
            {
                let (pos_in_bp, _) = pos_to_bp_and_difftype(pos);
                let diff_in_bp = pos_in_bp.max(picked_pos_in_bp) - pos_in_bp.min(picked_pos_in_bp);
                if diff_in_bp < MASK_LENGTH {
                    trace!("REMOVE\tPERM\t{pos_in_bp}\t{picked_pos_in_bp}\t{diff_in_bp}");
                    *selected = 2;
                } else {
                    let sok_sim = sokal_michener(profiles, picked_pos, pos);
                    let cos_sim = cosine_similarity(profiles, picked_pos, pos);
                    if 0.8 < sok_sim || 0.8 < cos_sim.abs() {
                        *selected = 3;
                    }
                }
                // let sok_sim = sokal_michener(profiles, picked_pos, pos);
                // let cos_sim = cosine_similarity(profiles, picked_pos, pos);
                // if 0.99 < sok_sim || 0.99 < cos_sim.abs() || diff_in_bp < MASK_LENGTH {
                //     trace!("REMOVE\tPERM\t{pos_in_bp}\t{picked_pos_in_bp}\t{sok_sim:.3}\t{cos_sim:.3}");
                //     *selected = 2;
                // } else if 0.8 < cos_sim.abs() {
                //     trace!("REMOVE\tTEMP\t{pos_in_bp}\t{picked_pos_in_bp}\t{sok_sim:.3}\t{cos_sim:.3}");
                //     *selected = 3;
                // }
            }
        }
    }
    probes
        .iter()
        .zip(is_selected)
        .filter_map(|(x, y)| (y == 1).then_some(*x))
        .collect()
}

fn column_sum<T: std::borrow::Borrow<[f64]>>(profiles: &[T]) -> Vec<(f64, usize)> {
    let mut total_improvement = vec![(0.0, 0); profiles[0].borrow().len()];
    for prof in profiles.iter().map(|x| x.borrow()) {
        for (slot, &gain) in total_improvement.iter_mut().zip(prof) {
            if POS_THR < gain {
                slot.0 += gain;
                slot.1 += 1;
            }
        }
    }
    total_improvement
}

fn find_next_variants(probes: &[(usize, f64)], is_selected: &[u8]) -> Option<usize> {
    probes
        .iter()
        .map(|x| x.1)
        .zip(is_selected.iter())
        .enumerate()
        .filter(|&(_, (_, &flag))| flag == 0)
        .map(|(idx, (lk, _))| (idx, lk))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .map(|x| x.0)
}

fn cosine_similarity<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> f64 {
    let (inner_prod, isumsq, jsumsq) = profiles
        .iter()
        .map(|profile| profile.borrow())
        .map(|profile| (profile[i], profile[j]))
        .filter(|(x, y)| POS_THR < x.abs() && POS_THR < y.abs())
        .fold((0f64, 0f64, 0f64), |(ip, isumsq, jsumsq), (x, y)| {
            (ip + x * y, isumsq + x * x, jsumsq + y * y)
        });
    if isumsq == 0f64 {
        return 0f64;
    }
    inner_prod / isumsq.sqrt() / jsumsq.sqrt()
}

// Return Max-Sokal-Michener similarity of the i-th column and j-th (flipped) column.
fn sokal_michener<T: std::borrow::Borrow<[f64]>>(profiles: &[T], i: usize, j: usize) -> f64 {
    let (matches, mismatches) = profiles
        .iter()
        .map(|prof| (prof.borrow()[i], prof.borrow()[j]))
        .filter(|&(x, y)| POS_THR < x.abs() && POS_THR < y.abs())
        .fold((0, 0), |(mat, mism), (x, y)| match 0f64 < x * y {
            true => (mat + 1, mism),
            false => (mat, mism + 1),
        });
    let total = matches + mismatches;
    if total == 0 {
        0f64
    } else {
        mismatches.max(matches) as f64 / total as f64
    }
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

use rand::seq::IteratorRandom;

// Return the maximum likelihood.
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).
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
    let mut lks = vec![vec![LKCount::default(); data[0].len()]; k];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            lk.add(*x);
        }
    }
    let mut lk = get_lk(&lks, &clusters, &size_to_lk);
    let (mut max, mut argmax) = (lk, assign.to_vec());
    let total = 2000 * data.len();
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
    let mut lks = vec![vec![LKCount::default(); data[0].len()]; k];
    for (xs, &asn) in data.iter().zip(assign.iter()) {
        clusters[asn] += 1;
        for (lk, x) in lks[asn].iter_mut().zip(xs) {
            lk.add(*x);
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
    lks: &mut [Vec<LKCount>],
    clusters: &mut [usize],
) {
    let from = assign[idx];
    // Remove data[idx] from lks[from] and clusters[from.]
    clusters[from] -= 1;
    for (lks, x) in lks[from].iter_mut().zip(data[idx].iter()) {
        lks.sub(*x);
    }
    assign[idx] = to;
    clusters[to] += 1;
    for (lks, x) in lks[to].iter_mut().zip(data[idx].iter()) {
        lks.add(*x);
    }
}

fn get_lk(lks: &[Vec<LKCount>], clusters: &[usize], _size_to_lk: &[f64]) -> f64 {
    // If true, then that column would not be used.
    let use_columns = get_used_columns(lks);
    let mut lk: f64 = clusters.iter().map(|&size| _size_to_lk[size]).sum();
    for lks in lks.iter() {
        for (lkc, _) in lks.iter().zip(use_columns.iter()).filter(|x| *x.1) {
            lk += lkc.total_gain.max(0f64);
        }
    }
    lk
}

#[derive(Debug, Clone)]
struct LKCount {
    total_gain: f64,
    num_pos: usize,
    num_neg: usize,
    num_zero: usize,
}

impl std::default::Default for LKCount {
    fn default() -> Self {
        Self {
            total_gain: 0f64,
            num_pos: 0,
            num_neg: 0,
            num_zero: 0,
        }
    }
}
// TODO: Can we determine FRAC and IN_POS_RATIO from the data?
// These values are critical.
impl LKCount {
    fn is_informative(&self) -> bool {
        const POS_FRAC: f64 = 0.70;
        let cov = (self.num_pos + self.num_neg) as f64 + 0.0000001;
        0f64 < self.total_gain && POS_FRAC < (self.num_pos as f64) / cov
    }
    fn add(&mut self, x: f64) {
        self.total_gain += x;
        if POS_THR < x {
            self.num_pos += 1;
        } else if x < -POS_THR {
            self.num_neg += 1;
        } else {
            assert!(x.abs() < POS_THR);
            self.num_zero += 1;
        }
    }
    fn sub(&mut self, x: f64) {
        self.total_gain -= x;
        if POS_THR < x {
            self.num_pos -= 1;
        } else if x < -POS_THR {
            self.num_neg -= 1;
        } else {
            assert!(x.abs() < POS_THR);
            self.num_zero -= 1;
        }
    }
}

fn get_used_columns(lks: &[Vec<LKCount>]) -> Vec<bool> {
    let mut to_uses = vec![false; lks[0].len()];
    for lks in lks.iter() {
        assert_eq!(to_uses.len(), lks.len());
        for (to_use, lkc) in to_uses.iter_mut().zip(lks.iter()) {
            *to_use |= lkc.is_informative();
        }
    }
    const IN_POS_RATIO: f64 = 2f64;
    for (d, to_use) in to_uses.iter_mut().enumerate() {
        let column = lks.iter().map(|lks| &lks[d]);
        let pos_in_use: usize = column
            .clone()
            .filter_map(|x| (0f64 < x.total_gain).then_some(x.num_pos))
            .sum();
        let pos_in_neg: usize = column
            .clone()
            .filter_map(|x| (x.total_gain <= 0f64).then_some(x.num_pos))
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
