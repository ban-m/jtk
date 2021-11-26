//! A small K-means clustering algorithm.
#[allow(dead_code)]
const DIFF_SIZE: usize = 9;
const DEL_SIZE: usize = 4;
const REP_SIZE: usize = 4;
// Average LK gain for one read. If you increment the number of the cluster,
// you should gain AVERAGE_LK * coverage log-likelihood.
const AVERAGE_LK: f64 = 1.1;
// const DEL_LK: f64 = 3f64;
// return expected number of variants under the null hypothesis.
// fn expected_mis_num(cov: usize) -> f64 {
//     cov as f64 * 0.1f64 + 0.35
// }
// First and last `MASK_LENGTH` bases would not be considered in variant calling.
const MASK_LENGTH: usize = 5;
use rand::Rng;

// use crate::assemble::string_graph::consensus;
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

/// Usual "flat(non-recursive)" clustering. Return assignments, template, and LK.
// TODO:Maybe we need hmm in the configuration.
pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<(Vec<u8>, Vec<u8>, f64)> {
    let band = config.band_width;
    let cons_template = kiley::consensus(reads, rng.gen(), 10, band);
    let band_width = 100;
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let hmm = hmm.fit_banded(&cons_template, &reads, band_width);
    clustering_with_template(&cons_template, reads, rng, &hmm, config)
        .map(|(x, y)| (x, cons_template, y))
}

pub fn clustering_with_template<R: Rng, T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    rng: &mut R,
    hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
    config: &ClusteringConfig,
) -> Option<(Vec<u8>, f64)> {
    let ClusteringConfig {
        band_width,
        cluster_num,
        coverage,
    } = *config;
    let profiles = get_profiles(template, hmm, reads, band_width as isize)?;
    let cluster_num = cluster_num as usize;
    let selected_variants: Vec<Vec<_>> = {
        let probes = filter_profiles(&profiles, cluster_num, 3, coverage, template.len());
        if log_enabled!(log::Level::Trace) {
            // DUMP Hot columns.
            for (pos, lk) in probes.iter() {
                let (idx, t) = (pos / 9, pos % 9);
                if idx < template.len() {
                    trace!("POS\t{}\t{}\t{}\tED\t{:.3}", pos, idx, t, lk);
                } else {
                    let idx = pos - 9 * template.len();
                    if idx < (DEL_SIZE - 1) * (template.len() - DEL_SIZE) {
                        let (idx, len) = (idx / (DEL_SIZE - 1), idx % (DEL_SIZE - 1));
                        trace!("POS\t{}\t{}\t{}\tDEL\t{:.3}", pos, idx, len, lk);
                    } else {
                        let idx = idx - (DEL_SIZE - 1) * (template.len() - DEL_SIZE);
                        let (idx, len) = (idx / REP_SIZE, idx % REP_SIZE + 1);
                        trace!("POS\t{}\t{}\t{}\tCP\t{:.3}", pos, idx, len, lk);
                    };
                }
            }
        }
        profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect()
    };
    let num = 10;
    let init_cluster_num = cluster_num.max(4) - 3;
    let (assignments, score) = (init_cluster_num..=cluster_num)
        .filter_map(|k| {
            let (asn, score) = (0..num)
                .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
                .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
            trace!("LK\t{}\t{:.3}", k, score);
            let expected_gain = (k - 1) as f64 * AVERAGE_LK * coverage;
            Some((asn, score - expected_gain))
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    // let selected_variants = filter_suspicious_variants(&selected_variants, &assignments);
    // let (assignments, score) = (init_cluster_num..=cluster_num)
    //     .filter_map(|k| {
    //         let (asn, score) = (0..num)
    //             .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
    //             .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    //         trace!("LK\t{}\t{:.3}", k, score);
    //         let expected_gain = (k - 1) as f64 * AVERAGE_LK * coverage;
    //         Some((asn, score - expected_gain))
    //     })
    //     .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    trace!("MAXLK\t{:.3}", score);
    if log_enabled!(log::Level::Trace) {
        for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
            let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
            trace!("ASN\t{}\t{}\t{}\t{}", cluster_num, id, i, prf.join("\t"));
        }
    }
    Some((assignments, score))
}

pub fn clustering_dev<R: Rng, T: std::borrow::Borrow<[u8]>>(
    template: &[u8],
    reads: &[T],
    rng: &mut R,
    hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
    config: &ClusteringConfig,
) -> Option<(Vec<u8>, Vec<Vec<f64>>, f64)> {
    let ClusteringConfig {
        band_width,
        cluster_num,
        coverage,
    } = *config;
    let profiles = get_profiles(template, hmm, reads, band_width as isize)?;
    let cluster_num = cluster_num as usize;
    let selected_variants: Vec<Vec<_>> = {
        let probes = filter_profiles(&profiles, cluster_num, 3, coverage, template.len());
        if log_enabled!(log::Level::Trace) {
            // DUMP Hot columns.
            for (pos, lk) in probes.iter() {
                let (idx, t) = (pos / 9, pos % 9);
                if idx < template.len() {
                    trace!("POS\t{}\t{}\t{}\tED\t{:.3}", pos, idx, t, lk);
                } else {
                    let idx = pos - 9 * template.len();
                    if idx < (DEL_SIZE - 1) * (template.len() - DEL_SIZE) {
                        let (idx, len) = (idx / (DEL_SIZE - 1), idx % (DEL_SIZE - 1));
                        trace!("POS\t{}\t{}\t{}\tDEL\t{:.3}", pos, idx, len, lk);
                    } else {
                        let idx = idx - (DEL_SIZE - 1) * (template.len() - DEL_SIZE);
                        let (idx, len) = (idx / REP_SIZE, idx % REP_SIZE + 1);
                        trace!("POS\t{}\t{}\t{}\tCP\t{:.3}", pos, idx, len, lk);
                    };
                }
            }
        }
        profiles
            .iter()
            .map(|xs| probes.iter().map(|&(pos, _)| xs[pos]).collect())
            .collect()
    };
    let num = 10;
    let init_cluster_num = cluster_num.max(4) - 3;
    let (assignments, _score) = (init_cluster_num..=cluster_num)
        .filter_map(|k| {
            let (asn, score) = (0..num)
                .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
                .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
            trace!("LK\t{}\t{:.3}", k, score);
            let expected_gain = (k - 1) as f64 * AVERAGE_LK * coverage;
            Some((asn, score - expected_gain))
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    let selected_variants = filter_suspicious_variants(&selected_variants, &assignments);
    let (assignments, score) = (init_cluster_num..=cluster_num)
        .filter_map(|k| {
            let (asn, score) = (0..num)
                .map(|_| mcmc_clustering(&selected_variants, k, coverage, rng))
                .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
            trace!("LK\t{}\t{:.3}", k, score);
            let expected_gain = (k - 1) as f64 * AVERAGE_LK * coverage;
            Some((asn, score - expected_gain))
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())?;
    let posterior_probs = get_posterior_probability(&selected_variants, &assignments);
    trace!("MAXLK\t{:.3}", score);
    if log_enabled!(log::Level::Trace) {
        for (id, (i, prf)) in assignments.iter().zip(selected_variants.iter()).enumerate() {
            let prf: Vec<_> = prf.iter().map(|x| format!("{:.2}", x)).collect();
            trace!("ASN\t{}\t{}\t{}\t{}", cluster_num, id, i, prf.join("\t"));
        }
    }
    Some((assignments, posterior_probs, score))
}

fn get_posterior_probability(variants: &[Vec<f64>], assignments: &[u8]) -> Vec<Vec<f64>> {
    let max_asn = *assignments.iter().max().unwrap() as usize;
    let mut count = vec![0; max_asn + 1];
    let mut total_lk_gain = vec![vec![0f64; variants[0].len()]; max_asn + 1];
    for (vars, &asn) in variants.iter().zip(assignments.iter()) {
        count[asn as usize] += 1;
        for (total, var) in total_lk_gain[asn as usize].iter_mut().zip(vars.iter()) {
            *total += var;
        }
    }
    let is_used_position: Vec<Vec<_>> = total_lk_gain
        .iter()
        .map(|totals| totals.iter().map(|x| x.is_sign_positive()).collect())
        .collect();
    let len = variants.len() as f64;
    let fractions: Vec<_> = count.iter().map(|&x| (x as f64 / len).ln()).collect();
    variants
        .iter()
        .map(|vars| {
            let mut lks: Vec<_> = is_used_position
                .iter()
                .zip(fractions.iter())
                .map(|(ps, f)| {
                    let lk = vars
                        .iter()
                        .zip(ps)
                        .fold(0f64, |acc, (v, &p)| if p { acc + v } else { acc });
                    lk + f
                })
                .collect();
            let total = logsumexp(&lks);
            lks.iter_mut().for_each(|x| *x = (*x - total).exp());
            lks
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
#[allow(dead_code)]
fn filter_suspicious_variants(variants: &[Vec<f64>], assignments: &[u8]) -> Vec<Vec<f64>> {
    let max_asn = *assignments.iter().max().unwrap() as usize;
    let dim = variants[0].len();
    let mut count = vec![0; max_asn + 1];
    let mut neg_count = vec![vec![0; max_asn + 1]; dim];
    for (&asn, vars) in assignments.iter().zip(variants.iter()) {
        count[asn as usize] += 1;
        for (pos, var) in vars.iter().enumerate() {
            neg_count[pos][asn as usize] += var.is_sign_negative() as usize;
        }
    }
    // Return 1 - cumulative distribution of the binomial distribution Binom(x|N,P).
    // In other words, i-> sum_{k=i}^{N} Binom(k|N,p).
    fn cum_dist(n: usize, p: f64) -> Vec<f64> {
        // Return log(n!), log(0!) = log(1) = 0.
        let sumlog = |n| -> f64 { (1..n + 1).map(|x| (x as f64).ln()).sum() };
        let (_, mut cdf) = (0..=n)
            .map(|i| {
                sumlog(n) - sumlog(i) - sumlog(n - i)
                    + i as f64 * p.ln()
                    + (n - i) as f64 * (1f64 - p).ln()
            })
            .fold((0f64, vec![]), |(acc, mut xs), x| {
                xs.push(acc + x.exp());
                (acc + x.exp(), xs)
            });
        cdf.iter_mut().for_each(|x| *x = 1f64 - *x);
        cdf
    }
    // Error rate is 15%. (Every error operation would erase the variant information, in the worst case.)
    let err = 0.15;
    let cum_distibution: Vec<_> = count.iter().map(|&n| cum_dist(n, err)).collect();
    let is_informative: Vec<bool> = neg_count
        .iter()
        .map(|negs| {
            let max_p_value = negs
                .iter()
                .zip(cum_distibution.iter())
                .map(|(&neg, pvalues)| pvalues[neg])
                .fold(std::f64::NEG_INFINITY, |x, y| x.max(y));
            0.01 < max_p_value
        })
        .collect();
    variants
        .iter()
        .map(|vars| {
            vars.iter()
                .zip(is_informative.iter())
                .filter_map(|(&v, &b)| b.then(|| v))
                .collect()
        })
        .collect()
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
    template_len: usize,
    profiles: &[T],
    _cluster_num: usize,
    coverage: f64,
    rng: &mut R,
) -> (usize, Vec<u8>) {
    // Filter profiles by doubling.
    // let var_num = 2 * (cluster_num - 1).max(1) as usize;
    let cluster_num = 2;
    let selected_variants: Vec<Vec<_>> = {
        let probes = filter_profiles(profiles, cluster_num, 3, coverage, template_len);
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
                let (sub_cluster_size, asns) =
                    recursive_clustering(template_len, &profiles, 2, coverage, rng);
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
    hmm: &kiley::gphmm::GPHMM<kiley::gphmm::Cond>,
    reads: &[T],
    band_width: isize,
) -> Option<Vec<Vec<f64>>> {
    use kiley::gphmm::*;
    let template = kiley::padseq::PadSeq::new(template);
    let mut profiles = vec![];
    for read in reads.iter() {
        let read = kiley::padseq::PadSeq::new(read.borrow());
        let prof = match banded::ProfileBanded::new(&hmm, &template, &read, band_width) {
            Some(res) => res,
            None => {
                for read in reads.iter() {
                    let read = String::from_utf8(read.borrow().to_vec()).unwrap();
                    error!("READ\tKMEANS\t{}", read);
                }
                let template = String::from_utf8(template.clone().into()).unwrap();
                error!("TEMPLATE\tKMEANS\t{}", template);
                error!("{}", hmm);
                return None;
            }
        };
        let lk = prof.lk();
        // Modif.
        let mut modif_table = prof.to_modification_table();
        modif_table.truncate(9 * template.len());
        // Del table
        let del_table = prof.to_deletion_table(DEL_SIZE);
        assert_eq!(
            del_table.len(),
            (DEL_SIZE - 1) * (template.len() - DEL_SIZE)
        );
        modif_table.extend(del_table);
        // Copy Table.
        let copy_table = prof.to_copy_table(REP_SIZE);
        assert_eq!(copy_table.len(), REP_SIZE * (template.len() - REP_SIZE));
        modif_table.extend(copy_table);
        // Normalize.
        modif_table.iter_mut().for_each(|x| *x -= lk);
        profiles.push(modif_table);
    }
    Some(profiles)
}

// Select round * (cluster_num-1) variants
// TODO: FIXME
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
    let base_position = |pos: usize| match pos / (9 * template_len) {
        0 => pos / 9,
        _ => {
            let pos = pos - 9 * template_len;
            match pos / ((DEL_SIZE - 1) * (template_len - DEL_SIZE)) {
                0 => pos / (DEL_SIZE - 1),
                _ => (pos - (DEL_SIZE - 1) * (template_len - DEL_SIZE)) / REP_SIZE,
            }
        }
    };
    let in_mask = |pos: usize| {
        let pos = base_position(pos);
        pos < MASK_LENGTH || (template_len - MASK_LENGTH) < pos
    };
    trace!("LKSUM\tpos\tbp\tlen\ttype\tlkdelta");
    let probes: Vec<(usize, f64)> = total_improvement
        .into_iter()
        .enumerate()
        .map(|(pos, (maxgain, count))| {
            let max_lk = (1..cluster_num + 1)
                .map(|k| poisson_lk(count, coverage * k as f64))
                .max_by(|x, y| x.partial_cmp(y).unwrap())
                .unwrap_or_else(|| panic!("{}", cluster_num));
            let (idx, len, ed_type) = if pos / 9 < template_len {
                (pos / 9, pos % 9, "Sn")
            } else {
                let ofs_pos = pos - 9 * template_len;
                if ofs_pos < (DEL_SIZE - 1) * (template_len - DEL_SIZE) {
                    (ofs_pos / (DEL_SIZE - 1), ofs_pos % (DEL_SIZE - 1) + 2, "Dl")
                } else {
                    let ofs_pos = ofs_pos - (DEL_SIZE - 1) * (template_len - DEL_SIZE);
                    (ofs_pos / REP_SIZE, ofs_pos % REP_SIZE + 1, "Cp")
                }
            };
            let total_lk = max_lk + maxgain;
            trace!(
                "LKSUM\t{}\t{}\t{}\t{}\t{}",
                pos,
                idx,
                len,
                ed_type,
                total_lk / profiles.len() as f64,
            );
            (pos, maxgain + max_lk)
        })
        .filter(|&(pos, _)| !in_mask(pos))
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
#[allow(dead_code)]
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
    poisson_prob < binom_prob
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
    // let sokal_norml = mismatches as f64 / (matches + mismatches) as f64;
    // let sokal_flip = matches as f64 / (mismatches + matches) as f64;
    // sokal_norml.max(sokal_flip)
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
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
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
    if k == 1 || data.iter().all(|xs| xs.is_empty()) || data.len() <= k {
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
        let idx = match indices.choose_weighted(rng, |&idx| dists[idx]) {
            Ok(res) => *res,
            Err(why) => {
                for d in data.iter() {
                    error!("{:?}", d);
                }
                panic!("{:?}", why);
            }
        };
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

#[allow(dead_code)]
fn em_clustering<R: Rng>(data: &[Vec<f64>], k: usize, _cov: f64, rng: &mut R) -> (Vec<u8>, f64) {
    if k == 1 || data.iter().all(|xs| xs.is_empty()) || data.len() <= k {
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
        let idx = match indices.choose_weighted(rng, |&idx| dists[idx]) {
            Ok(res) => *res,
            Err(why) => {
                for d in data.iter() {
                    error!("{:?}", d);
                }
                panic!("{:?}", why);
            }
        };
        centers.push(&data[idx]);
    }
    let mut weights: Vec<_> = data
        .iter()
        .map(|xs| {
            let mut ws: Vec<_> = centers
                .iter()
                .map(|center| (1f64 + euclid_norm_f64(center, xs)).recip())
                .collect();
            let sum: f64 = ws.iter().copied().sum();
            for w in ws.iter_mut() {
                *w /= sum;
            }
            ws
        })
        .collect();
    let lk = em_cl(data, &mut weights, k);
    let choises: Vec<_> = (0..k).collect();
    let assignments: Vec<_> = weights
        .iter()
        .map(|ws| *choises.choose_weighted(rng, |&k| ws[k]).unwrap() as u8)
        .collect();
    //for (ws, asn) in weights.iter().zip(assignments.iter()) {
    //     let ws: Vec<_> = ws.iter().map(|x| format!("{:.3}", x)).collect();
    // }
    (assignments, lk)
}

// Return the P(n_1,...,n_k|theta), where the model behide is Chinese Restaurant Process.
// Note that the zero-sized cluster would be supressed.
#[allow(dead_code)]
fn log_partition(clusters: &[usize], theta: f64) -> f64 {
    // Return log theta^K * Prod((n_k-1)!) / theta / (theta+1) / ... / (theta+n-1)
    // where K= # non-zero elements
    let log_theta: f64 = theta.ln();
    fn log_fact(n: usize) -> f64 {
        match n {
            0 => 0f64,
            _ => (1..n + 1).map(|x| (x as f64).ln()).sum(),
        }
    }
    let numerator: f64 = clusters
        .iter()
        .filter(|&&x| x != 0)
        .map(|&x| log_theta + log_fact(x - 1))
        .sum();
    let total: usize = clusters.iter().copied().sum();
    let denominator: f64 = (0..total).map(|i| (theta + i as f64).ln()).sum();
    numerator - denominator
}

// Return the maximum likelihood.
// In this function, we sample the assignemnts of the data, then evaluate the
// likelihood of the probability, P(X|Z) = max_O P(X|Z,O). Here, the parameters O can be analytically resolvable,
// and we can get the maximum.
// Since we have P(Z|X,O) = P(X|Z,O) * P(Z|O) / P(X|O) ~ P(X|Z,O) * P(Z|O), (not precise?)
// We can flip Z->Z' by comparing P(Z'|X)/P(Z|X).
// Finally, the retuned value would be P(Z,X) = P(Z) P(X|Z).
fn mcmc<R: Rng>(data: &[Vec<f64>], assign: &mut [u8], k: usize, cov: f64, rng: &mut R) -> f64 {
    let binom_lks: Vec<_> = (0..=data.len())
        .map(|x| max_binom_lk(x, data.len(), cov / data.len() as f64))
        .collect();
    let get_partition_lk =
        |clusters: &[usize]| -> f64 { clusters.iter().map(|&x| binom_lks[x]).sum::<f64>() };
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
    // Current Likelihood of the data.
    let data_lk: f64 = lks
        .iter()
        .map(|xs| xs.iter().map(|x| x.max(0f64)).sum::<f64>())
        .sum();
    let partition_lk = get_partition_lk(&clusters);
    let mut lk = data_lk + partition_lk;
    let total = 100_000 * k as usize;
    // MAP estimation.
    let (mut max, mut argmax) = (f64::NEG_INFINITY, vec![]);
    for _ in 0..total {
        let (idx, new) = (rng.gen_range(0..data.len()), rng.gen_range(0..k));
        let (old, xs) = (assign[idx], &data[idx]);
        if old == new as u8 {
            continue;
        }
        let old = old as usize;
        // Change the assignment.
        assign[idx] = new as u8;
        clusters[old] -= 1;
        clusters[new] += 1;
        for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
            *lk -= x;
        }
        for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
            *lk += x;
        }
        // Calc likelihoods.
        let new_data_lk: f64 = lks
            .iter()
            .map(|xs| xs.iter().map(|x| x.max(0f64)).sum::<f64>())
            .sum();
        let new_part_lk: f64 = get_partition_lk(&clusters);
        let new_lk = new_data_lk + new_part_lk;
        let prob = (new_lk - lk).exp().min(1f64);
        if rng.gen_bool(prob) {
            lk = new_lk;
            if max < lk {
                max = lk;
                argmax = assign.to_vec();
            }
        } else {
            // Flip fail. Roll back the the previous state.
            assign[idx] = old as u8;
            clusters[old] += 1;
            clusters[new] -= 1;
            for (lk, x) in lks[old].iter_mut().zip(xs.iter()) {
                *lk += x;
            }
            for (lk, x) in lks[new].iter_mut().zip(xs.iter()) {
                *lk -= x;
            }
        }
    }
    assign.iter_mut().zip(argmax).for_each(|(x, y)| *x = y);
    max
}

// Return the maximum likelihood.
// Traditional EM clustering.
fn em_cl(data: &[Vec<f64>], weights: &mut [Vec<f64>], k: usize) -> f64 {
    fn model_fraction(
        data: &[Vec<f64>],
        weights: &[Vec<f64>],
        k: usize,
    ) -> (Vec<Vec<bool>>, Vec<f64>) {
        let mut fractions = vec![0f64; k];
        for ws in weights.iter() {
            for (c, w) in fractions.iter_mut().zip(ws.iter()) {
                *c += w;
            }
        }
        for w in fractions.iter_mut() {
            *w /= data.len() as f64;
        }
        // Current (un-modified) likelihoods.
        let mut lks = vec![vec![0f64; data[0].len()]; k];
        for (xs, ws) in data.iter().zip(weights.iter()) {
            for (w, lks) in ws.iter().zip(lks.iter_mut()) {
                for (lk, x) in lks.iter_mut().zip(xs.iter()) {
                    *lk += x * w;
                }
            }
        }
        let models: Vec<Vec<bool>> = lks
            .iter()
            .map(|lks| lks.iter().map(|lk| lk.is_sign_positive()).collect())
            .collect();
        (models, fractions)
    }
    fn get_lk(data: &[Vec<f64>], models: &[Vec<bool>], fractions: &[f64]) -> f64 {
        data.iter()
            .map(|xs| -> f64 {
                models
                    .iter()
                    .map(|ms| -> f64 {
                        xs.iter()
                            .zip(ms.iter())
                            .filter_map(|(x, &b)| b.then(|| x))
                            .sum()
                    })
                    .zip(fractions.iter())
                    .map(|(lk, w)| lk + w.ln())
                    .sum()
            })
            .sum()
    }
    fn update_weight(weights: &mut [f64], xs: &[f64], models: &[Vec<bool>], fractions: &[f64]) {
        for ((m, f), w) in models.iter().zip(fractions.iter()).zip(weights.iter_mut()) {
            *w = xs.iter().zip(m).filter_map(|(&x, &b)| b.then(|| x)).sum();
            *w += f.ln();
        }
        let total = logsumexp(weights);
        for w in weights.iter_mut() {
            *w = (*w - total).exp();
        }
    }
    // how many instance in a cluster.
    let (mut models, mut fractions) = model_fraction(data, &weights, k);
    let mut lk = get_lk(data, &models, &fractions);
    loop {
        for (ws, xs) in weights.iter_mut().zip(data.iter()) {
            update_weight(ws, xs, &models, &fractions);
        }
        let (new_models, new_fractions) = model_fraction(data, &weights, k);
        let new_lk = get_lk(data, &new_models, &new_fractions);
        // debug!("LK\t{:.3}\t{:.3}", lk, new_lk);
        if lk + 0.00001 < new_lk {
            lk = new_lk;
            models = new_models;
            fractions = new_fractions;
        } else {
            break;
        }
    }
    lk
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
