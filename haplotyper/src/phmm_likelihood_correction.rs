use definitions::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Default)]
pub struct CorrectionConfig {}
pub trait AlignmentCorrection {
    fn correct_clustering(&mut self, config: &CorrectionConfig);
    fn correct_clustering_selected(&mut self, selection: &HashSet<u64>, config: &CorrectionConfig);
}
impl AlignmentCorrection for DataSet {
    fn correct_clustering(&mut self, config: &CorrectionConfig) {
        let present: HashSet<_> = self
            .encoded_reads
            .iter()
            .flat_map(|r| r.nodes.iter())
            .map(|n| n.unit)
            .collect();
        self.selected_chunks.retain(|c| present.contains(&c.id));
        let selections: HashSet<_> = self
            .selected_chunks
            .iter()
            .filter(|c| 1 < c.cluster_num)
            .map(|c| c.id)
            .collect();
        self.correct_clustering_selected(&selections, config);
        self.sanity_check();
    }
    fn correct_clustering_selected(
        &mut self,
        selections: &HashSet<u64>,
        config: &CorrectionConfig,
    ) {
        let copy_numbers = estimate_copy_number_of_cluster(self);
        let corrected_clusterings: Vec<_> = self
            .selected_chunks
            .par_iter()
            .filter(|c| 1 < c.cluster_num && selections.contains(&c.id))
            .map(|c| (c.id, (c.cluster_num, c.copy_num)))
            .map(|(id, cluster_copy)| correct_unit(self, id, cluster_copy, &copy_numbers, config))
            .collect();
        let protected = get_protected_clusterings(self);
        let corrected_clustering_on_read = {
            let mut chunks_mut_ref: HashMap<_, _> =
                self.selected_chunks.iter_mut().map(|c| (c.id, c)).collect();
            let mut corrected_clustering_on_read: HashMap<_, Vec<_>> = HashMap::new();
            for (correcteds, (uid, cluster_num)) in corrected_clusterings {
                let chunk = chunks_mut_ref.get_mut(&uid).unwrap();
                let (score, prev) = (chunk.score, chunk.cluster_num);
                if cluster_num == 1 && protected.contains(&uid) {
                    debug!("SQUISHED\t{uid}\t{cluster_num}\t{prev}\t{score}\tP");
                    continue;
                }
                assert!(cluster_num <= chunk.copy_num);
                chunk.cluster_num = cluster_num;
                for (id, idx, asn) in correcteds {
                    corrected_clustering_on_read
                        .entry(id)
                        .or_default()
                        .push((idx, asn));
                }
            }
            corrected_clustering_on_read
        };
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = corrected_clustering_on_read.get(&read.id) {
                for &(pos, asn) in corrected.iter() {
                    read.nodes[pos].cluster = asn;
                }
            }
        }
    }
}

fn get_protected_clusterings(ds: &mut DataSet) -> HashSet<u64> {
    let mut coverage: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *coverage.entry(node.unit).or_default() += 1;
    }
    if ds.model_param.is_none() {
        crate::model_tune::update_model(ds);
    }
    let hmm = crate::model_tune::get_model(ds).unwrap();
    // ds.error_rate()
    let gain = crate::local_clustering::estimate_minimum_gain(&hmm);
    debug!("POLISHED\tMinGain\t{gain:.3}");
    ds.selected_chunks
        .iter()
        .filter(|c| coverage.contains_key(&c.id))
        .filter_map(|c| (coverage[&c.id] as f64 * gain / 1.5 < c.score).then(|| c.id))
        .collect()
}

fn estimate_copy_number_of_cluster(ds: &DataSet) -> Vec<Vec<f64>> {
    let (copy_number, obs_copy_number) = {
        let max = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
        let mut cps = vec![0; max as usize + 1];
        let mut cls = vec![0; max as usize + 1];
        for u in ds.selected_chunks.iter() {
            cps[u.id as usize] = u.copy_num;
            cls[u.id as usize] = u.cluster_num;
        }
        let mut obs_counts: Vec<Vec<f64>> = cls.into_iter().map(|k| vec![0f64; k]).collect();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            let unit = node.unit as usize;
            // We need to normalize the posterior....
            let total = logsumexp(&node.posterior);
            obs_counts[unit]
                .iter_mut()
                .zip(node.posterior.iter())
                .for_each(|(o, p)| *o += (p - total).exp());
        }
        // We already have coverage, at least I assume so.
        let cov = ds.coverage.unwrap();
        // Normalize the coverage.
        obs_counts.iter_mut().flatten().for_each(|x| *x /= cov);
        (cps, obs_counts)
    };
    obs_copy_number
        .par_iter()
        .zip(copy_number.par_iter())
        .map(|(obs_cp, &total_cp)| {
            let mut est_cp = vec![1f64; obs_cp.len()];
            // Give a copy number to the most needed one...
            for _ in (obs_cp.len()).min(total_cp)..total_cp {
                let (idx, _) = obs_cp
                    .iter()
                    .zip(est_cp.iter())
                    .enumerate()
                    .map(|(i, (&o, &e))| (i, 2f64 * (o - e) - 1f64))
                    .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
                    .unwrap();
                est_cp[idx] += 1f64;
            }
            // let obs: Vec<_> = obs_cp.iter().map(|x| format!("{x:.1}")).collect();
            // let est: Vec<_> = est_cp.iter().map(|x| format!("{x:.0}")).collect();
            // debug!("[{}]\t[{}]\t{total_cp}", obs.join(","), est.join(","));
            est_cp
        })
        .collect()
}

type CorrectionResult = (Vec<(u64, usize, u64)>, (u64, usize));
fn correct_unit(
    ds: &DataSet,
    unit_id: u64,
    cluster_and_copynum: (usize, usize),
    copy_numbers: &[Vec<f64>],
    config: &CorrectionConfig,
) -> CorrectionResult {
    //) -> Vec<(u64, usize, Vec<f64>)> {
    let mut reads = vec![];
    for read in ds.encoded_reads.iter() {
        for (idx, node) in read.nodes.iter().enumerate() {
            if node.unit == unit_id {
                reads.push((idx, read));
            }
        }
    }
    reads.sort_by_cached_key(|&(idx, read)| read.nodes[idx].cluster);
    trace!("CENTER\t{:?}", copy_numbers[unit_id as usize]);
    let len = reads.len();
    trace!("Correction\t{unit_id}\t{cluster_and_copynum:?}\t{len}",);
    let (assignments, k) = clustering(&reads, cluster_and_copynum, copy_numbers, unit_id, config);
    assert_eq!(assignments.len(), reads.len());
    let assignments: Vec<_> = reads
        .into_iter()
        .zip(assignments)
        .map(|((idx, read), asn)| (read.id, idx, asn as u64))
        .collect();
    (assignments, (unit_id, k))
}

fn logsumexp(xs: &[f64]) -> f64 {
    match xs.len() {
        0 => 0f64,
        _ => {
            let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
            let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
            max + sum
        }
    }
}

fn clustering(
    reads: &[(usize, &EncodedRead)],
    (k, _upper_k): (usize, usize),
    copy_numbers: &[Vec<f64>],
    id: u64,
    _: &CorrectionConfig,
) -> (Vec<usize>, usize) {
    let contexts: Vec<_> = reads
        .iter()
        .map(|&(idx, read)| {
            let center = &read.nodes[idx];
            let (up, tail) = match center.is_forward {
                true => {
                    let up: Vec<_> = read.nodes[..idx]
                        .iter()
                        .map(|n| (n.unit, n.posterior.as_slice()))
                        .rev()
                        .collect();
                    let tail: Vec<_> = read.nodes[idx + 1..]
                        .iter()
                        .map(|n| (n.unit, n.posterior.as_slice()))
                        .collect();
                    (up, tail)
                }
                false => {
                    let tail: Vec<_> = read.nodes[..idx]
                        .iter()
                        .map(|n| (n.unit, n.posterior.as_slice()))
                        .rev()
                        .collect();
                    let up: Vec<_> = read.nodes[idx + 1..]
                        .iter()
                        .map(|n| (n.unit, n.posterior.as_slice()))
                        .collect();
                    (up, tail)
                }
            };
            (up, center, tail)
        })
        .collect();
    // for (i, ctx) in contexts.iter().enumerate() {
    //     let cl = argmax(&ctx.1.posterior);
    //     trace!("DUMP\t{}\t{}\t{cl}", i, ctx.1.cluster);
    // }
    let sims: Vec<Vec<_>> = contexts
        .iter()
        .enumerate()
        .map(|(i, ctx)| {
            // let id = ids[i];
            contexts
                .iter()
                .enumerate()
                .map(|(j, dtx)| match i == j {
                    false => alignment(ctx, dtx, copy_numbers),
                    true => 0f64,
                })
                .collect()
        })
        .collect();
    if log_enabled!(log::Level::Trace) {
        let ids: Vec<_> = reads.iter().map(|x| x.1.id).collect();
        for (i, sm) in sims.iter().enumerate() {
            let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
            trace!("SIM\t{i}\t{}\t{}", ids[i], line.join("\t"));
        }
    }
    let laplacian = get_graph_laplacian(&sims);
    let (mut eigens, pick_k) = get_eigenvalues(&laplacian, k, id);
    for (eigen, (idx, read)) in eigens.iter_mut().zip(reads.iter()) {
        assert_eq!(eigen.len(), pick_k);
        let post = &read.nodes[*idx].posterior;
        let total = logsumexp(post);
        eigen.extend(post.iter().map(|x| (x - total).exp()));
    }
    {
        let mut sums: Vec<_> = eigens[0][pick_k..].iter().map(|x| x * x).collect();
        for eig in eigens.iter().skip(1) {
            sums.iter_mut()
                .zip(eig[pick_k..].iter())
                .for_each(|(s, e)| *s += e * e);
        }
        sums.iter_mut().for_each(|x| *x = x.sqrt());
        for eigen in eigens.iter_mut() {
            let posts = eigen.iter_mut().skip(pick_k).zip(sums.iter());
            posts.for_each(|(e, s)| *e /= s);
        }
    }
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    let mut rng = Xoroshiro128PlusPlus::seed_from_u64(id * k as u64);
    // TODO:Maybe we can tune the number of the cluster here, again?
    let cluster_num = k.min(pick_k);
    let asn = (0..10)
        .map(|_| kmeans(&eigens, cluster_num, &mut rng))
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .map(|x| x.0)
        .unwrap();
    // for k in k..=_upper_k {
    //     let (_, dist) = (0..10)
    //         .map(|_| kmeans(&eigens, k, &mut rng))
    //         .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
    //         .unwrap();
    //     debug!("KMEANS\t{id}\t{k}\t{dist}");
    // }
    if log_enabled!(log::Level::Trace) {
        for (i, sm) in eigens.iter().enumerate() {
            let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
            let cls = contexts[i].1.cluster;
            trace!("EIG\t{id}\t{i}\t{}\t{}\t{}", cls, asn[i], line.join("\t"));
        }
    }
    assert!(asn.iter().all(|&x| x <= cluster_num));
    (asn, cluster_num)
}

// Return the *normalized* graph lap.
// I - D^{-1}W, where I is identity matrix, D is the diagonal "summed" matrix and
// W is the similarity matrix itself..
fn get_graph_laplacian(sims: &[Vec<f64>]) -> Vec<Vec<f64>> {
    sims.iter()
        .enumerate()
        .map(|(i, row)| {
            let div: f64 = row.iter().sum::<f64>().recip();
            row.iter()
                .enumerate()
                .map(|(t, w)| if t == i { 1f64 - w * div } else { -w * div })
                .collect()
        })
        .collect()
}

const EIGEN_THR: f64 = 0.25;
fn get_eigenvalues(matrix: &[Vec<f64>], _k: usize, id: u64) -> (Vec<Vec<f64>>, usize) {
    let datalen = matrix.len();
    if datalen == 0 {
        panic!("{}", id)
    }
    if matrix.iter().any(|x| x.is_empty()) {
        panic!("{}", datalen);
    }
    let rows: Vec<_> = matrix
        .iter()
        .map(|row| nalgebra::RowDVector::from(row.to_vec()))
        .collect();
    let matrix = nalgebra::DMatrix::from_rows(&rows);
    let eigens = matrix.clone().symmetric_eigen();
    let mut eigen_and_eigenvec: Vec<_> = eigens
        .eigenvectors
        .column_iter()
        .zip(eigens.eigenvalues.iter())
        .collect();
    eigen_and_eigenvec.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
    // Random thoughts:
    // In this code, I try to determine the *optimal* number of the cluster.
    // One metric is the # of the eigenvalues below threshould. It works fine usually.
    // However, suppose two varaints separated more than 50Kbp. In this case,
    // a few reads span them. Thus, the graph is densely connected.
    // The smallest eigenvalue is zero, but the second smallest eigenvalues
    // would be larger than THR usually.
    let opt_k = eigen_and_eigenvec
        .iter()
        .take_while(|&(_, &lam)| lam < EIGEN_THR)
        .count();
    let pick_k = opt_k;
    if pick_k == 0 {
        for row in matrix.row_iter() {
            eprintln!("{row:?}");
        }
        panic!("{}", opt_k);
    }
    let features: Vec<Vec<_>> = (0..datalen)
        .map(|i| (0..pick_k).map(|j| eigen_and_eigenvec[j].0[i]).collect())
        .collect();
    (features, pick_k)
}

use rand::Rng;

// use crate::stats::Stats;
fn kmeans<R: Rng>(xs: &[Vec<f64>], k: usize, rng: &mut R) -> (Vec<usize>, f64) {
    fn update_assignment(data: &[Vec<f64>], centers: &[Vec<f64>], asn: &mut [usize]) {
        asn.iter_mut().zip(data).for_each(|(asn, datum)| {
            *asn = centers
                .iter()
                .enumerate()
                .map(|(idx, center)| {
                    let dist: f64 = datum
                        .iter()
                        .zip(center)
                        .map(|(d, c)| (d - c) * (d - c))
                        .sum();
                    (idx, dist)
                })
                .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
                .map(|x| x.0)
                .unwrap();
        });
    }
    fn update_centers(data: &[Vec<f64>], centers: &mut [Vec<f64>], asn: &[usize], k: usize) {
        let mut counts = vec![0; k];
        centers.iter_mut().flatten().for_each(|x| *x = 0f64);
        for (datum, &asn) in data.iter().zip(asn.iter()) {
            counts[asn] += 1;
            centers[asn]
                .iter_mut()
                .zip(datum)
                .for_each(|(c, d)| *c += d);
        }
        centers
            .iter_mut()
            .zip(counts.iter())
            .filter(|&(_, &count)| count > 0)
            .for_each(|(center, count)| {
                center.iter_mut().for_each(|c| *c /= *count as f64);
            });
    }
    fn compute_distance(data: &[Vec<f64>], centers: &[Vec<f64>], asn: &[usize]) -> f64 {
        data.iter()
            .zip(asn.iter())
            .map(|(datum, &asn)| -> f64 {
                centers[asn]
                    .iter()
                    .zip(datum.iter())
                    .map(|(c, d)| (c - d) * (c - d))
                    .sum()
            })
            .sum()
    }
    let mut assignments: Vec<_> = xs.iter().map(|_| rng.gen_range(0..k)).collect();
    let dim = xs[0].len();
    let mut centers = vec![vec![0f64; dim]; k];
    let mut dist = compute_distance(xs, &centers, &assignments);
    loop {
        update_centers(xs, &mut centers, &assignments, k);
        update_assignment(xs, &centers, &mut assignments);
        let new_dist = compute_distance(xs, &centers, &assignments);
        if dist - new_dist < 0.0001 {
            break;
        } else {
            dist = new_dist;
        }
    }
    (assignments, dist)
}

type Context<'a> = (Vec<(u64, &'a [f64])>, &'a Node, Vec<(u64, &'a [f64])>);
const LK_CAP: f64 = 15f64;
fn alignment<'a>(
    (up1, center1, down1): &Context<'a>,
    (up2, center2, down2): &Context<'a>,
    copy_numbers: &[Vec<f64>],
) -> f64 {
    assert_eq!(center1.unit, center2.unit);
    let center_copy_num = &copy_numbers[center1.unit as usize];
    let center = sim(&center1.posterior, &center2.posterior, center_copy_num);
    let up_aln = align_swg(up1, up2, copy_numbers);
    let down_aln = align_swg(down1, down2, copy_numbers);
    // This is probably the correct formulation of the similarity score, a.k.a the
    // likelihood ratio of the identity-by-haplotype. Exp(ln(From the same Hap) - ln(1-From the same hap)).
    // However, likelihood ratio test is usually goes to extreme and we want to examine the difference between exp(150) vs exp(151) into
    // similarity matrix.
    // Anyway, use Max(0, log likleihood ratio) instead.
    let likelihood_ratio = up_aln + down_aln + center;
    match likelihood_ratio <= LK_CAP {
        true => likelihood_ratio.exp(),
        false => LK_CAP.exp() + (likelihood_ratio - LK_CAP),
    }
    // let similarity = (up_aln + down_aln + center).min(LK_CAP).exp();
    // assert!(0f64 <= similarity);
    // if log_enabled!(log::Level::Trace) {
    //     fn to_line<'a>(vector: &[(u64, &'a [f64])]) -> Vec<String> {
    //         vector
    //             .iter()
    //             .map(|(u, p)| format!("{u}-{}", argmax(p)))
    //             .collect()
    //     }
    //     trace!("ALN\t\t{}", to_line(up1).join(" "));
    //     trace!("ALN\t\t{}", to_line(up2).join(" "));
    //     trace!("ALN\t\t{}", to_line(down1).join(" "));
    //     trace!("ALN\t\t{}", to_line(down2).join(" "));
    //     let (cl1, cl2) = (argmax(&center1.posterior), argmax(&center2.posterior));
    //     trace!("ALN\t{up_aln:.3}\t{center:.3}\t{down_aln:.3}\t{similarity:.3}\t{cl1}\t{cl2}");
    // }
    // similarity
}

// Align by SWG.
fn align_swg<'a>(
    arm1: &[(u64, &'a [f64])],
    arm2: &[(u64, &'a [f64])],
    copy_numbers: &[Vec<f64>],
) -> f64 {
    const GAP_OPEN: f64 = -0.5f64;
    const GAP_EXTEND: f64 = -100f64;
    const MISM: f64 = -100f64;
    let (len1, len2) = (arm1.len(), arm2.len());
    let lower = (len1 + len2 + 2) as f64 * MISM;
    // Match,Del on arm2, Del on arm1.
    let mut dp = vec![vec![(lower, lower, lower); len2 + 1]; len1 + 1];
    for (i, row) in dp.iter_mut().enumerate().skip(1) {
        row[0].2 = GAP_OPEN + (i - 1) as f64 * GAP_EXTEND;
    }
    for j in 1..len2 + 1 {
        dp[0][j].1 = GAP_OPEN + (j - 1) as f64 * GAP_EXTEND;
    }
    dp[0][0].0 = 0f64;
    fn max((a, b, c): (f64, f64, f64)) -> f64 {
        a.max(b).max(c)
    }
    for (i, (u1, p1)) in arm1.iter().enumerate() {
        let i = i + 1;
        for (j, (u2, p2)) in arm2.iter().enumerate() {
            let j = j + 1;
            let match_score = match u1 == u2 {
                true => sim(p1, p2, &copy_numbers[*u1 as usize]),
                false => MISM,
            };
            // Match
            let mat = max(dp[i - 1][j - 1]) + match_score;
            let del2 = {
                let (mat, d2, d1) = dp[i][j - 1];
                (mat + GAP_OPEN).max(d2 + GAP_EXTEND).max(d1 + GAP_OPEN)
            };
            let del1 = {
                let (mat, d2, d1) = dp[i - 1][j];
                (mat + GAP_OPEN).max(d2 + GAP_OPEN).max(d1 + GAP_EXTEND)
            };
            dp[i][j] = (mat, del2, del1);
        }
    }
    let row_last = dp[len1].iter().map(|&x| max(x));
    let column_last = dp.iter().filter_map(|l| l.last()).map(|&x| max(x));
    row_last
        .chain(column_last)
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap()
}

#[allow(dead_code)]
fn align<'a>(
    arm1: &[(u64, &'a [f64])],
    arm2: &[(u64, &'a [f64])],
    copy_numbers: &[Vec<f64>],
) -> f64 {
    // Allow gap with small penalty.
    // To treat the deletion error, or "encoding error" due to the
    // errors in the sequnece,
    // we treat one-length deletion with small penalty,
    // and deletion longer than one with high penalty.
    // To do this, we, currently, evaluate the alignemnt in post-hoc.
    // I know it is not optimal, nor correct algorithm, but it works (should we hoep more?)
    // TODO: implement NWG algorithm
    const MISM: f64 = -10000f64;
    const GAP: f64 = -0.5f64;
    let mut dp = vec![vec![0f64; arm2.len() + 1]; arm1.len() + 1];
    for (i, _) in arm1.iter().enumerate() {
        dp[i + 1][0] = dp[i][0] + GAP;
    }
    for (j, _) in arm2.iter().enumerate() {
        dp[0][j + 1] = dp[0][j] + GAP;
    }
    for (i, (u1, p1)) in arm1.iter().enumerate() {
        let i = i + 1;
        for (j, (u2, p2)) in arm2.iter().enumerate() {
            let j = j + 1;
            let match_score = match u1 == u2 {
                true => sim(p1, p2, &copy_numbers[*u1 as usize]),
                false => MISM,
            };
            dp[i][j] = (dp[i - 1][j - 1] + match_score)
                .max(dp[i - 1][j] + GAP)
                .max(dp[i][j - 1] + GAP);
        }
    }
    let last_row = dp
        .last()
        .unwrap()
        .iter()
        .enumerate()
        .map(|(j, s)| (arm1.len(), j, s));
    let last_column = dp
        .iter()
        .filter_map(|x| x.last())
        .enumerate()
        .map(|(i, s)| (i, arm2.len(), s));
    let (mut i, mut j, score) = last_row
        .chain(last_column)
        .max_by(|x, y| x.2.partial_cmp(y.2).unwrap())
        .unwrap();
    // 0->Deletion or insertion,1->Match
    let mut ops = vec![];
    while 0 < i && 0 < j {
        let current = dp[i][j];
        let (u1, p1) = arm1[i - 1];
        let (u2, p2) = arm2[j - 1];
        let match_score = match u1 == u2 {
            true => sim(p1, p2, &copy_numbers[u1 as usize]),
            false => MISM,
        };
        assert!(!match_score.is_nan());
        if (dp[i - 1][j - 1] + match_score - current).abs() < 0.00001 {
            ops.push(1);
            i -= 1;
            j -= 1;
        } else if (dp[i - 1][j] + GAP - current).abs() < 0.000001 {
            ops.push(0);
            i -= 1;
        } else {
            assert!((dp[i][j - 1] + GAP - current).abs() < 0.000001);
            ops.push(0);
            j -= 1;
        }
    }
    const LONG_GAP: f64 = -80f64;
    let gap_pens: f64 = ops
        .split(|&x| x == 1)
        .map(|gaps| (gaps.len().max(1) - 1) as f64 * LONG_GAP)
        .sum();
    *score + gap_pens
}

const MOCK_CP: f64 = 1.5;
pub fn sim(xs: &[f64], ys: &[f64], cps: &[f64]) -> f64 {
    assert_eq!(xs.len(), cps.len());
    assert_eq!(xs.len(), ys.len());
    if cps.len() == 1 {
        let total: f64 = cps.iter().sum();
        return -(total.max(MOCK_CP) - 1f64).ln();
    }
    let iter = xs
        .iter()
        .zip(ys.iter())
        .zip(cps.iter())
        .map(|((x, y), z)| x + y - z.ln());
    let logp = logsumexp_str(iter);
    let logit = logit_from_lnp(logp);
    assert!(!logit.is_infinite(), "{},{}", logp, logit);
    logit
}

// log(p) -> log(p) - log(1-p)
// TODO:Is this OK?
fn logit_from_lnp(lnp: f64) -> f64 {
    assert!(lnp <= 0f64);
    const LOWER_CUT: f64 = -80f64;
    const UPPER_CUT: f64 = 80f64;
    // The same as -exp(-UPPER_CUT)
    const UPPER_THR: f64 = -1.8e-35f64;
    if lnp < LOWER_CUT {
        LOWER_CUT
    } else if UPPER_THR < lnp {
        UPPER_CUT
    } else {
        lnp - f64::ln_1p(-lnp.exp())
    }
}

fn logsumexp_str<I: Iterator<Item = f64>>(xs: I) -> f64 {
    let (mut max, mut accum, mut count) = (std::f64::NEG_INFINITY, 0f64, 0);
    for x in xs {
        count += 1;
        if x <= max {
            accum += (x - max).exp();
        } else {
            accum = (max - x).exp() * accum + 1f64;
            max = x;
        }
    }
    match count {
        1 => max,
        _ => accum.ln() + max,
    }
}

#[allow(dead_code)]
fn argmax(xs: &[f64]) -> usize {
    xs.iter()
        .enumerate()
        .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
        .unwrap()
        .0
}
