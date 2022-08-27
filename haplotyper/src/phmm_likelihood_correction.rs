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
        let supress_cluster = supress_threshold(&corrected_clusterings);
        debug!("SUPRESS\t{supress_cluster:.3}");
        let corrected_clustering_on_read = {
            let mut chunks_mut_ref: HashMap<_, _> =
                self.selected_chunks.iter_mut().map(|c| (c.id, c)).collect();
            let mut corrected_clustering_on_read: HashMap<_, Vec<_>> = HashMap::new();
            for (correcteds, ari, (uid, cluster_num)) in corrected_clusterings {
                let chunk = chunks_mut_ref.get_mut(&uid).unwrap();
                let (score, prev) = (chunk.score, chunk.cluster_num);
                assert!(cluster_num <= chunk.copy_num);
                let supress = cluster_num == 1 || (ari < supress_cluster);
                if supress && protected.contains(&uid) {
                    debug!("PROTECT\t{uid}\t{cluster_num}\t{prev}\t{score}");
                    continue;
                } else if supress {
                    chunk.cluster_num = 1;
                    for (id, idx, _) in correcteds {
                        corrected_clustering_on_read
                            .entry(id)
                            .or_default()
                            .push((idx, 0));
                    }
                } else {
                    chunk.cluster_num = cluster_num;
                    for (id, idx, asn) in correcteds {
                        corrected_clustering_on_read
                            .entry(id)
                            .or_default()
                            .push((idx, asn));
                    }
                }
            }
            corrected_clustering_on_read
        };
        let cluster_num: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = corrected_clustering_on_read.get(&read.id) {
                for &(pos, asn) in corrected.iter() {
                    let node = read.nodes.get_mut(pos).unwrap();
                    node.cluster = asn;
                    let cl = cluster_num[&node.unit];
                    node.posterior.clear();
                    node.posterior.extend(std::iter::repeat(-10000f64).take(cl));
                    node.posterior[asn as usize] = 0f64;
                }
            }
        }
    }
}

const ADJ_RAND_QUANTILE: f64 = 0.05;
fn supress_threshold(clusterings: &[CorrectionResult]) -> f64 {
    let mut adj_rand_indicies: Vec<_> = clusterings.iter().map(|x| x.1).collect();
    adj_rand_indicies.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let pick = (adj_rand_indicies.len() as f64 * ADJ_RAND_QUANTILE).ceil() as usize;
    adj_rand_indicies[pick]
}

const PROTECT_FACTOR: f64 = 1f64;
fn get_protected_clusterings(ds: &mut DataSet) -> HashSet<u64> {
    let mut coverage: HashMap<_, u32> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *coverage.entry(node.unit).or_default() += 1;
    }
    if ds.model_param.is_none() {
        crate::model_tune::update_model(ds);
    }
    let hmm = crate::model_tune::get_model(ds).unwrap();
    let gain = crate::likelihood_gains::estimate_minimum_gain(&hmm) * PROTECT_FACTOR;
    debug!("POLISHED\tMinGain\t{gain:.3}");
    ds.selected_chunks
        .iter()
        .filter_map(|c| {
            let cov = *coverage.get(&c.id)? as f64;
            let cluster_num = c.cluster_num as f64;
            let improve_frac = (cluster_num - 1f64) / cluster_num;
            (cov * improve_frac * gain < c.score).then(|| c.id)
        })
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
            let total = crate::misc::logsumexp(&node.posterior);
            obs_counts[unit]
                .iter_mut()
                .zip(node.posterior.iter())
                .for_each(|(o, p)| *o += (p - total).exp());
        }
        // Normalize the coverage.
        // obs_counts.iter_mut().flatten().for_each(|x| *x /= cov);
        (cps, obs_counts)
    };
    let cov = ds.coverage.unwrap();
    obs_copy_number
        .par_iter()
        .zip(copy_number.par_iter())
        .map(|(obs_cp, &total_cp)| {
            let mut est_cp: Vec<_> = obs_cp
                .iter()
                .map(|obs| (obs / cov).round().max(1f64))
                .collect();
            let sum: usize = est_cp.iter().sum::<f64>().round() as usize;
            for _ in (sum).min(total_cp)..total_cp {
                if let Some((_, est)) = obs_cp
                    .iter()
                    .zip(est_cp.iter_mut())
                    .map(|(obs, est)| {
                        let now = (obs - *est * cov).powi(2);
                        let next = (obs - (*est + 1f64) * cov).powi(2);
                        (now - next, est)
                    })
                    .max_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
                {
                    *est += 1.0;
                }
            }
            // let mut est_cp = vec![1f64; obs_cp.len()];
            // Give a copy number to the most needed one...
            // for _ in (obs_cp.len()).min(total_cp)..total_cp {
            // let (idx, _) = obs_cp
            //     .iter()
            //     .zip(est_cp.iter())
            //     .enumerate()
            //     .map(|(i, (&o, &e))| {
            //         let now = (o - e * cov).powi(2);
            //         let next = (o - (e + 1f64) * cov).powi(2);
            //         (i, now - next)
            //     })
            //     .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            //     .unwrap();
            // let (idx, _) = obs_cp
            //     .iter()
            //     .zip(est_cp.iter())
            //     .enumerate()
            //     .map(|(i, (&o, &e))| (i, 2f64 * (o - e) - 1f64))
            //     .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
            //     .unwrap();
            //     est_cp[idx] += 1f64;
            // }
            est_cp
        })
        .collect()
}

type CorrectionResult = (Vec<(u64, usize, u64)>, f64, (u64, usize));
fn correct_unit(
    ds: &DataSet,
    unit_id: u64,
    cluster_and_copynum: (usize, usize),
    copy_numbers: &[Vec<f64>],
    config: &CorrectionConfig,
) -> CorrectionResult {
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
    let unit = ds.selected_chunks.iter().find(|u| u.id == unit_id).unwrap();
    let (assignments, k) = clustering(&reads, cluster_and_copynum, copy_numbers, unit, config);
    let (adj_rand_index, raw) = adj_rand_on_biased(&reads, &assignments);
    debug!("ARI\t{}\t{k}\t{adj_rand_index:.3}\t{raw:.3}", unit.id);
    assert_eq!(assignments.len(), reads.len());
    let assignments: Vec<_> = reads
        .into_iter()
        .zip(assignments)
        .map(|((idx, read), asn)| (read.id, idx, asn as u64))
        .collect();
    (assignments, adj_rand_index, (unit_id, k))
}

fn adj_rand_on_biased(reads: &[(usize, &EncodedRead)], asns: &[usize]) -> (f64, f64) {
    const BIAS_THR: f64 = 0.2;
    assert_eq!(asns.len(), reads.len());
    let prev: Vec<_> = reads
        .iter()
        .map(|&(i, ref r)| r.nodes[i].cluster as usize)
        .collect();
    let adj_raw = crate::misc::adjusted_rand_index(&prev, &asns);
    let (prev, asns): (Vec<_>, Vec<_>) = std::iter::zip(reads, asns)
        .filter_map(|(&(i, ref read), &asn)| {
            let n = &read.nodes[i];
            n.is_biased(BIAS_THR).then(|| (n.cluster as usize, asn))
        })
        .unzip();
    let adj = crate::misc::adjusted_rand_index(&prev, &asns);
    if adj.is_nan() {
        (1f64, adj_raw)
    } else {
        (adj, adj_raw)
    }
}

type Context<'a> = (Vec<(u64, &'a [f64])>, &'a Node, Vec<(u64, &'a [f64])>);
fn to_context<'a>(&(idx, read): &(usize, &'a EncodedRead)) -> Context<'a> {
    let center = &read.nodes[idx];
    fn get_tuple(n: &Node) -> (u64, &[f64]) {
        (n.unit, n.posterior.as_slice())
    }
    let (up, tail) = match center.is_forward {
        true => {
            let up: Vec<_> = read.nodes[..idx].iter().map(get_tuple).rev().collect();
            let tail: Vec<_> = read.nodes[idx + 1..].iter().map(get_tuple).collect();
            (up, tail)
        }
        false => {
            let tail: Vec<_> = read.nodes[..idx].iter().map(get_tuple).rev().collect();
            let up: Vec<_> = read.nodes[idx + 1..].iter().map(get_tuple).collect();
            (up, tail)
        }
    };
    (up, center, tail)
}

// fn format_ctx(&(ref up, ref center, ref down): &Context) -> (Vec<String>, String, Vec<String>) {
//     fn vec_to_str(xs: &[f64]) -> String {
//         let xs: Vec<_> = xs.iter().map(|x| format!("{x:.2}")).collect();
//         xs.join(",")
//     }
//     let up: Vec<_> = up
//         .iter()
//         .map(|(u, post)| format!("[{u}({})]", vec_to_str(post)))
//         .collect();
//     let donw: Vec<_> = down
//         .iter()
//         .map(|(u, post)| format!("[{u}({})]", vec_to_str(post)))
//         .collect();
//     let center = format!("{}({})", center.unit, vec_to_str(&center.posterior));
//     (up, center, donw)
// }

fn clustering(
    reads: &[(usize, &EncodedRead)],
    (k, _upper_k): (usize, usize),
    copy_numbers: &[Vec<f64>],
    unit: &Unit,
    _: &CorrectionConfig,
) -> (Vec<usize>, usize) {
    let id = unit.id;
    let contexts: Vec<_> = reads.iter().map(to_context).collect();
    let sims: Vec<Vec<_>> = contexts
        .iter()
        .enumerate()
        .map(|(i, ctx)| {
            contexts
                .iter()
                .enumerate()
                .map(|(j, dtx)| match i == j {
                    false => alignment(ctx, dtx, copy_numbers, i, j),
                    true => 0f64,
                })
                .collect()
        })
        .collect();
    if log_enabled!(log::Level::Trace) {
        let ids: Vec<_> = reads.iter().map(|x| x.1.id).collect();
        for (i, lap) in sims.iter().enumerate() {
            let line: Vec<_> = lap.iter().map(|x| format!("{x:.2}")).collect();
            trace!("LAP\t{i}\t{}\t{}", ids[i], line.join("\t"));
        }
    }
    let cov_per_copy = reads.len() - reads.len() / unit.copy_num / 4;
    let sims = filter_similarity(sims, cov_per_copy);
    let (rowsum, laplacian) = get_graph_laplacian(&sims);
    let (mut eigens, pick_k) = get_eigenvalues(&laplacian, &rowsum, id);
    append_posterior_probability(&mut eigens, pick_k, reads);
    normalize_columns(&mut eigens);
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    let mut rng = Xoroshiro128PlusPlus::seed_from_u64(id * k as u64);
    let cluster_num = k.min(pick_k);
    let asn = (0..20)
        .map(|_| crate::misc::kmeans(&eigens, cluster_num, &mut rng))
        .min_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
        .map(|x| x.1)
        .unwrap();
    if log_enabled!(log::Level::Trace) {
        let ids: Vec<_> = reads.iter().map(|x| x.1.id).collect();
        for (i, sm) in sims.iter().enumerate() {
            let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
            trace!("SIM\t{i}\t{}\t{}", ids[i], line.join("\t"));
        }
        for (i, (up, _, tail)) in contexts.iter().enumerate() {
            trace!("LENS\t{i}\t{}\t{}", ids[i], up.len() + tail.len());
        }
        for (i, sm) in eigens.iter().enumerate() {
            let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
            let line = line.join("\t");
            let cls = contexts[i].1.cluster;
            let len = reads[i].1.nodes.len();
            let asn = asn[i];
            trace!("EIG\t{id}\t{i}\t{cls}\t{len}\t{asn}\t{line}",);
        }
    }
    assert!(asn.iter().all(|&x| x <= cluster_num));
    (asn, cluster_num)
}

fn filter_similarity(mut sims: Vec<Vec<f64>>, len: usize) -> Vec<Vec<f64>> {
    const SMALL: f64 = 0.0000000000000001;
    const MIN_REQ: f64 = 0.51;
    let mut to_retain: Vec<Vec<_>> = sims.iter().map(|xs| vec![false; xs.len()]).collect();
    for (i, xs) in sims.iter().enumerate() {
        let threshold = select_nth(xs, len).max(MIN_REQ);
        for (j, _) in xs.iter().enumerate().filter(|x| threshold <= *x.1) {
            to_retain[i][j] = true;
            to_retain[j][i] = true;
        }
    }
    for (xs, retains) in sims.iter_mut().zip(to_retain) {
        for (x, _) in xs.iter_mut().zip(retains).filter(|x| !x.1) {
            *x = SMALL;
        }
    }
    sims
}

fn select_nth(sims: &[f64], pivot: usize) -> f64 {
    assert!(pivot <= sims.len());
    let mut sims = sims.to_vec();
    sims.sort_by(|x, y| x.partial_cmp(y).unwrap());
    sims[pivot]
}

fn append_posterior_probability(
    eigens: &mut [Vec<f64>],
    k: usize,
    reads: &[(usize, &EncodedRead)],
) {
    for (eigen, (idx, read)) in eigens.iter_mut().zip(reads.iter()) {
        assert_eq!(eigen.len(), k);
        let post = &read.nodes[*idx].posterior;
        let total = crate::misc::logsumexp(post);
        eigen.extend(post.iter().map(|x| (x - total).exp()));
    }
}

fn normalize_columns(eigens: &mut [Vec<f64>]) {
    let mut sums = vec![0f64; eigens[0].len()];
    for eig in eigens.iter() {
        sums.iter_mut()
            .zip(eig.iter())
            .for_each(|(s, e)| *s += e * e);
    }
    sums.iter_mut().for_each(|x| *x = x.sqrt());
    for eigen in eigens.iter_mut() {
        let posts = eigen.iter_mut().zip(sums.iter());
        posts.for_each(|(e, s)| *e /= s);
    }
}

// Return the *normalized* graph lap, which is I - D^{-1/2}WD^{-1/2}
// also return D.
fn get_graph_laplacian(sims: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let rowsum: Vec<f64> = sims.iter().map(|xs| xs.iter().sum()).collect();
    let sq_inv: Vec<_> = rowsum.iter().map(|xs| xs.recip().sqrt()).collect();
    let lap: Vec<Vec<_>> = sims
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, w)| match j == i {
                    true => 1f64,
                    false => -w * sq_inv[i] * sq_inv[j],
                })
                .collect()
        })
        .collect();
    (rowsum, lap)
}

// Return the *normalized* graph lap.
// I - D^{-1}W, where I is identity matrix, D is the diagonal "summed" matrix and
// W is the similarity matrix itself..
// fn get_graph_laplacian(sims: &[Vec<f64>]) -> Vec<Vec<f64>> {
//     sims.iter()
//         .enumerate()
//         .map(|(i, row)| {
//             let sum: f64 = row.iter().sum();
//             assert!(0.000001 < sum);
//             let div: f64 = row.iter().sum::<f64>().recip();
//             row.iter()
//                 .enumerate()
//                 .map(|(t, w)| if t == i { 1f64 } else { -w * div })
//                 .collect()
//         })
//         .collect()
// }

// const EIGEN_THR: f64 = 0.25;
const EIGEN_THR: f64 = 0.2;
fn get_eigenvalues(matrix: &[Vec<f64>], rowsum: &[f64], id: u64) -> (Vec<Vec<f64>>, usize) {
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
    if log_enabled!(log::Level::Trace) {
        let eigs: Vec<_> = eigen_and_eigenvec
            .iter()
            .map(|x| format!("{:.3}", x.1))
            .take(opt_k + 3)
            .collect();
        trace!("EIGVALUE\t{id}\t{}", eigs.join("\t"));
    }
    if pick_k == 0 {
        for row in matrix.row_iter() {
            let sum: f64 = row.iter().sum();
            let row: Vec<_> = row.iter().map(|x| format!("{:.2}", x)).collect();
            error!("LAP\t{sum:.1}\t{}", row.join("\t"));
        }
        let eigens: Vec<_> = eigen_and_eigenvec.iter().take(4).map(|x| x.1).collect();
        panic!("{:?},{},{}", eigens, id, opt_k);
    }
    let sqrt_inv = rowsum.iter().map(|x| x.recip().sqrt());
    let features: Vec<Vec<_>> = (0..datalen)
        .zip(sqrt_inv)
        .map(|(i, d)| {
            (0..pick_k)
                .map(|j| eigen_and_eigenvec[j].0[i] * d)
                .collect()
        })
        .collect();
    (features, opt_k)
}

fn alignment<'a>(
    (up1, center1, down1): &Context<'a>,
    ctx2: &Context<'a>,
    copy_numbers: &[Vec<f64>],
    _i: usize,
    _j: usize,
) -> f64 {
    let (up2, center2, down2) = ctx2;
    assert_eq!(center1.unit, center2.unit);
    let up_aln = align_swg(up1, up2, copy_numbers);
    let down_aln = align_swg(down1, down2, copy_numbers);
    let center_copy_num = &copy_numbers[center1.unit as usize];
    let center = sim(&center1.posterior, &center2.posterior, center_copy_num);
    let likelihood_ratio = up_aln + down_aln + center;
    (1f64 + (-likelihood_ratio).exp()).recip()
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
    let logp = crate::misc::logsumexp_str(iter);
    let logit = logit_from_lnp(logp);
    assert!(!logit.is_infinite(), "{},{}", logp, logit);
    logit
}

// log(p) -> log(p) - log(1-p)
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
