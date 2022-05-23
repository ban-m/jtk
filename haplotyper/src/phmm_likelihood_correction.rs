use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Default)]
pub struct CorrectionConfig {}
pub trait AlignmentCorrection {
    fn correct_clustering(&mut self, config: &CorrectionConfig);
    // fn correct_clustering_selected(&mut self, selection: &HashSet<u64>, config: &CorrectionConfig);
}
impl AlignmentCorrection for DataSet {
    // TODO: After calling this function, correspondance between
    // the posterior distributions and the assignments would be broken.
    // shoule we fix
    fn correct_clustering(&mut self, config: &CorrectionConfig) {
        let selections: Vec<_> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .filter(|&(_, k)| 1 < k)
            // .filter(|&(id, _)| id == 659)
            .collect();
        let copy_numbers = estimate_copy_number_of_cluster(self);
        // let posterior_distributions: Vec<_> = selections
        let corrected_clusterings: Vec<_> = selections
            .par_iter()
            .map(|&(id, cluster_num)| correct_unit(self, id, cluster_num, &copy_numbers, &config))
            .collect();
        let corrected_clustering_on_read = {
            let mut corrected_clustering_on_read: HashMap<_, Vec<_>> = HashMap::new();
            for correcteds in corrected_clusterings {
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
        // let mut result: HashMap<u64, Vec<(usize, Vec<f64>)>> = {
        //     let mut result: HashMap<u64, Vec<(usize, Vec<f64>)>> = HashMap::new();
        //     for post_dist_on_chunk in posterior_distributions {
        //         for (id, pos, posterior_dist) in post_dist_on_chunk {
        //             result.entry(id).or_default().push((pos, posterior_dist));
        //         }
        //     }
        //     result
        // };
        // for read in self.encoded_reads.iter_mut() {
        //     if let Some(corrected) = result.remove(&read.id) {
        // for (pos, post) in corrected {
        //     read.nodes[pos].cluster = argmax(&post) as u64;
        // }
        //     }
        // }
    }
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

fn correct_unit(
    ds: &DataSet,
    unit_id: u64,
    k: usize,
    copy_numbers: &[Vec<f64>],
    config: &CorrectionConfig,
) -> Vec<(u64, usize, u64)> {
    //) -> Vec<(u64, usize, Vec<f64>)> {
    let mut reads = vec![];
    for read in ds.encoded_reads.iter() {
        for (idx, node) in read.nodes.iter().enumerate() {
            if node.unit == unit_id {
                reads.push((idx, read));
            }
        }
    }
    // let id2desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    // reads.sort_by_cached_key(|&(idx, ref read)| {
    //     let answer = id2desc[&read.id].contains("000251v2") as usize;
    //     let cluster = read.nodes[idx].cluster;
    //     (answer, cluster)
    // });
    // for (i, &(idx, ref read)) in reads.iter().enumerate() {
    //     let ans = id2desc[&read.id].contains("000251v2") as usize;
    //     let cls = read.nodes[idx].cluster;
    //     debug!("ANSWER\t{i}\t{cls}\t{ans}");
    // }
    // debug!("CENTER\t{:?}", copy_numbers[unit_id as usize]);
    // debug!("Correction\t{unit_id}\t{k}\t{}", reads.len());
    let assignments = clustering(&reads, k, copy_numbers, config);
    assert_eq!(assignments.len(), reads.len());
    reads
        .into_iter()
        .zip(assignments)
        .map(|((idx, read), asn)| (read.id, idx, asn as u64))
        .collect()
    // let mut clusters: Vec<Vec<_>> = vec![vec![]; k];
    // for (&asn, elm) in assignments.iter().zip(reads) {
    //     clusters[asn].push(elm);
    // }
    // clusters
    //     .into_iter()
    //     .enumerate()
    //     .flat_map(|(i, cluster)| {
    //         // Here I want to compute p1p2...pn for each cluster,
    //         // then normalize them.ln(p1p2...pn/norm) -> (lnp1 + lnp2 ... ) - ln(sum(exp(lnp1+...)))
    //         let mut summed: Vec<_> =
    //             cluster
    //                 .iter()
    //                 .fold(vec![0f64; k], |mut sums, &(idx, read)| {
    //                     let node = &read.nodes[idx];
    //                     assert_eq!(node.posterior.len(), sums.len());
    //                     sums.iter_mut()
    //                         .zip(node.posterior.iter())
    //                         .for_each(|(x, y)| *x += y);
    //                     sums
    //                 });
    //         let total = logsumexp(&summed);
    //         summed.iter_mut().for_each(|x| *x = *x - total);
    //         let line: Vec<_> = summed.iter().map(|x| format!("{x:.2}")).collect();
    //         debug!("CLU\t{i}\t{}", line.join("\t"));
    //         cluster
    //             .iter()
    //             .map(|&(idx, read)| (read.id, idx, summed.clone()))
    //             .collect::<Vec<_>>()
    //     })
    //     .collect()
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
    k: usize,
    copy_numbers: &[Vec<f64>],
    _: &CorrectionConfig,
) -> Vec<usize> {
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
    //let ids: Vec<_> = reads.iter().map(|x| x.1.id).collect();
    let sims: Vec<Vec<_>> = contexts
        .iter()
        .enumerate()
        .map(|(_i, ctx)| {
            // let id = ids[i];
            contexts
                .iter()
                .enumerate()
                .map(|(_j, dtx)| {
                    // debug!("-------------");
                    // debug!("ALN\t{id}\t{i}\t{j}");
                    alignment(ctx, dtx, copy_numbers) + 0.00000001
                })
                .collect()
        })
        .collect();
    // for (i, sm) in sims.iter().enumerate() {
    //     let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
    //     debug!("SIM\t{i}\t{}\t{}", ids[i], line.join("\t"));
    // }
    let laplacian = get_graph_laplacian(&sims);
    // for (i, sm) in laplacian.iter().enumerate() {
    //     let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
    //     debug!("LAP\t{i}\t{}\t{}", ids[i], line.join("\t"));
    // }
    let eigens = get_eigenvalues(&laplacian, k);
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    let mut rng = Xoroshiro128PlusPlus::seed_from_u64(394203);
    let asn = (0..10)
        .map(|_| kmeans(&eigens, k, &mut rng))
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .map(|x| x.0)
        .unwrap();
    // for (i, sm) in eigens.iter().enumerate() {
    //     let line: Vec<_> = sm.iter().map(|x| format!("{x:.2}")).collect();
    //     let cls = contexts[i].1.cluster;
    //     debug!("EIG\t{i}\t{}\t{}\t{}", cls, asn[i], line.join("\t"));
    // }
    asn
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

fn get_eigenvalues(matrix: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let datalen = matrix.len();
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
    // eigen_and_eigenvec.reverse();
    let top_k_eigenvec = eigen_and_eigenvec[..k]
        .iter()
        .flat_map(|(v, _)| v.iter().copied());
    let top_k_eigenvec = nalgebra::DMatrix::from_iterator(datalen, k, top_k_eigenvec);
    top_k_eigenvec
        .row_iter()
        .map(|row| row.iter().map(|&x| x).collect())
        .collect()
    // let pca_vectors = &eigen_and_eigenvec[..k];
    // matrix
    //     .column_iter()
    //     .map(|x| pca_vectors.iter().map(|(v, _)| x.dot(v)).collect())
    //     .collect()
}

use rand::Rng;
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
fn alignment<'a>(
    (up1, center1, down1): &Context<'a>,
    (up2, center2, down2): &Context<'a>,
    copy_numbers: &[Vec<f64>],
) -> f64 {
    assert_eq!(center1.unit, center2.unit);
    let center_copy_num = &copy_numbers[center1.unit as usize];
    let center = sim(&center1.posterior, &center2.posterior, center_copy_num);
    let up_aln = align(up1, up2, &copy_numbers);
    let down_aln = align(down1, down2, &copy_numbers);
    // let (up_aln, _up_len) = align(up1, up2, &copy_numbers);
    // let (down_aln, _down_len) = align(down1, down2, &copy_numbers);
    // let similarity = ((up_aln + down_aln + center) / (up_len + down_len + 1) as f64).exp();
    //let similarity = (up_aln + down_aln + center).exp();
    let similarity = (up_aln + down_aln + center).max(0f64);
    assert!(0f64 <= similarity);
    // let line1: Vec<_> = up1
    //     .iter()
    //     .map(|(u, p)| (*u, argmax(p)))
    //     .chain(std::iter::once((center1.unit, center1.cluster as usize)))
    //     .chain(down1.iter().map(|(u, p)| (*u, argmax(p))))
    //     .map(|(u, c)| format!("{u}-{c}"))
    //     .collect();
    // let line2: Vec<_> = up2
    //     .iter()
    //     .map(|(u, p)| (*u, argmax(p)))
    //     .chain(std::iter::once((center2.unit, center2.cluster as usize)))
    //     .chain(down2.iter().map(|(u, p)| (*u, argmax(p))))
    //     .map(|(u, c)| format!("{u}-{c}"))
    //     .collect();
    // debug!("ALN\t{}", line1.join("\t"));
    // debug!("ALN\t{}", line2.join("\t"));
    // debug!("ALN\t{up_aln:.3}\t{center:.3}\t{down_aln:.3}\t{similarity:.3}");
    similarity
}

fn align<'a>(
    arm1: &[(u64, &'a [f64])],
    arm2: &[(u64, &'a [f64])],
    copy_numbers: &[Vec<f64>],
) -> f64 {
    // const GAP: f64 = -1f64;
    // Allow gap without any penalty.
    // It might seem very peculer, but it is OK,
    // because if this deletion is indeed true one,
    // other (unit,cluster) pair would be informative.
    // Otherwise, or the deletion is just random error,
    // we should ignore them from the score calulation,
    // as we do usual seuqnecing errors.
    const GAP: f64 = 0f64;
    const MISM: f64 = -10000f64;
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
    let last_row = dp.last().unwrap().iter();
    let last_column = dp.iter().filter_map(|x| x.last());
    *last_row
        .chain(last_column)
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap()
    // let last_row = dp
    //     .last()
    //     .unwrap()
    //     .iter()
    //     .enumerate()
    //     .map(|(j, s)| (arm1.len(), j, s));
    // let last_column = dp
    //     .iter()
    //     .filter_map(|x| x.last())
    //     .enumerate()
    //     .map(|(i, s)| (i, arm2.len(), s));
    // let (mut i, mut j, score) = last_row
    //     .chain(last_column)
    //     .max_by(|x, y| x.2.partial_cmp(&y.2).unwrap())
    //     .unwrap();
    // let mut alnlen = 0;
    // while 0 < i && 0 < j {
    //     alnlen += 1;
    //     let current = dp[i][j];
    //     let (u1, p1) = arm1[i - 1];
    //     let (u2, p2) = arm2[j - 1];
    //     let match_score = match u1 == u2 {
    //         true => sim(p1, p2, &copy_numbers[u1 as usize]),
    //         false => MISM,
    //     };
    //     assert!(!match_score.is_nan());
    //     if (dp[i - 1][j - 1] + match_score - current).abs() < 0.00001 {
    //         i -= 1;
    //         j -= 1;
    //     } else if (dp[i - 1][j] + GAP - current).abs() < 0.000001 {
    //         i -= 1;
    //     } else {
    //         assert!((dp[i][j - 1] + GAP - current).abs() < 0.000001);
    //         j -= 1;
    //     }
    // }
    // alnlen += j + i;
    // (*score, alnlen)
}

const MOCK_CP: f64 = 1.5;
pub fn sim(xs: &[f64], ys: &[f64], cps: &[f64]) -> f64 {
    assert_eq!(xs.len(), cps.len());
    assert_eq!(xs.len(), ys.len());
    // TODO: .ln before call this function!
    if cps.len() == 1 {
        let total: f64 = cps.iter().sum();
        return (total.max(MOCK_CP) - 1f64).recip().ln();
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
