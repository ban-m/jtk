use definitions::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
const SMALL_WEIGHT: f64 = 0.000_000_001;
mod config;
pub mod naive_variant_calling;
pub use config::*;
pub mod create_model;
pub mod eread;
pub mod variant_calling;
use create_model::*;
use eread::*;
use variant_calling::*;
pub mod kmeans;
pub trait LocalClustering {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        self,
        c: &ClusteringConfig<F>,
    ) -> Self;
}

impl LocalClustering for DataSet {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        mut self,
        c: &ClusteringConfig<F>,
    ) -> Self {
        let max_unit = self.selected_chunks.iter().map(|u| u.id).max().unwrap() + 1;
        let mut pileups: Vec<Vec<(_, _, &mut Node)>> = (0..max_unit).map(|_| Vec::new()).collect();
        for read in self.encoded_reads.iter_mut() {
            for (idx, node) in read.nodes.iter_mut().enumerate() {
                pileups[node.unit as usize].push((read.id, idx, node));
            }
        }
        let selected_chunks = self.selected_chunks.clone();
        let id_to_name: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
        let id_to_desc: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
        pileups
            .par_iter_mut()
            .enumerate()
            .for_each(|(unit_id, mut units)| {
                debug!("The {}-th pileup", unit_id);
                let ref_unit = match selected_chunks.iter().find(|u| u.id == unit_id as u64) {
                    Some(res) => res,
                    None => return,
                };
                assert!(units.iter().all(|n| n.2.unit == unit_id as u64));
                unit_clustering(&mut units, c, ref_unit);
                if log_enabled!(log::Level::Debug) {
                    for cl in 0..c.cluster_num {
                        let cl = cl as u64;
                        for (id, _, _) in units.iter().filter(|&(_, _, n)| n.cluster == cl) {
                            let name = id_to_name[&id];
                            match id_to_desc.get(&id) {
                                Some(d) => debug!("{}\t{}\t{}\t{}\t{}", unit_id, cl, id, name, d),
                                None => debug!("{}\t{}\t{}\t{}", unit_id, cl, id, name),
                            }
                        }
                    }
                }
            });
        self
    }
}

pub fn unit_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    units: &mut [(u64, usize, &mut Node)],
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
) {
    let len = if ref_unit.seq().len() % c.subchunk_length == 0 {
        ref_unit.seq().len() / c.subchunk_length
    } else {
        ref_unit.seq().len() / c.subchunk_length + 1
    };
    let (min_length, max_length) = (c.subchunk_length * 9 / 10, c.subchunk_length * 11 / 10);
    let mut data: Vec<_> = units
        .iter()
        .map(|&(_, _, ref node)| {
            let mut chunks = node_to_subchunks(node, c.subchunk_length);
            chunks.retain(|chunk| min_length < chunk.seq.len() && chunk.seq.len() < max_length);
            ChunkedUnit { cluster: 0, chunks }
        })
        .collect();
    let _m = clustering_by_kmeans(&mut data, len, &c, ref_unit.id);
    for ((_, _, ref mut n), d) in units.iter_mut().zip(data.iter()) {
        n.cluster = d.cluster as u64;
    }
}

pub fn node_to_subchunks(node: &Node, len: usize) -> Vec<Chunk> {
    let (mut q_pos, mut r_pos) = (0, 0);
    let mut target = len;
    let mut ops = node.cigar.clone();
    ops.reverse();
    let mut chunk_position = vec![0];
    let mut scores = vec![0];
    loop {
        // Pop until target.
        let mut score = 0;
        let last_op = {
            loop {
                match ops.pop() {
                    Some(Op::Del(l)) => {
                        r_pos += l;
                        score += -2 * l as i64;
                        if target <= r_pos {
                            break Op::Del(l);
                        }
                    }
                    Some(Op::Ins(l)) => {
                        score += -2 * l as i64;
                        q_pos += l;
                    }
                    Some(Op::Match(l)) => {
                        r_pos += l;
                        q_pos += l;
                        // This confuses match and mismatch.
                        // Should be tuned.
                        score += l as i64;
                        if target <= r_pos {
                            break Op::Match(l);
                        }
                    }
                    None => unreachable!(),
                }
            }
        };
        // Push back
        if target < r_pos {
            let overflow = r_pos - target;
            match last_op {
                Op::Del(_) => {
                    ops.push(Op::Del(overflow));
                    r_pos -= overflow;
                    score += 2 * overflow as i64;
                }
                Op::Match(_) => {
                    ops.push(Op::Match(overflow));
                    r_pos -= overflow;
                    q_pos -= overflow;
                    score -= overflow as i64;
                }
                _ => unreachable!(),
            }
        }
        chunk_position.push(q_pos);
        scores.push(score);
        // Move to next iteration
        let rest = ops
            .iter()
            .map(|op| match op {
                Op::Del(l) | Op::Match(l) => *l,
                _ => 0,
            })
            .sum::<usize>();
        let query_rest = ops
            .iter()
            .map(|op| match op {
                Op::Ins(l) | Op::Match(l) => *l,
                _ => 0,
            })
            .sum::<usize>();
        if rest < len && query_rest > 0 {
            let score = ops
                .iter()
                .map(|op| match op {
                    Op::Ins(l) | Op::Del(l) => -2 * *l as i64,
                    Op::Match(l) => *l as i64,
                })
                .sum::<i64>();
            scores.push(score);
            chunk_position.push(node.seq().len());
            break;
        } else if rest < len && query_rest == 0 {
            break;
        } else {
            target += len;
        }
    }
    assert_eq!(scores.len(), chunk_position.len());
    chunk_position
        .windows(2)
        .zip(scores.windows(2))
        .enumerate()
        .filter(|(_, (_, s))| s[1] > 0)
        .map(|(idx, (w, _))| Chunk {
            pos: idx,
            seq: node.seq()[w[0]..w[1]].to_vec(),
        })
        .collect()
}

fn update_by_precomputed_lks<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &mut [ChunkedUnit],
    chain_len: usize,
    beta: f64,
    pos: &[bool],
    lks: &[Vec<f64>],
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> u32 {
    let fractions: Vec<_> = (0..c.cluster_num)
        .map(|cl| 1 + data.iter().filter(|d| d.cluster == cl).count())
        .collect();
    let choices: Vec<_> = (0..c.cluster_num).collect();
    data.iter_mut()
        .enumerate()
        .map(|(idx, d)| {
            let mut weights: Vec<_> = (0..c.cluster_num)
                .map(|l| {
                    (0..c.cluster_num)
                        .map(|k| {
                            if l == k {
                                0.
                            } else {
                                let frac = (fractions[k] as f64 / fractions[l] as f64).ln();
                                let modellk = (0..chain_len)
                                    .filter(|&p| pos[p])
                                    .map(|p| {
                                        lks[idx][p + k * chain_len] - lks[idx][p + l * chain_len]
                                    })
                                    .sum::<f64>();
                                modellk + beta * frac
                            }
                        })
                        .map(|lkdiff| lkdiff.exp())
                        .sum::<f64>()
                        .recip()
                })
                .collect();
            let sum = weights.iter().sum::<f64>();
            weights.iter_mut().for_each(|x| *x /= sum);
            let (new_assignment, max) = weights
                .iter()
                .enumerate()
                .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or((0, &0.));
            use rand::seq::SliceRandom;
            let new_assignment = *choices.choose_weighted(rng, |&k| weights[k]).unwrap();
            let diff = weights
                .iter()
                .filter(|&x| (x - max).abs() > 0.001)
                .map(|x| max - x)
                .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .unwrap_or(0.);
            // trace!("{}\t{:.3}\t{:.3}\t{:.3}", idx, weights[0], weights[1], diff);
            if diff > 0.3 && d.cluster != new_assignment {
                d.cluster = new_assignment;
                1
            } else {
                0
            }
        })
        .sum::<u32>()
}

pub fn clustering_by_kmeans<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    ref_unit: u64,
) -> f64 {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(ref_unit);
    let mut count = 0;
    let rng = &mut rng;
    let start = std::time::Instant::now();
    let mut coef = 1.;
    let stable_thr = match c.read_type {
        ReadType::CCS => (data.len() as f64 * 0.02).max(1.).floor() as u32,
        ReadType::ONT => (data.len() as f64 * 0.08).max(3.).floor() as u32,
        ReadType::CLR => (data.len() as f64 * 0.08).max(3.).floor() as u32,
    };
    let (betas, pos, lks, margin) = (0..20)
        .map(|_| {
            data.iter_mut().for_each(|cs| {
                cs.cluster = rng.gen_range(0, c.cluster_num);
            });
            get_variants(&data, chain_len, rng, c)
        })
        .max_by(|x, y| (x.3).partial_cmp(&y.3).unwrap())
        .unwrap();
    trace!("Init Margin:{}", margin);
    let beta = 1.;
    // let beta = (c.initial_beta * coef).min(c.max_beta);
    for (idx, d) in data.iter().enumerate() {
        trace!("{}\t{}", idx, d.cluster);
    }
    let lks = average_correction(lks, chain_len, c.cluster_num);
    update_by_precomputed_lks(data, chain_len, beta, &pos, &lks, rng, c);
    loop {
        // coef *= c.beta_increase;
        // let beta = (c.initial_beta * coef).min(c.max_beta);
        let (betas, pos, lks, margin) = get_variants(&data, chain_len, rng, c);
        let lks = average_correction(lks, chain_len, c.cluster_num);
        let changed_num = update_by_precomputed_lks(data, chain_len, beta, &pos, &lks, rng, c);
        count += (changed_num <= stable_thr) as u32;
        count *= (changed_num <= stable_thr) as u32;
        report(ref_unit, &data, count, margin, coef, changed_num);
        if (std::time::Instant::now() - start).as_secs() > c.limit {
            info!("Break by timelimit:{}", c.limit);
            break margin;
        } else if count >= c.stable_limit {
            break margin;
        }
    }
}

fn average_correction(
    mut lks: Vec<Vec<f64>>,
    chain_len: usize,
    cluster_num: usize,
) -> Vec<Vec<f64>> {
    let mut count = vec![0.; chain_len * cluster_num];
    // To avoid zero-division.
    let mut sum = vec![0.001; chain_len * cluster_num];
    for read in lks.iter() {
        for (i, lk) in read.iter().enumerate() {
            if lk.is_sign_negative() {
                count[i] += 1.;
                sum[i] += lk;
            }
        }
    }
    let mean_vec: Vec<_> = sum.iter().zip(count).map(|(&s, c)| s / c).collect();
    for matrix in lks.iter_mut() {
        assert_eq!(mean_vec.len(), matrix.len());
        matrix
            .iter_mut()
            .zip(mean_vec.iter())
            .for_each(|(x, &y)| *x -= y);
    }
    lks
}

fn report(id: u64, data: &[ChunkedUnit], count: u32, lk: f64, beta: f64, c: u32) {
    if !log_enabled!(log::Level::Trace) {
        return;
    }
    let counts: Vec<_> = {
        let mut count: std::collections::HashMap<usize, usize> = HashMap::new();
        for read in data {
            *count.entry(read.cluster).or_default() += 1;
        }
        let mut count: Vec<(usize, usize)> = count.into_iter().collect();
        count.sort_by_key(|e| e.0);
        count.iter().map(|(_, cnt)| format!("{}", cnt)).collect()
    };
    let counts = counts.join("\t");
    trace!(
        "Summary\t{}\t{}\t{:.3}\t{:.3}\t{}\t{}",
        id,
        count,
        lk,
        beta,
        c,
        counts
    );
}

fn update_assignment<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &mut [ChunkedUnit],
    c: &ClusteringConfig<F>,
    picked: usize,
    betas: &[Vec<Vec<f64>>],
    beta: f64,
    pos: &[bool],
    chain_len: usize,
    rng: &mut R,
) -> u32 {
    let mut weights = get_weights(data, chain_len, rng, &pos, picked, beta, betas, c);
    let sum = weights.iter().sum::<f64>();
    weights.iter_mut().for_each(|x| *x /= sum);
    let (new_assignment, max) = weights
        .iter()
        .enumerate()
        .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((0, &0.));
    let diff = weights
        .iter()
        .filter(|&x| (x - max).abs() > 0.001)
        .map(|x| max - x)
        .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(0.);
    trace!(
        "{}\t{:.3}\t{:.3}\t{:.3}",
        picked,
        weights[0],
        weights[1],
        diff,
    );
    let read = data.get_mut(picked).unwrap();
    if diff > 0.4 && read.cluster != new_assignment {
        read.cluster = new_assignment;
        1
    } else {
        0
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

fn get_weights<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: rand::Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    pos: &[bool],
    picked: usize,
    beta: f64,
    betas: &[Vec<Vec<f64>>],
    c: &ClusteringConfig<F>,
) -> Vec<f64> {
    let fractions: Vec<_> = (0..c.cluster_num)
        .map(|cl| data.iter().filter(|d| d.cluster == cl).count() + 1)
        .collect();
    let models: Vec<_> = get_models(&data, chain_len, rng, c, &pos, Some(picked));
    let mut likelihoods: Vec<(usize, Vec<_>)> = data[picked]
        .chunks
        .iter()
        .filter(|chunk| pos[chunk.pos])
        .map(|chunk| {
            let lks: Vec<_> = models
                .iter()
                .map(|ms| vec![ms[chunk.pos].forward(&chunk.seq, &c.poa_config)])
                .collect();
            (chunk.pos, lks)
        })
        .collect();
    for _ in 0..c.repeat_num - 1 {
        let models: Vec<_> = get_models(&data, chain_len, rng, c, &pos, Some(picked));
        for (chunk, (pos, lks)) in data[picked]
            .chunks
            .iter()
            .filter(|chunk| pos[chunk.pos])
            .zip(likelihoods.iter_mut())
        {
            assert_eq!(*pos, chunk.pos);
            lks.iter_mut()
                .zip(models.iter())
                .for_each(|(xs, ms)| xs.push(ms[chunk.pos].forward(&chunk.seq, &c.poa_config)));
        }
    }
    let likelihoods: Vec<(usize, Vec<_>)> = likelihoods
        .into_iter()
        .map(|(pos, lks)| {
            let lks: Vec<_> = lks
                .iter()
                .map(|lks| logsumexp(lks) - (c.repeat_num as f64).ln())
                .collect();
            // let dump: Vec<_> = lks
            //     .iter()
            //     .map(|x| lks.iter().map(|y| (y - x).exp()).sum::<f64>().recip())
            //     .map(|x| format!("{}", x))
            //     .collect();
            // debug!("DUMP\t{}\t{}", pos, dump.join("\t"));
            (pos, lks)
        })
        .collect();
    (0..c.cluster_num)
        .map(|l| {
            (0..c.cluster_num)
                .map(|k| {
                    if k == l {
                        return 0.;
                    }
                    let (i, j) = (l.max(k), l.min(k));
                    let prior = (fractions[k] as f64 / fractions[l] as f64).ln();
                    beta * prior
                        + likelihoods
                            .iter()
                            .map(|&(pos, ref lks)| betas[i][j][pos] * (lks[k] - lks[l]))
                            .sum::<f64>()
                })
                .map(|lkdiff| lkdiff.exp())
                .sum::<f64>()
                .recip()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #[derive(Clone, Copy, Debug)]
    struct TestConfig {
        cl: usize,
        num: usize,
        fail: f64,
        skip: f64,
        max_len: usize,
        min_len: usize,
        unit_len: usize,
    }
    use super::*;
    #[test]
    fn works() {}
    #[allow(dead_code)]
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> Vec<ERead> {
        let TestConfig {
            cl,
            num,
            fail,
            skip,
            max_len,
            min_len,
            unit_len,
        } = conf;
        use rand::seq::SliceRandom;
        let _unit_map: HashMap<usize, u64> = {
            let units: Vec<_> = (0..unit_len).collect();
            units
                .choose_multiple(r, unit_len)
                .enumerate()
                .map(|(idx, &u)| (idx, u as u64))
                .collect()
        };
        (0..num)
            .map(|i| {
                let id = i as u64;
                let cluster = i % cl;
                let len = r.gen::<usize>() % (max_len - min_len) + min_len;
                let path: Vec<_> = if r.gen_bool(0.5) {
                    let start = r.gen::<usize>() % (unit_len - len);
                    (start..start + len)
                        .map(|unit| unit as u64)
                        .filter_map(|unit| {
                            if r.gen_bool(skip) {
                                None
                            } else if r.gen_bool(fail) {
                                let cluster = r.gen::<usize>() % cl;
                                Some(Elm { unit, cluster })
                            } else {
                                Some(Elm { unit, cluster })
                            }
                        })
                        .collect()
                } else {
                    let start = r.gen::<usize>() % (unit_len - len) + len;
                    (start - len..start)
                        .rev()
                        .map(|unit| unit as u64)
                        .filter_map(|unit| {
                            if r.gen_bool(skip) {
                                None
                            } else if r.gen_bool(fail) {
                                let cluster = r.gen::<usize>() % cl;
                                Some(Elm { unit, cluster })
                            } else {
                                Some(Elm { unit, cluster })
                            }
                        })
                        .collect()
                };
                ERead { id, cluster, path }
            })
            .collect()
    }
}
