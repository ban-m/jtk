use definitions::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
mod config;
pub use config::*;
pub mod create_model;
pub mod eread;
pub mod variant_calling;
const REPEAT_NUM: usize = 8;

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
        // If true, the unit is clustered.
        let mut clustered_units: HashMap<_, _> =
            self.selected_chunks.iter().map(|u| (u.id, false)).collect();
        let mut unit_num = clustered_units.values().filter(|&&x| !x).count();
        for i in 0..c.retry_limit {
            debug!("NUM\t{}\t{}\tUnits to be clsutered.", i, unit_num);
            debug!("==================================");
            local_clustering_on_selected(&mut self, c, &clustered_units, i);
            clustered_units = super::unit_correlation::select_uninformative_units(&self, c.p_value);
            let next_num = clustered_units.values().filter(|&&x| !x).count();
            let diff = unit_num.max(next_num) - unit_num.min(next_num);
            unit_num = next_num;
            if next_num == 0 || diff == 0 {
                debug!("DONE\t{}\t{}", next_num, diff);
                break;
            }
        }
        debug!("Remaining {} unresolved cluster.", unit_num);
        self
    }
}

fn local_clustering_on_selected<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    ds: &mut DataSet,
    c: &ClusteringConfig<F>,
    positions: &HashMap<u64, bool>,
    iteration: u64,
) {
    // If positions[idx] is *false*, the position would be used.
    let mut pileups: HashMap<u64, Vec<(_, _, &mut Node)>> = HashMap::new();
    for read in ds.encoded_reads.iter_mut() {
        for (idx, node) in read.nodes.iter_mut().enumerate() {
            if let Some(&res) = positions.get(&node.unit) {
                // This is correct.
                if !res {
                    pileups
                        .entry(node.unit)
                        .or_default()
                        .push((read.id, idx, node));
                }
            }
        }
    }
    let selected_chunks = ds.selected_chunks.clone();
    let id_to_name: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
    let id_to_desc: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
    pileups.par_iter_mut().for_each(|(&unit_id, mut units)| {
        debug!("RECORD\t{}\t{}\tStart", unit_id, iteration);
        let ref_unit = match selected_chunks.iter().find(|u| u.id == unit_id as u64) {
            Some(res) => res,
            None => return,
        };
        // If true, this position is already clustered.
        if positions[&unit_id] {
            return;
        }
        assert!(units.iter().all(|n| n.2.unit == unit_id as u64));
        unit_clustering(&mut units, c, ref_unit, iteration);
        debug!("RECORD\t{}\t{}\tEnd", unit_id, iteration);
        if log_enabled!(log::Level::Debug) {
            for cl in 0..c.cluster_num.max(ref_unit.cluster_num) {
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
}

pub fn unit_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    units: &mut [(u64, usize, &mut Node)],
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
    iteration: u64,
) {
    if ref_unit.cluster_num == 1 {
        units.iter_mut().for_each(|x| x.2.cluster = 0);
        return;
    }
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
            let cluster = if c.retain_current_clustering {
                node.cluster as usize
            } else {
                0
            };
            ChunkedUnit { cluster, chunks }
        })
        .collect();
    let seed = ref_unit.id * iteration;
    clustering_by_kmeans(&mut data, len, &c, ref_unit, seed);
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

fn update_by_precomputed_lks<R: Rng>(
    data: &mut [ChunkedUnit],
    (cluster_num, chain_len): (usize, usize),
    pick: f64,
    pos: &[bool],
    lks: &[Vec<f64>],
    rng: &mut R,
) -> u32 {
    let fractions: Vec<_> = (0..cluster_num)
        .map(|cl| 1 + data.iter().filter(|d| d.cluster == cl).count())
        .collect();
    let choices: Vec<_> = (0..cluster_num).collect();
    data.iter_mut()
        .enumerate()
        .map(|(idx, d)| {
            let mut weights: Vec<_> = (0..cluster_num)
                .map(|l| {
                    (0..cluster_num)
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
                                modellk + frac
                            }
                        })
                        .map(|lkdiff| lkdiff.exp())
                        .sum::<f64>()
                        .recip()
                })
                .collect();
            let sum = weights.iter().sum::<f64>();
            weights.iter_mut().for_each(|x| *x /= sum);
            use rand::seq::SliceRandom;
            let new_assignment = *choices
                .choose_weighted(rng, |&k| weights[k] + pick)
                .unwrap();
            if d.cluster != new_assignment {
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
    ref_unit: &Unit,
    seed: u64,
) -> f64 {
    let dim = (c.cluster_num.max(ref_unit.cluster_num), chain_len);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut count: u32 = 0;
    let rng = &mut rng;
    let start = std::time::Instant::now();
    let stable_thr = match c.read_type {
        ReadType::CCS => (data.len() as f64 * 0.02).max(1.).floor() as u32,
        ReadType::ONT => (data.len() as f64 * 0.1).max(4.).floor() as u32,
        ReadType::CLR => (data.len() as f64 * 0.1).max(4.).floor() as u32,
    };
    let repeat_num = if c.retain_current_clustering {
        1
    } else {
        REPEAT_NUM
    };
    let (pos, lks, margin) = (0..repeat_num)
        .map(|_| {
            if !c.retain_current_clustering {
                data.iter_mut().for_each(|cs| {
                    cs.cluster = rng.gen_range(0, dim.0);
                });
            }
            let vn = 3 * c.variant_num;
            let (_, pos, lks, margin) = get_variants(&data, dim, rng, c, vn);
            (pos, lks, margin)
        })
        .max_by(|x, y| (x.2).partial_cmp(&y.2).unwrap())
        .unwrap();
    trace!("Init Margin:{}", margin);
    let lks = average_correction(lks, dim);
    update_by_precomputed_lks(data, dim, 0., &pos, &lks, rng);
    loop {
        let vn = 3 * c.variant_num - c.variant_num * count as usize / (c.stable_limit as usize - 1);
        let (_, pos, lks, margin) = get_variants(&data, dim, rng, c, vn);
        let lks = average_correction(lks, dim);
        let changed_num = update_by_precomputed_lks(data, dim, 0., &pos, &lks, rng);
        count += (changed_num <= stable_thr) as u32;
        count *= (changed_num <= stable_thr) as u32;
        report(ref_unit.id, &data, count, margin, 1., changed_num);
        if (std::time::Instant::now() - start).as_secs() > c.limit {
            info!("Break by timelimit:{}", c.limit);
            break margin;
        } else if count >= c.stable_limit {
            break margin;
        }
    }
}

#[allow(dead_code)]
fn likelihoods<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    (cluster_num, chain_len): (usize, usize),
    c: &ClusteringConfig<F>,
    ref_unit: u64,
) -> f64 {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(ref_unit);
    let dim = (cluster_num, chain_len);
    let position = vec![true; chain_len];
    let poa = create_model::get_models(data, dim, &mut rng, c, &position, None);
    let mut fraction = vec![0.1; cluster_num];
    for read in data.iter() {
        fraction[read.cluster] += 1.;
    }
    fraction.iter_mut().for_each(|x| *x /= data.len() as f64);
    data.par_iter()
        .map(|read| {
            let lks: Vec<_> = (0..cluster_num)
                .map(|cl| {
                    fraction[cl].ln()
                        + read
                            .chunks
                            .iter()
                            .map(|chunk| poa[cl][chunk.pos].forward(&chunk.seq, &c.poa_config))
                            .sum::<f64>()
                })
                .collect();
            logsumexp(&lks)
        })
        .sum::<f64>()
}

#[allow(dead_code)]
fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

fn average_correction(
    mut lks: Vec<Vec<f64>>,
    (cluster_num, chain_len): (usize, usize),
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
