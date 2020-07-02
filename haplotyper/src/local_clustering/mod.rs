#![allow(dead_code)]
use definitions::*;
use poa_hmm::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::HashMap;
const SMALL_WEIGHT: f64 = 0.000_000_001;
mod config;
pub use config::*;
mod create_model;
mod eread;
mod variant_calling;
use create_model::*;
use eread::*;
use variant_calling::*;
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
        rayon::ThreadPoolBuilder::new()
            .num_threads(c.threads)
            .build_global()
            .unwrap();
        let max_unit = self.selected_chunks.iter().map(|u| u.id).max().unwrap() + 1;
        let mut pileups: Vec<Vec<(_, _, &mut Node)>> = (0..max_unit).map(|_| Vec::new()).collect();
        for read in self.encoded_reads.iter_mut() {
            for (idx, node) in read.nodes.iter_mut().enumerate() {
                pileups[node.unit as usize].push((read.id, idx, node));
            }
        }
        let selected_chunks = self.selected_chunks.clone();
        let id_to_name: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
        pileups
            .par_iter_mut()
            .enumerate()
            .for_each(|(unit_id, mut units)| {
                debug!("The {}-th pileup", unit_id);
                let ref_unit = match selected_chunks.iter().find(|u| u.id == unit_id as u64) {
                    Some(res) => res,
                    None => return,
                };
                unit_clustering(&mut units, c, ref_unit);
                if log_enabled!(log::Level::Debug) {
                    for cl in 0..c.cluster_num {
                        let cl = cl as u64;
                        for (id, _, _) in units.iter().filter(|&(_, _, n)| n.cluster == cl) {
                            debug!("{}\t{}\t{}", cl, id, id_to_name[&id]);
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
    let mut data: Vec<ChunkedUnit> = units
        .iter()
        .map(|(_, _, node)| node_to_subchunks(node, c.subchunk_length))
        .map(|chunks| ChunkedUnit { cluster: 0, chunks })
        .collect();
    clustering_by_kmeans(&mut data, len, c, ref_unit.id);
    assert_eq!(data.len(), units.len());
    for ((_, _, ref mut n), d) in units.iter_mut().zip(data.iter()) {
        n.cluster = d.cluster as u64;
    }
}

fn node_to_subchunks(node: &Node, len: usize) -> Vec<Chunk> {
    let (mut q_pos, mut r_pos) = (0, 0);
    let mut target = len;
    let mut ops = node.cigar.clone();
    ops.reverse();
    let mut chunk_position = vec![0];
    loop {
        // Pop until target.
        let last_op = {
            loop {
                match ops.pop() {
                    Some(Op::Del(l)) => {
                        r_pos += l;
                        if target <= r_pos {
                            break Op::Del(l);
                        }
                    }
                    Some(Op::Ins(l)) => {
                        q_pos += l;
                    }
                    Some(Op::Match(l)) => {
                        r_pos += l;
                        q_pos += l;
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
                }
                Op::Match(_) => {
                    ops.push(Op::Match(overflow));
                    r_pos -= overflow;
                    q_pos -= overflow;
                }
                _ => unreachable!(),
            }
        }
        chunk_position.push(q_pos);
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
            chunk_position.push(node.seq().len());
            break;
        } else if rest < len && query_rest == 0 {
            break;
        } else {
            target += len;
        }
    }
    chunk_position
        .windows(2)
        .enumerate()
        .map(|(idx, w)| Chunk {
            pos: idx,
            seq: node.seq()[w[0]..w[1]].to_vec(),
        })
        .collect()
}

fn clustering_by_kmeans<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    ref_unit: u64,
) {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(ref_unit * 101);
    data.iter_mut().for_each(|cs| {
        cs.cluster = rng.gen_range(0, c.cluster_num);
    });
    let id = ref_unit;
    let mut beta = c.initial_beta;
    let mut count = 0;
    let mut lk;
    let rng = &mut rng;
    let start = std::time::Instant::now();
    let stable_thr = (data.len() as f64 * c.sample_rate / 2.).max(2.) as u32;
    while count < c.stable_limit {
        let (betas, pos, next_lk) = get_variants(&data, chain_len, rng, c, 1.);
        beta = (beta * c.beta_increase).min(c.max_beta);
        lk = next_lk;
        let changed_num = (0..c.sample_rate.recip().ceil() as usize / 2)
            .map(|_| {
                let update_data: Vec<_> = (0..data.len())
                    .map(|_| rng.gen_bool(c.sample_rate))
                    .collect::<Vec<_>>();
                let ms = get_models(&data, chain_len, rng, c, &pos, &update_data);
                update_assignment(data, chain_len, c, &update_data, &betas, beta, &ms)
            })
            .sum::<u32>();
        report(id, &data, count, lk, beta, changed_num);
        count += (changed_num <= stable_thr) as u32;
        count *= (changed_num <= stable_thr) as u32;
        let elapsed = (std::time::Instant::now() - start).as_secs();
        if elapsed > c.limit && count < c.stable_limit / 2 {
            info!("Break by timelimit:{:?}", elapsed);
            break;
        }
    }
}

fn report(id: u64, data: &[ChunkedUnit], count: u32, lk: f64, beta: f64, c: u32) {
    if !log_enabled!(log::Level::Trace) {
        return;
    }
    trace!("{}\t{}\t{}\t{}\t{}", id, count, lk, beta, c);
    let mut count: std::collections::HashMap<usize, usize> = HashMap::new();
    for read in data {
        *count.entry(read.cluster).or_default() += 1;
    }
    let mut count: Vec<(usize, usize)> = count.into_iter().collect();
    count.sort_by_key(|e| e.0);
    let count: Vec<_> = count.iter().map(|(_, cnt)| format!("{}", cnt)).collect();
    trace!("{}", count.join("\t"));
}

fn get_fraction_on_position(
    data: &[ChunkedUnit],
    chain_len: usize,
    cluster_num: usize,
) -> Vec<Vec<f64>> {
    let mut total_count = vec![0; chain_len];
    let mut counts = vec![vec![0; chain_len]; cluster_num];
    for read in data.iter() {
        for chunk in read.chunks.iter() {
            total_count[chunk.pos] += 1;
            counts[read.cluster][chunk.pos] += 1;
        }
    }
    counts
        .iter()
        .map(|cs| {
            cs.iter()
                .zip(total_count.iter())
                .map(|(&c, &t)| c as f64 / t as f64 + SMALL_WEIGHT)
                .collect()
        })
        .collect()
}

fn get_argmax(ws: &[f64]) -> Option<usize> {
    ws.iter()
        .enumerate()
        .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Less))
        .map(|e| e.0)
}

fn update_assignment<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut [ChunkedUnit],
    chain_len: usize,
    c: &ClusteringConfig<F>,
    update_data: &[bool],
    betas: &[Vec<Vec<f64>>],
    beta: f64,
    models: &[Vec<POA>],
) -> u32 {
    let fraction_on_positions = get_fraction_on_position(data, chain_len, c.cluster_num);
    data.iter_mut()
        .zip(update_data.iter())
        .filter(|&(_, &b)| b)
        .map(|(read, _)| {
            let likelihoods: Vec<(usize, Vec<_>)> = read
                .chunks
                .iter()
                .map(|chunk| {
                    let lks = models
                        .iter()
                        .map(|ms| ms[chunk.pos].forward(&chunk.seq, &c.poa_config))
                        .collect();
                    (chunk.pos, lks)
                })
                .collect();
            let ws: Vec<_> = fraction_on_positions
                .iter()
                .map(|ws| {
                    read.chunks.iter().map(|c| ws[c.pos]).sum::<f64>() / read.chunks.len() as f64
                })
                .collect();
            let weights: Vec<_> = (0..c.cluster_num)
                .map(|l| {
                    (0..c.cluster_num)
                        .map(|k| {
                            if k == l {
                                return 0.;
                            }
                            let (i, j) = (l.max(k), l.min(k));
                            let prior = beta * (ws[k].ln() - ws[l].ln());
                            prior
                                + likelihoods
                                    .iter()
                                    .map(|&(pos, ref lks)| betas[i][j][pos] * (lks[k] - lks[l]))
                                    .sum::<f64>()
                        })
                        .map(|lkdiff| lkdiff.exp())
                        .sum::<f64>()
                        .recip()
                })
                .collect();
            let argmax = get_argmax(&weights).unwrap_or(0);
            let is_the_same = if read.cluster == argmax { 0 } else { 1 };
            read.cluster = argmax;
            is_the_same
        })
        .sum::<u32>()
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
