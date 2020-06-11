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
use create_model::*;
mod variant_calling;
use variant_calling::*;
mod eread;
use eread::*;

pub trait Clustering {
    fn clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(self, c: &ClusteringConfig<F>) -> Self;
    // fn serialize_subchunk<F: Fn(u8, u8) -> i32>(
    //     &self,
    //     c: &ClusteringConfig<F>,
    // ) -> (usize, Vec<ERead>);
}

impl Clustering for DataSet {
    fn clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        mut self,
        c: &ClusteringConfig<F>,
    ) -> Self {
        rayon::ThreadPoolBuilder::new()
            .num_threads(c.threads)
            .build_global()
            .unwrap();
        let max_unit = self.selected_chunks.iter().map(|u| u.id).max().expect("u");
        let mut pileups: Vec<(_, _, &Node)> = vec![vec![]; max_unit as usize + 1];
        for read in self.encoded_reads.iter() {
            for (idx, node) in read.nodes.iter() {
                pileup[node.unit as usize].push((read.id, idx, node));
            }
        }
        let reads_as_paths:Vec<_> = {
            use std::collections::HashMap;
            let mut clustered_chunks: HashMap<u64, Vec<(usize, Elm)>> = HashMap::new();
            for pileup in 
        };
        let assignments = clustering_by_kmeans(data, chain_len, c);
        self.assignments = assignments;
        self
    }
    // fn serialize_subchunk<F: Fn(u8, u8) -> i32>(
    //     &self,
    //     c: &ClusteringConfig<F>,
    // ) -> (usize, Vec<ERead>) {
    //     // (UnitID, subchunkID) -> Serialized Position
    //     let locationmap: HashMap<(u64, usize), usize> = self
    //         .selected_chunks
    //         .iter()
    //         .flat_map(|u| {
    //             let len = if u.seq().len() % c.subchunk_length == 0 {
    //                 u.seq().len() / c.subchunk_length
    //             } else {
    //                 u.seq().len() / c.subchunk_length + 1
    //             };
    //             (0..len).map(|i| (u.id, i)).collect::<Vec<_>>()
    //         })
    //         .enumerate()
    //         .map(|(x, y)| (y, x))
    //         .collect();
    //     let chain_length = *locationmap.values().max().expect("chain length") + 1;
    //     let data = self
    //         .encoded_reads
    //         .iter()
    //         .map(|read| split_into_subchunks(read, c.subchunk_length, &locationmap))
    //         .collect();
    //     (chain_length, data)
    // }
}

fn split_into_subchunks(
    read: &EncodedRead,
    len: usize,
    map: &HashMap<(u64, usize), usize>,
) -> ERead {
    let chunks: Vec<_> = read
        .nodes
        .iter()
        .flat_map(|node| node_to_subchunks(node, len, map))
        .collect();
    ERead {
        id: read.id,
        chunks,
        cluster: 0,
    }
}

fn node_to_subchunks(node: &Node, len: usize, map: &HashMap<(u64, usize), usize>) -> Vec<Chunk> {
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
        .filter_map(|(idx, w)| {
            let pos = *map.get(&(node.unit, idx))?;
            let chunk = Chunk {
                pos,
                seq: node.seq()[w[0]..w[1]].to_vec(),
            };
            Some(chunk)
        })
        .collect()
}

fn clustering_by_kmeans<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    mut data: Vec<ERead>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
) -> Vec<Assignment> {
    // Construct Pileups.
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(c.seed);
    data.iter_mut().for_each(|r| {
        r.cluster = rng.gen_range(0, c.cluster_num);
    });
    let mut beta = c.initial_beta;
    let mut count = 0;
    let mut lk = std::f64::NEG_INFINITY;
    let rng = &mut rng;
    let start = std::time::Instant::now();
    let stable_thr = (data.len() as f64 * c.sample_rate / 2.).max(4.) as u32;
    while count < c.stable_limit {
        debug!("Calc Variants");
        let (betas, pos, next_lk) = get_variants(&data, chain_len, rng, c, 2.);
        beta = (beta * c.beta_increase).min(c.max_beta);
        report(c.id, &data, count, lk, beta);
        lk = next_lk;
        let changed_num = (0..c.sample_rate.recip().ceil() as usize / 2)
            .map(|_| {
                let update_data: Vec<_> = (0..data.len())
                    .map(|_| rng.gen_bool(c.sample_rate))
                    .collect::<Vec<_>>();
                let ms = get_models(&data, chain_len, rng, c, &pos, &update_data);
                update_assignment(&mut data, chain_len, c, &update_data, &betas, beta, &ms)
            })
            .sum::<u32>();
        debug!("CHANGENUM\t{}", changed_num);
        count += (changed_num <= stable_thr) as u32;
        let elapsed = (std::time::Instant::now() - start).as_secs();
        if elapsed > c.limit && count < c.stable_limit / 2 {
            info!("Break by timelimit:{:?}", elapsed);
            break;
        }
    }
    data.iter()
        .map(|d| Assignment::new(d.id, d.cluster as u32))
        .collect()
}

fn report(id: u64, data: &[ERead], count: u32, lk: f64, beta: f64) {
    debug!("{}\t{}\t{}\t{}", id, count, lk, beta);
    let mut count: std::collections::HashMap<usize, usize> = HashMap::new();
    for read in data {
        *count.entry(read.cluster).or_default() += 1;
    }
    let mut count: Vec<(usize, usize)> = count.into_iter().collect();
    count.sort_by_key(|e| e.0);
    for (cl, cnt) in count {
        debug!("{}:{}", cl, cnt);
    }
}

fn get_fraction_on_position(data: &[ERead], chain_len: usize, cluster_num: usize) -> Vec<Vec<f64>> {
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
    data: &mut [ERead],
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
                .par_iter()
                .map(|chunk| {
                    let lks = models
                        .par_iter()
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
