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
        let max_unit = self.selected_chunks.iter().map(|u| u.id).max().unwrap() + 1;
        let mut pileups: Vec<Vec<(_, _, &Node)>> = vec![vec![]; max_unit as usize];
        for read in self.encoded_reads.iter() {
            for (idx, node) in read.nodes.iter().enumerate() {
                pileups[node.unit as usize].push((read.id, idx, node));
            }
        }
        let reads_as_paths: Vec<_> = {
            let mut clustered_chunks: HashMap<u64, Vec<(usize, Elm)>> = HashMap::new();
            for (unit_id, units) in pileups.into_iter().enumerate() {
                // debug!("The {}-th pileup", unit_id);
                let ref_unit = match self.selected_chunks.iter().find(|u| u.id == unit_id as u64) {
                    Some(res) => res,
                    None => continue,
                };
                if false {
                    let _assignments = unit_clustering(&units, c, ref_unit);
                }
                let mut rng: Xoshiro256StarStar = rand::SeedableRng::seed_from_u64(100);
                let assignments: Vec<_> = units
                    .iter()
                    .map(|_| rng.gen_range(0, c.cluster_num))
                    .collect();
                if log_enabled!(log::Level::Trace) {
                    debug!("Dump result....");
                    let id_to_name: HashMap<_, _> =
                        self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
                    for cl in 0..c.cluster_num {
                        for (_, (readid, _, _)) in assignments
                            .iter()
                            .zip(units.iter())
                            .filter(|&(&a, _)| a == cl)
                        {
                            debug!("{}\t{}\t{}", cl, readid, id_to_name[&readid]);
                        }
                    }
                }
                for (cluster, (readid, pos, node)) in assignments.into_iter().zip(units) {
                    let elm = Elm::new(node.unit, cluster);
                    clustered_chunks.entry(readid).or_default().push((pos, elm))
                }
            }
            clustered_chunks
                .into_iter()
                .map(|(id, mut units)| {
                    units.sort_by_key(|e| e.0);
                    let path: Vec<_> = units.into_iter().map(|e| e.1).collect();
                    let cluster = 0;
                    ERead { id, path, cluster }
                })
                .collect()
        };
        let assignments = path_clustering(reads_as_paths, c);
        if log_enabled!(log::Level::Debug) {
            let id_to_name: HashMap<_, _> =
                self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
            for cl in 0..c.cluster_num {
                for asn in assignments.iter().filter(|asn| asn.cluster == cl) {
                    debug!("{}\t{}\t{}", cl, asn.id, id_to_name[&asn.id]);
                }
            }
        }
        self.assignments = assignments;
        self
    }
}
// Clustering units.
pub fn unit_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    units: &[(u64, usize, &Node)],
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
) -> Vec<usize> {
    let len = if ref_unit.seq().len() % c.subchunk_length == 0 {
        ref_unit.seq().len() / c.subchunk_length
    } else {
        ref_unit.seq().len() / c.subchunk_length + 1
    };
    let data: Vec<ChunkedUnit> = units
        .into_iter()
        .map(|(_, _, node)| node_to_subchunks(node, c.subchunk_length))
        .map(|chunks| ChunkedUnit { cluster: 0, chunks })
        .collect();
    clustering_by_kmeans(data, len, c, c.seed * ref_unit.id)
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
        .filter_map(|(idx, w)| {
            let chunk = Chunk {
                pos: idx,
                seq: node.seq()[w[0]..w[1]].to_vec(),
            };
            Some(chunk)
        })
        .collect()
}

fn clustering_by_kmeans<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    mut data: Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    seed: u64,
) -> Vec<usize> {
    // Construct Pileups.
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    data.iter_mut().for_each(|cs| {
        cs.cluster = rng.gen_range(0, c.cluster_num);
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
        assert_eq!(chain_len, pos.len());
        for beta in betas.iter() {
            for bs in beta.iter() {
                assert_eq!(bs.len(), chain_len);
            }
        }
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
    data.iter().map(|d| d.cluster).collect()
}

fn report(id: u64, data: &[ChunkedUnit], count: u32, lk: f64, beta: f64) {
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

struct Graph {
    nodes: Vec<GNode>,
    node_num: usize,
}

impl Graph {
    fn new(reads_as_paths: &[ERead], cluster: usize, cluster_num: usize) -> Self {
        let node_num = reads_as_paths
            .iter()
            .flat_map(|r| r.path.iter())
            .map(|elm| elm.unit)
            .max()
            .unwrap() as usize
            + 1;
        let nodes: Vec<_> = (0..cluster_num * node_num).map(|_| GNode::new()).collect();
        let graph = Graph { nodes, node_num };
        reads_as_paths
            .iter()
            .filter(|r| r.cluster == cluster)
            .fold(graph, |g, read| g.register(read))
    }
    fn register(mut self, read: &ERead) -> Self {
        for w in read.path.windows(2) {
            let from = w[0].unit as usize + w[0].cluster * self.node_num;
            let to = w[1].unit as usize + w[1].cluster * self.node_num;
            self.nodes[from].register(to);
        }
        self
    }
    fn score(&self, path: &ERead) -> f64 {
        path.path
            .windows(2)
            .map(|w| {
                let from = w[0].unit as usize + w[0].cluster * self.node_num;
                let to = w[1].unit as usize + w[1].cluster * self.node_num;
                let edg = &self.nodes[from];
                edg.edges
                    .iter()
                    .find(|e| e.to == to)
                    .map(|e| e.weight as f64 / edg.total as f64)
                    .unwrap_or(0.)
            })
            .sum::<f64>()
    }
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let node_num = self.nodes.len();
        let edge_num = self.nodes.iter().map(|n| n.edges.len()).count();
        write!(f, "NodeNum:{}\tEdgeNum:{}", node_num, edge_num)
    }
}

struct GNode {
    edges: Vec<GEdge>,
    total: u64,
}

impl GNode {
    fn new() -> Self {
        Self {
            edges: vec![],
            total: 0,
        }
    }
    fn register(&mut self, to: usize) {
        self.total += 1;
        if let Some(edge) = self.edges.iter_mut().find(|e| e.to == to) {
            edge.weight += 1;
        } else {
            self.edges.push(GEdge { to, weight: 1 });
        }
    }
}

struct GEdge {
    to: usize,
    weight: u64,
}

fn path_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    mut reads_as_paths: Vec<ERead>,
    c: &ClusteringConfig<F>,
) -> Vec<Assignment> {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(c.seed);
    reads_as_paths.iter_mut().for_each(|cs| {
        cs.cluster = rng.gen_range(0, c.cluster_num);
    });
    debug!("Start path clustering");
    let mut count = 0;
    let start = std::time::Instant::now();
    let stable_thr = (reads_as_paths.len() / 100).max(4) as u32;
    while count < c.stable_limit {
        let models: Vec<_> = (0..c.cluster_num)
            .map(|cl| Graph::new(&reads_as_paths, cl, c.cluster_num))
            .collect();
        let change_num = reads_as_paths
            .iter_mut()
            .map(|read| {
                let scores: Vec<_> = models.iter().map(|m| m.score(read)).collect();
                let line: Vec<_> = scores.iter().map(|e| format!("{:.3}", e)).collect();
                if scores.iter().all(|&s| s < 0.001) {
                    debug!("{:?}", read);
                }
                debug!("{}\t{}", read.id, line.join(","));
                let (argmax, _) = scores
                    .iter()
                    .enumerate()
                    .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
                    .unwrap();
                debug!("{}", argmax);
                let changed = 1 - (read.cluster == argmax) as u32;
                read.cluster = argmax;
                changed
            })
            .sum::<u32>();
        count += (change_num < stable_thr) as u32;
        debug!("ChangeNum:{},{}", change_num, stable_thr);
        let elapsed = (std::time::Instant::now() - start).as_secs();
        if elapsed > c.limit && count < c.stable_limit / 2 {
            info!("Break by timelimit:{:?}", elapsed);
            break;
        }
    }
    reads_as_paths
        .iter()
        .map(|read| Assignment::new(read.id, read.cluster))
        .collect()
}
