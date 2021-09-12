use definitions::*;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
mod config;
mod var_clustering;
pub use config::*;
pub mod clustering_by_assemble;
pub mod create_model;
pub mod eread;
pub mod pca;
pub mod variant_calling;
use eread::*;
use variant_calling::*;
pub mod kmeans;

/// Return rand index.
pub fn rand_index(label: &[u8], pred: &[u8]) -> f64 {
    assert_eq!(label.len(), pred.len());
    let mut both_same_pair = 0;
    let mut both_diff_pair = 0;
    for (i, (l1, p1)) in label.iter().zip(pred.iter()).enumerate() {
        for (l2, p2) in label.iter().zip(pred.iter()).take(i) {
            if l1 == l2 && p1 == p2 {
                both_same_pair += 1;
            } else if l1 != l2 && p1 != p2 {
                both_diff_pair += 1;
            }
        }
    }
    let len = label.len();
    (both_same_pair + both_diff_pair) as f64 / (len * (len - 1) / 2) as f64
}

pub trait LocalClustering {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        self,
        c: &ClusteringConfig<F>,
    ) -> Self;
}

impl LocalClustering for DataSet {
    fn local_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        mut self,
        _c: &ClusteringConfig<F>,
    ) -> Self {
        local_clustering_all(&mut self);
        self
    }
}

pub fn local_clustering_selected(ds: &mut DataSet, selection: &HashSet<u64>) {
    let mut pileups: HashMap<u64, Vec<&mut Node>> =
        selection.iter().map(|&id| (id, vec![])).collect();
    let cluster_num: HashMap<u64, u8> = ds
        .selected_chunks
        .iter()
        .filter(|c| selection.contains(&c.id))
        .map(|c| (c.id, c.cluster_num as u8))
        .collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some(bucket) = pileups.get_mut(&node.unit) {
            bucket.push(node);
        }
    }
    let coverage = ds.coverage;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .map(|(&unit_id, units)| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 25);
            let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
            let coverage = match coverage {
                Some(res) => res,
                None => {
                    debug!("No coverage estimation. Use adhoc coverage");
                    (units.len() / 2) as f64
                }
            };
            let config = kmeans::ClusteringConfig::new(100, cluster_num[&unit_id], coverage);
            let start = std::time::Instant::now();
            let cov = seqs.len();
            let (asn, consensus, score) = kmeans::clustering(&seqs, &mut rng, &config).unwrap();
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            let end = std::time::Instant::now();
            let elapsed = (end - start).as_secs();
            let len = consensus.len();
            debug!(
                "RECORD\t{}\t{}\t{}\t{:.3}\t{}",
                unit_id, elapsed, len, score, cov
            );

            // let (prevcl, cl) = (cluster_num[&unit_id], config.cluster_num);
            // debug!(
            //     "RECORD\t{}\t{}\t{}\t{}\t{:.3}\t{}",
            //     unit_id, elapsed, prevcl, cl, score, cov
            // );
            (unit_id, (consensus, score, config.cluster_num))
        })
        .collect();
    for unit in ds.selected_chunks.iter_mut() {
        if let Some((consensus, score, cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = String::from_utf8(consensus.to_vec()).unwrap();
            unit.cluster_num = *cluster_num as usize;
            unit.score = *score;
        }
    }
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some((cons, _, _)) = consensus_and_clusternum.get(&node.unit) {
            let band_size = (cons.len() / 10).max(5);
            let (_, cigar) =
                kiley::bialignment::global_banded(cons, node.seq(), 2, -2, -4, -2, band_size);
            node.cigar = crate::encode::compress_kiley_ops(&cigar);
        }
    }
}

pub fn local_clustering_all(ds: &mut DataSet) {
    let mut pileups: HashMap<u64, Vec<&mut Node>> = HashMap::new();
    let cluster_num: HashMap<u64, u8> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num as u8))
        .collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        pileups.entry(node.unit).or_default().push(node);
    }
    let coverage = ds.coverage;
    let consensus_and_clusternum: HashMap<_, _> = pileups
        .par_iter_mut()
        .map(|(&unit_id, units)| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(unit_id * 23);
            let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
            let coverage = match coverage {
                Some(res) => res,
                None => {
                    debug!("No coverage estimation. Use adhoc coverage");
                    (units.len() / 2) as f64
                }
            };
            let config = kmeans::ClusteringConfig::new(100, cluster_num[&unit_id], coverage);
            let start = std::time::Instant::now();
            let (asn, consensus, score) = kmeans::clustering(&seqs, &mut rng, &config).unwrap();
            for (node, asn) in units.iter_mut().zip(asn) {
                node.cluster = asn as u64;
            }
            let end = std::time::Instant::now();
            let elapsed = (end - start).as_secs();
            let (prevcl, cl) = (cluster_num[&unit_id], config.cluster_num);
            debug!(
                "RECORD\t{}\t{}\t{}\t{}\t{:.3}",
                unit_id, elapsed, prevcl, cl, score
            );
            (unit_id, (consensus, score, config.cluster_num))
        })
        .collect();
    for unit in ds.selected_chunks.iter_mut() {
        if let Some((consensus, score, cluster_num)) = consensus_and_clusternum.get(&unit.id) {
            unit.seq = String::from_utf8(consensus.to_vec()).unwrap();
            unit.cluster_num = *cluster_num as usize;
            unit.score = *score;
        }
    }
    // Modify alignment so that the alignment would be valid.
    // We would do this very rough, because anyway we do not like to use this alignment in the future.
    // TODO: Make this banded procedure into SWG with banded mode.
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some((cons, _, _)) = consensus_and_clusternum.get(&node.unit) {
            let band_size = (cons.len() / 10).max(5);
            let (_, cigar) =
                kiley::bialignment::global_banded(cons, node.seq(), 2, -2, -4, -2, band_size);
            node.cigar = crate::encode::compress_kiley_ops(&cigar);
        }
    }
}

pub fn unit_clustering<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    units: &mut [(u64, usize, &mut Node)],
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
    iteration: u64,
) -> usize {
    if ref_unit.cluster_num == 1 {
        units.iter_mut().for_each(|x| x.2.cluster = 0);
        return ref_unit.cluster_num;
    }
    let len = if ref_unit.seq().len() % c.subchunk_length == 0 {
        ref_unit.seq().len() / c.subchunk_length
    } else {
        ref_unit.seq().len() / c.subchunk_length + 1
    };
    let (min_length, max_length) = ((c.subchunk_length * 9) / 10, c.subchunk_length * 11 / 10);
    let mut data: Vec<_> = units
        .iter()
        .map(|&(_, _, ref node)| {
            let mut chunks = node_to_subchunks(node, c.subchunk_length);
            chunks.retain(|chunk| min_length < chunk.seq.len() && chunk.seq.len() < max_length);
            assert!(chunks.iter().all(|x| !x.seq.is_empty()), "{:?}", chunks);
            let cluster = if c.retain_current_clustering {
                node.cluster as usize
            } else {
                0
            };
            ChunkedUnit { cluster, chunks }
        })
        .collect();
    let seed = ref_unit.id * iteration;
    clustering_by_kmeans(&mut data, len, c, ref_unit, seed);
    for ((_, _, ref mut n), d) in units.iter_mut().zip(data.iter()) {
        n.cluster = d.cluster as u64;
    }
    ref_unit.cluster_num
}

fn to_pileup(node: &Node, unit: &Unit) -> Vec<u8> {
    let mut slots = vec![0; unit.seq().len()];
    let qseq = node.seq();
    let (mut qpos, mut rpos) = (0, 0);
    for op in node.cigar.iter() {
        match *op {
            Op::Del(l) => rpos += l,
            Op::Ins(l) => qpos += l,
            Op::Match(l) => {
                slots[rpos..rpos + l].clone_from_slice(&qseq[qpos..qpos + l]);
                // for p in 0..l {
                //     slots[rpos + p] = qseq[qpos + p];
                // }
                rpos += l;
                qpos += l;
            }
        }
    }
    slots
}

fn clustering_variant_vector(vars: &[Vec<i8>], k: usize) -> Vec<u64> {
    use var_clustering::clustering_variant_to;
    let (mut asn, mut aic) = (vec![], std::f64::NEG_INFINITY);
    for cl in k..2 * k {
        for s in 0..10 {
            let (new_asn, new_aic) = clustering_variant_to(vars, cl, s + cl as u64);
            if new_aic > aic {
                asn = new_asn;
                aic = new_aic;
            }
        }
    }
    asn
}

pub fn unit_clustering_ccs_kmervec<T: std::borrow::Borrow<[u8]>>(
    xs: &[T],
    cluster_num: u8,
    k: u8,
    initial_assignments: &[u8],
) -> Vec<u8> {
    let mut assignments = initial_assignments.to_vec();
    let mut feature_vectors: Vec<_> = xs
        .iter()
        .map(|x| {
            let mut slots = vec![0f64; 4usize.pow(k as u32)];
            for w in x.borrow().windows(k as usize) {
                let location = w
                    .iter()
                    .map(|&x| match x {
                        b'A' | b'a' => 0,
                        b'C' | b'c' => 1,
                        b'G' | b'g' => 2,
                        b'T' | b't' => 3,
                        _ => 3,
                    })
                    .fold(0, |acc, x| (acc << 2) | x);
                slots[location] += 1f64;
            }
            slots
        })
        .collect();
    // Centering;
    let mut center = vec![0f64; 4usize.pow(k as u32)];
    for fv in feature_vectors.iter() {
        for (x, &y) in center.iter_mut().zip(fv.iter()) {
            *x += y;
        }
    }
    // fn usize_to_kmer(x: usize, k: u8) -> Vec<u8> {
    //     (0..k)
    //         .fold((vec![], x), |(mut kmer, rest), _| {
    //             match rest & 0b11 {
    //                 0 => kmer.push(b'A'),
    //                 1 => kmer.push(b'C'),
    //                 2 => kmer.push(b'G'),
    //                 _ => kmer.push(b'T'),
    //             }
    //             (kmer, rest >> 2)
    //         })
    //         .0
    // }
    // center
    //     .iter_mut()
    //     .for_each(|x| *x /= feature_vectors.len() as f64);
    // for (idx, count) in center.iter().enumerate() {
    //     if 2f64 < count * feature_vectors.len() as f64 {
    //         let kmer = String::from_utf8(usize_to_kmer(idx, k)).unwrap();
    //         println!("CENTER\t{}\t{}", kmer, count);
    //     }
    // }
    // feature_vectors
    //     .iter_mut()
    //     .for_each(|fv| fv.iter_mut().zip(center.iter()).for_each(|(x, y)| *x -= y));
    // for (idx, fv) in feature_vectors.iter().enumerate() {
    //     for (kmer, val) in fv.iter().enumerate() {
    //         if 0.25 < val.abs() {
    //             let kmer = String::from_utf8(usize_to_kmer(kmer, k)).unwrap();
    //             println!("FV\t{}\t{}\t{}", idx, kmer, val);
    //         }
    //     }
    // }
    let mask_kmer: Vec<_> = center.iter().map(|&count| 2f64 < count).collect();
    feature_vectors.iter_mut().for_each(|fv| {
        fv.iter_mut().zip(mask_kmer.iter()).for_each(|(x, &y)| {
            if y {
                *x = 0f64;
            }
        })
    });
    for i in 0..100 {
        let mut centers = vec![vec![0f64; 4usize.pow(k as u32)]; cluster_num as usize];
        let mut counts = vec![0; cluster_num as usize];
        for (fv, &cl) in feature_vectors.iter().zip(assignments.iter()) {
            for (x, y) in centers[cl as usize].iter_mut().zip(fv.iter()) {
                *x += y;
            }
            counts[cl as usize] += 1;
        }
        for (center, &count) in centers.iter_mut().zip(counts.iter()) {
            center.iter_mut().for_each(|x| *x /= count as f64);
        }
        let new_assignments: Vec<_> = feature_vectors
            .iter()
            .map(|fv| {
                centers
                    .iter()
                    .enumerate()
                    .map(|(cl, center)| {
                        let dist = center
                            .iter()
                            .zip(fv.iter())
                            .map(|(x, y)| (x - y).powi(2))
                            .sum::<f64>();
                        (cl as u8, dist)
                    })
                    .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
                    .unwrap()
                    .0
            })
            .collect();
        if new_assignments == assignments {
            break;
        } else {
            assignments = new_assignments;
        }
        let nums: Vec<_> = (0..cluster_num)
            .map(|cl| {
                let count = bytecount::count(&assignments, cl);
                // let count = assignments.iter().filter(|&&x| x == cl).count();
                format!("{}", count)
            })
            .collect();
        println!("ITER\t{}\t{}", i, nums.join("\t"));
    }
    assignments
}

pub fn unit_clustering_ccs<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    units: &mut [(u64, usize, &mut Node)],
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
    _iteration: u64,
) -> usize {
    let pileup: Vec<Vec<u8>> = units
        .iter()
        .map(|(_, _, unit)| to_pileup(unit, ref_unit))
        .collect();
    let variants: Vec<(usize, u8, u8)> = {
        let mut maf_fractions: Vec<_> = (0..ref_unit.seq().len())
            .filter_map(|pos| {
                let bases: Vec<u8> = pileup
                    .iter()
                    .filter(|pi| pi[pos] != 0)
                    .map(|pi| pi[pos])
                    .collect();
                let coverage = bases.len();
                if coverage * 3 < pileup.len() {
                    return None;
                }
                let mut count = vec![0; 4];
                for base in bases.iter() {
                    match *base {
                        b'A' => count[0] += 1,
                        b'C' => count[1] += 1,
                        b'G' => count[2] += 1,
                        b'T' => count[3] += 1,
                        x => panic!("{}\t{}", x, line!()),
                    }
                }
                let (major, _) = count.iter().enumerate().max_by_key(|x| x.1).unwrap();
                let (minor, minor_count) = count
                    .iter()
                    .enumerate()
                    .filter(|&(i, _)| i != major)
                    .max_by_key(|x| x.1)
                    .unwrap();
                let frac = *minor_count as f64 / coverage as f64;
                let major = match major {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => panic!(),
                };
                let minor = match minor {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => panic!(),
                };
                Some((pos, frac, major, minor))
            })
            .collect();
        maf_fractions.sort_by(|x, y| {
            let x_var = (x.1) * (1f64 - x.1);
            let y_var = (y.1) * (1f64 - y.1);
            (x_var).partial_cmp(&y_var).unwrap()
        });
        maf_fractions
            .into_iter()
            .rev()
            .take(10)
            .map(|(pos, _, major, minor)| (pos, major, minor))
            .collect()
    };
    for &(pos, major, minor) in variants.iter() {
        let (major, minor) = (major as char, minor as char);
        debug!("VAR\t{}\t{}\t{}\t{}", ref_unit.id, pos, major, minor);
    }
    let variant_vector: Vec<Vec<_>> = pileup
        .iter()
        .map(|pl| {
            variants
                .iter()
                .map(|&(pos, major, minor)| match pl[pos] {
                    x if x == major => 1,
                    x if x == minor => 0,
                    _ => -1,
                })
                .collect()
        })
        .collect();
    for (i, var) in variant_vector.iter().enumerate() {
        let var: String = var
            .iter()
            .map(|x| match *x {
                1 => '1',
                0 => '0',
                _ => '-',
            })
            .collect();
        debug!("{}\t{}\t{}", ref_unit.id, i, var);
    }
    let clusters = clustering_variant_vector(&variant_vector, c.cluster_num);
    for (cl, (_, _, unit)) in clusters.iter().zip(units.iter_mut()) {
        unit.cluster = *cl;
    }
    c.cluster_num
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

pub fn clustering_by_gibbs<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    _ref_unit: &Unit,
    seed: u64,
) -> f64 {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let rng = &mut rng;
    let len = data.len();
    data.iter_mut().for_each(|d| d.cluster = 0);
    let use_position = vec![true; chain_len];
    let prior = create_model::get_models(data, (1, chain_len), rng, c, &use_position, None)
        .pop()
        .unwrap();
    let alpha = 0.5;
    let (mut num_data, mut num_cluster) = (vec![data.len()], 1);
    for t in 0..100 {
        debug!("The {}-th iteration.", t);
        assert_eq!(data.len(), num_data.iter().sum::<usize>());
        assert_eq!(num_data.len(), num_cluster);
        for i in 0..len {
            // Remove the i-th reads from cluster.
            num_data[data[i].cluster] -= 1;
            if num_data[data[i].cluster] == 0 {
                // Remove this cluste.
                let removed = data[i].cluster;
                let last = num_data.pop().unwrap();
                num_data[removed] = last;
                data.iter_mut()
                    .filter(|d| d.cluster == last)
                    .for_each(|d| d.cluster = removed);
            }
            let models = create_model::get_models(
                data,
                (num_cluster, chain_len),
                rng,
                c,
                &use_position,
                Some(i),
            );
            let mut posterior: Vec<_> = num_data
                .iter()
                .zip(models.iter())
                .map(|(&n, model)| {
                    let lk = data[i]
                        .chunks
                        .iter()
                        .map(|ch| model[ch.pos].forward(&ch.seq, &c.poa_config))
                        .sum::<f64>();
                    let dirichlet = (n as f64).ln() - (data.len() as f64 - 1. + alpha).ln();
                    lk + dirichlet
                })
                .collect();
            // We can pre-compute this value, right?
            let new_lk = data[i]
                .chunks
                .iter()
                .map(|ch| prior[ch.pos].forward(&ch.seq, &c.poa_config))
                .sum::<f64>();
            let new_cluster = alpha.ln() - (data.len() as f64 - 1. + alpha).ln();
            posterior.push(new_lk + new_cluster);
            use rand::seq::SliceRandom;
            // debug!("{:?}", posterior);
            let sum = logsumexp(&posterior);
            posterior.iter_mut().for_each(|x| *x = (*x - sum).exp());
            // debug!("{:?}", posterior);
            let picked = {
                let choises: Vec<_> = (0..num_cluster + 1).collect();
                *choises.choose_weighted(rng, |&k| posterior[k]).unwrap()
            };
            if num_cluster < picked {
                num_cluster += 1;
                num_data.push(0);
            }
            data[i].cluster = picked;
            num_data[data[i].cluster] += 1;
        }
    }
    0.
}

pub fn clustering_by_kmeans_em<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    ref_unit: &Unit,
    seed: u64,
) -> usize {
    let mut labels = vec![vec![]; data.len()];
    let tn = 3;
    for i in 0..tn {
        clustering_by_kmeans(data, chain_len, c, ref_unit, seed + i);
        for (i, d) in data.iter().enumerate() {
            labels[i].push(d.cluster);
        }
    }
    // for (i, d) in labels.iter().enumerate() {
    //     trace!("{}\t{:?}", i, d);
    // }
    assert!(labels.iter().all(|xs| xs.len() == tn as usize));
    let (asn, _) = (0..4)
        .map(|s| {
            let param = (s, tn as usize, ref_unit.cluster_num);
            em_clustering(&labels, ref_unit.cluster_num, param)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();

    // let (asn, _) = ((ref_unit.cluster_num / 2).max(1)..=ref_unit.cluster_num * 2)
    //     .flat_map(|cl| {
    //         (0..4)
    //             .map(|s| {
    //                 let param = (cl as u64 + s, tn as usize, ref_unit.cluster_num);
    //                 em_clustering(&labels, cl, param)
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
    //     .unwrap();
    let cluster_num = *asn.iter().max().unwrap();
    for (a, x) in asn.into_iter().zip(data.iter_mut()) {
        x.cluster = a;
    }
    cluster_num
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
        ReadType::None => (data.len() as f64 * 0.1).max(4.).floor() as u32,
    };
    let (pos, lks, margin) = // (0..repeat_num).map(|r|
    {
        if !c.retain_current_clustering {
            initial_clustering(data, ref_unit, dim, seed);
        }
        let vn = 3 * c.variant_num;
        let (_, pos, lks, margin) = get_variants(data, dim, rng, c, vn);
        (pos, lks, margin)
    };
    trace!("Init Margin:{}", margin);
    let lks = average_correction(lks, dim);
    update_by_precomputed_lks(data, dim, 0., &pos, &lks, rng);
    loop {
        let vn = 3 * c.variant_num - c.variant_num * count as usize / (c.stable_limit as usize - 1);
        let (_, pos, lks, margin) = get_variants(data, dim, rng, c, vn);
        let lks = average_correction(lks, dim);
        let changed_num = update_by_precomputed_lks(data, dim, 0., &pos, &lks, rng);
        count += (changed_num <= stable_thr) as u32;
        count *= (changed_num <= stable_thr) as u32;
        report(ref_unit.id, data, count, margin, 1., changed_num);
        if (std::time::Instant::now() - start).as_secs() > c.limit {
            info!("Break by timelimit:{}", c.limit);
            break margin;
        } else if count >= c.stable_limit {
            break margin;
        }
    }
}

pub fn initial_clustering(
    data: &mut [ChunkedUnit],
    _ref_unit: &Unit,
    (cluster, chain_len): (usize, usize),
    seed: u64,
) {
    let maf = convert_to_maf(data, cluster, chain_len, seed);
    let length = maf[0].len();
    let (assignments, _) = (5..10)
        .map(|i| em_clustering(&maf, cluster, (seed + i, length, 2)))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    for (d, a) in data.iter_mut().zip(assignments) {
        d.cluster = a;
    }
}

fn convert_to_maf(
    data: &[ChunkedUnit],
    _cluster: usize,
    chain_len: usize,
    seed: u64,
) -> Vec<Vec<usize>> {
    let mut chunks = vec![vec![]; chain_len];
    for d in data.iter() {
        for c in d.chunks.iter() {
            if !c.seq.is_empty() {
                chunks[c.pos].push(c.seq.as_slice());
            }
        }
    }
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let consensus: Vec<_> = chunks
        .iter_mut()
        .map(|cs| {
            if cs.is_empty() {
                vec![]
            } else {
                let cs: Vec<_> = (0..10)
                    .map(|_| {
                        use rand::seq::SliceRandom;
                        cs.shuffle(&mut rng);
                        let cs = &cs[..cs.len().min(15)];
                        poa_hmm::POA::from_slice_default(cs).consensus()
                    })
                    .collect();
                poa_hmm::POA::from_vec_default(&cs).consensus()
            }
        })
        .collect();
    use bio::alignment::pairwise::Aligner;
    if data
        .iter()
        .any(|d| d.chunks.iter().any(|c| c.pos >= consensus.len()))
    {
        debug!("Error");
        debug!("{}\t{}", consensus.len(), chain_len);
        for (i, d) in data.iter().enumerate() {
            let c: Vec<_> = d.chunks.iter().map(|c| (c.pos, c.seq.len())).collect();
            debug!("{}\t{:?}", i, c);
        }
        for (i, c) in chunks.iter().enumerate() {
            debug!("CHUNK\t{}\t{}", i, c.len());
        }
        panic!();
    }
    let mut aligner = Aligner::new(-3, -1, |x, y| if x == y { 2 } else { -5 });
    let maf: Vec<HashMap<usize, Vec<u8>>> = data
        .iter()
        .map(|d| {
            d.chunks
                .iter()
                .map(|c| (c.pos, alignment(&mut aligner, &c.seq, &consensus[c.pos])))
                .collect()
        })
        .collect();
    let mut summary: Vec<Vec<[u32; 5]>> = {
        let mut lens = vec![vec![]; chain_len];
        for aln in maf.iter() {
            for (&pos, c) in aln.iter() {
                lens[pos].push(c.len());
            }
        }
        lens.iter()
            .map(|lens| {
                let max = *lens.iter().max().unwrap_or(&0);
                vec![[0; 5]; max]
            })
            .collect()
    };
    for aln in maf.iter() {
        for (&pos, c) in aln.iter() {
            for (i, &x) in c.iter().enumerate() {
                let x = match x {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => 4,
                };
                summary[pos][i][x] += 1;
            }
        }
    }
    let mut fractions: Vec<_> = summary
        .iter()
        .flat_map(|smry| {
            smry.iter().map(|clm| {
                let tot = clm.iter().sum::<u32>();
                let (argmax, _) = clm.iter().enumerate().max_by_key(|x| x.1).unwrap();
                let (_, next) = clm
                    .iter()
                    .enumerate()
                    .filter(|x| x.0 != argmax)
                    .max_by_key(|x| x.1)
                    .unwrap();
                (next + 1) as f64 / (tot + 1) as f64
            })
        })
        .collect();
    fractions.sort_by(|x, y| x.partial_cmp(y).unwrap());
    fractions.reverse();
    // eprintln!("{:?}", &fractions[..10]);
    let thr = fractions[fractions.len() / 1000];
    let summary: HashMap<(usize, usize), u8> = summary
        .into_iter()
        .enumerate()
        .flat_map(|(idx, smry)| {
            smry.into_iter()
                .enumerate()
                .filter_map(|(pos, column)| {
                    let (argmax, _) = column.iter().enumerate().max_by_key(|x| x.1).unwrap();
                    let (_, next) = column
                        .iter()
                        .enumerate()
                        .filter(|x| x.0 != argmax)
                        .max_by_key(|x| x.1)
                        .unwrap();
                    let tot = column.iter().sum::<u32>();
                    let frac = (next + 1) as f64 / (tot + 1) as f64;
                    let base = match argmax {
                        0 => b'A',
                        1 => b'C',
                        2 => b'G',
                        3 => b'T',
                        _ => b'-',
                    };
                    if frac > thr {
                        Some(((idx, pos), base))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();
    maf.iter()
        .map(|aln| {
            summary
                .iter()
                .map(|(&(chunk, pos), mj)| match aln.get(&chunk) {
                    Some(aln) if aln[pos] == *mj => 1,
                    _ => 0,
                })
                .collect()
        })
        .collect()
}

fn alignment<F: bio::alignment::pairwise::MatchFunc>(
    aligner: &mut bio::alignment::pairwise::Aligner<F>,
    query: &[u8],
    template: &[u8],
) -> Vec<u8> {
    let aln = aligner.global(query, template);
    let mut seq = vec![];
    let mut ins = vec![b'-'; template.len()];
    let (mut rpos, mut qpos) = (0, 0);
    let mut prev = None;
    for &op in aln.operations.iter() {
        use bio::alignment::AlignmentOperation::*;
        match op {
            Del => {
                seq.push(b'-');
                rpos += 1;
            }
            Ins => {
                if prev != Some(Ins) && rpos < ins.len() {
                    ins[rpos] = query[qpos];
                }
                qpos += 1;
            }
            Subst | Match => {
                seq.push(query[qpos]);
                rpos += 1;
                qpos += 1;
            }
            _ => panic!(),
        }
        prev = Some(op);
    }
    let mut result = vec![];
    for (x, y) in seq.into_iter().zip(ins) {
        result.push(x);
        result.push(y);
    }
    result
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
    let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
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

pub fn em_clustering(
    xs: &[Vec<usize>],
    cluster: usize,
    (seed, trial_number, cluster_num): (u64, usize, usize),
) -> (Vec<usize>, f64) {
    let seed = 11 * seed + cluster as u64;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut weights: Vec<_> = xs
        .iter()
        .map(|_| {
            let mut ws = vec![0.; cluster];
            ws[rng.gen::<usize>() % cluster] = 1.;
            ws
        })
        .collect();
    let mut lk = std::f64::NEG_INFINITY;
    let mut parameters: Vec<Vec<Vec<f64>>>;
    let mut fractions: Vec<f64>;
    loop {
        // Update parameters.
        let cluster_weights: Vec<_> = (0..cluster)
            .map(|cl| weights.iter().map(|ws| ws[cl]).sum::<f64>() + 0.001)
            .collect();
        fractions = cluster_weights
            .iter()
            .map(|w| w / xs.len() as f64)
            .collect();
        parameters = (0..cluster)
            .map(|cl| {
                (0..trial_number)
                    .map(|tr| {
                        (0..cluster_num)
                            .map(|t| {
                                let sum = xs
                                    .iter()
                                    .zip(weights.iter())
                                    .filter(|(x, _)| x[tr] == t)
                                    .map(|(_, ws)| ws[cl])
                                    .sum::<f64>();
                                //(sum + 1.) / (cluster_weights[cl] + k as f64)
                                (sum + 0.000001) / cluster_weights[cl]
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        // Update weights
        let mut next_lk = 0.;
        weights = xs
            .iter()
            .map(|x| {
                let log_prob: Vec<_> = (0..cluster)
                    .map(|cl| {
                        fractions[cl].ln()
                            + x.iter()
                                .zip(parameters[cl].iter())
                                .map(|(&x, ps)| ps[x].ln())
                                .sum::<f64>()
                    })
                    .collect();
                let tot = logsumexp(&log_prob);
                next_lk += tot;
                log_prob.iter().map(|p| (p - tot).exp()).collect()
            })
            .collect();
        // Check log-LK.
        assert!(next_lk + 0.00001 > lk, "{}=>{}", lk, next_lk);
        if next_lk - lk < 0.00001 {
            break;
        }
        lk = next_lk;
    }
    let num_parameters = cluster * (1 + trial_number * (cluster_num - 1)) - 1;
    let neg_aic = lk - num_parameters as f64;
    trace!("{}\t{:.2}\t{}\t{:.2}", cluster, lk, num_parameters, neg_aic);
    // let fs: Vec<_> = fractions.iter().map(|x| format!("{:.2}", x)).collect();
    // trace!("{}", fs.join(","));
    let assignments: Vec<_> = weights
        .iter()
        .filter_map(|ws| {
            ws.iter()
                .enumerate()
                .max_by(|x, y| x.1.partial_cmp(y.1).unwrap())
        })
        .map(|x| x.0)
        .collect();
    let mut cluster: HashMap<_, usize> = HashMap::new();
    for a in assignments.iter() {
        if !cluster.contains_key(a) {
            cluster.insert(*a, cluster.len());
        }
    }
    let assignments: Vec<_> = assignments.iter().map(|x| cluster[x]).collect();
    (assignments, neg_aic)
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn rand_index_test() {
        let pred = [0, 0, 0, 1, 1, 1];
        let answ = [0, 0, 1, 1, 2, 2];
        assert!((0.6666 - rand_index(&pred, &answ)).abs() < 0.0001);
    }
}
