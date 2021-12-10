//! Read clustering module.
use definitions::{DataSet, EncodedRead};
use rayon::prelude::*;
use std::collections::HashMap;
pub trait ReadClustering {
    fn read_clustering(&mut self, config: &ReadClusteringConfig);
}

#[derive(Debug, Clone)]
pub struct ReadClusteringConfig {
    repeat_num: usize,
    coverage_thr: usize,
}

impl std::default::Default for ReadClusteringConfig {
    fn default() -> Self {
        Self {
            repeat_num: 20,
            coverage_thr: 4,
        }
    }
}

impl ReadClusteringConfig {
    pub fn new(repeat_num: usize, coverage_thr: usize) -> Self {
        Self {
            repeat_num,
            coverage_thr,
        }
    }
}

impl ReadClustering for DataSet {
    fn read_clustering(&mut self, config: &ReadClusteringConfig) {
        for chunk in self.selected_chunks.iter() {
            assert!(chunk.cluster_num > 0, "{},{}", chunk.id, chunk.cluster_num);
        }
        let cl_num: HashMap<u64, usize> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            assert_eq!(node.posterior.len(), cl_num[&node.unit]);
        }
        let posterior_distributions: Vec<_> = self
            .selected_chunks
            .par_iter()
            .map(|c| (c.id, c.cluster_num))
            .map(|(id, cluster_num)| correct_unit(self, id, cluster_num, config))
            .collect();
        let mut result: HashMap<u64, Vec<(usize, Vec<f64>)>> = {
            let mut result: HashMap<u64, Vec<(usize, Vec<f64>)>> = HashMap::new();
            for post_dist_on_chunk in posterior_distributions {
                for (id, pos, posterior_dist) in post_dist_on_chunk {
                    result.entry(id).or_default().push((pos, posterior_dist));
                }
            }
            result
        };
        let argmax = |xs: &[f64]| {
            xs.iter()
                .enumerate()
                .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
                .unwrap()
                .0
        };
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = result.remove(&read.id) {
                for (pos, post) in corrected {
                    read.nodes[pos].cluster = argmax(&post) as u64;
                }
            }
        }
    }
}

pub fn correct_unit(
    ds: &DataSet,
    unit_id: u64,
    k: usize,
    config: &ReadClusteringConfig,
) -> Vec<(u64, usize, Vec<f64>)> {
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
        .collect();
    if reads.is_empty() {
        return vec![];
    }
    let start = std::time::Instant::now();
    let cov = ds.coverage.unwrap_or(10f64);
    let (mut new_clustering, _lk, new_k) = (1..=k)
        .map(|k| clustering(&reads, unit_id, k, cov, config))
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    let pad_len = k.saturating_sub(new_k);
    for (_, _, prob) in new_clustering.iter_mut() {
        prob.extend(std::iter::repeat(0f64).take(pad_len));
    }
    let end = std::time::Instant::now();
    trace!("ReadClustering\t{}\t{}", unit_id, (end - start).as_secs());
    new_clustering
}

use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

fn clustering(
    reads: &[&EncodedRead],
    unit: u64,
    cluster_num: usize,
    cov: f64,
    config: &ReadClusteringConfig,
) -> (Vec<(u64, usize, Vec<f64>)>, f64, usize) {
    let seed = (unit + 1) * cluster_num as u64;
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    unit_counts.retain(|_, count| config.coverage_thr < *count);
    // TODO:Here maybe we can serialize the units so that it follow numerical order, 0,1,2,...
    let contexts: Vec<_> = {
        let mut buffer = vec![];
        for read in reads.iter() {
            for index in 0..read.nodes.len() {
                if read.nodes[index].unit == unit {
                    buffer.push(Context::new(read, index, &unit_counts));
                }
            }
        }
        buffer
    };
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let take_len = match cluster_num {
        1 => 1,
        _ => config.repeat_num,
    };
    let (asn, likelihood) = (0..take_len)
        .map(|_| simple_clustering_inner(&contexts, cluster_num, cov, config, &mut rng))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    (asn, likelihood, cluster_num)
}

#[derive(Debug, Clone)]
struct Context<'a> {
    // The ID of the original encoded read.
    id: u64,
    // The original index of this context.
    index: usize,
    unit: u64,
    cluster: &'a [f64],
    // Only upstream/downstream information is preversed, the
    // ordering is arbitrary.
    upstream: Vec<(u64, &'a [f64])>,
    downstream: Vec<(u64, &'a [f64])>,
}

impl std::fmt::Display for Context<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let upstream: Vec<_> = self
            .upstream
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .collect();
        let downstream: Vec<_> = self
            .downstream
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .collect();
        write!(
            f,
            "{}\t({},{:?})\t{}",
            upstream.join("-"),
            self.unit,
            self.cluster,
            downstream.join("-")
        )
    }
}

impl<'a> Context<'a> {
    fn new(read: &'a EncodedRead, index: usize, unit_counts: &HashMap<u64, usize>) -> Self {
        let (unit, cluster) = (
            read.nodes[index].unit,
            read.nodes[index].posterior.as_slice(),
        );
        let nodes = read.nodes.iter();
        let upstream: Vec<_> = nodes
            .clone()
            .skip(index + 1)
            .map(|n| (n.unit, n.posterior.as_slice()))
            .filter(|n| unit_counts.contains_key(&n.0))
            .collect();
        let downstream: Vec<_> = nodes
            .clone()
            .take(index)
            .rev()
            .map(|n| (n.unit, n.posterior.as_slice()))
            .filter(|n| unit_counts.contains_key(&n.0))
            .collect();
        let id = read.id;
        match read.nodes[index].is_forward {
            true => Self {
                id,
                index,
                unit,
                cluster,
                upstream,
                downstream,
            },
            false => Self {
                id,
                index,
                unit,
                cluster,
                upstream: downstream,
                downstream: upstream,
            },
        }
    }
}

fn max_poisson_lk(obs: usize, coverage: f64) -> f64 {
    const ERROR_FRAC: f64 = 0.25;
    fn poisson(obs: usize, copy_num: usize, coverage: f64) -> f64 {
        let lambda = match copy_num {
            0 => coverage * ERROR_FRAC,
            _ => copy_num as f64 * coverage,
        };
        let denom: f64 = (1..obs + 1).map(|i| (i as f64).ln()).sum();
        (obs as f64 * lambda.ln() - lambda - denom).exp()
    }
    let copy_num = (obs as f64 / coverage).floor() as usize;
    (copy_num.saturating_sub(1)..=copy_num + 1)
        .map(|x| poisson(obs, x, coverage))
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap()
}

fn simple_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    cov: f64,
    _config: &ReadClusteringConfig,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64) {
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    // Current (un-modified) likelihoods.
    let mut partition = ReadPartition::new(&contexts, k, cov, rng);
    let total = 100_000 * k as usize;
    // MAP estimation.
    let (mut max, mut argmax) = (partition.lk(), partition.clone());
    let mut no_improve = 0;
    for t in 0..total {
        match partition.flip(rng) {
            Some(new_lk) if max < new_lk => {
                trace!("ReadPart\t{}\t{}\t{:.3}", id, t, new_lk);
                max = new_lk;
                argmax = partition.clone();
                no_improve = 0;
            }
            _ => no_improve += 1,
        }
        if 100 * contexts.len() * k < no_improve {
            break;
        }
    }
    let lk_gain = argmax.read_gains();
    (lk_gain, max)
}

#[derive(Debug, Clone)]
pub struct ReadPartition<'a> {
    // assignments
    assignments: Vec<usize>,
    // reads(contexts).
    reads: Vec<Context<'a>>,
    // Current likelihood bugs.
    likelihood_vectors: Vec<LikelihoodVector>,
    // Poisson likleihood(memoized)
    poisson_lks: Vec<f64>,
    // binomal likelihoos(memoized),
    // binom_lks: Vec<f64>,
    // # of cluster.
    cluster_num: usize,
}

impl<'a> std::fmt::Display for ReadPartition<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for lkv in self.likelihood_vectors.iter() {
            writeln!(f, "{}", lkv)?;
        }
        write!(f, "{:?}", self.assignments)
    }
}

impl<'a> ReadPartition<'a> {
    fn new<R: Rng>(contexts: &[Context<'a>], k: usize, cov: f64, rng: &mut R) -> Self {
        let assignments: Vec<_> = contexts.iter().map(|_| rng.gen_range(0..k)).collect();
        let max_occ = {
            let mut occs: HashMap<_, usize> = HashMap::new();
            for ctx in contexts.iter() {
                for (unit, _) in ctx.upstream.iter() {
                    *occs.entry(*unit).or_default() += 1;
                }
                for (unit, _) in ctx.downstream.iter() {
                    *occs.entry(*unit).or_default() += 1;
                }
            }
            (*occs.values().max().unwrap_or(&0)).max(contexts.len())
        };
        // let binom_lks: Vec<_> = (0..=max_occ)
        //     .map(|x| max_binom_lk(x, contexts.len(), cov / contexts.len() as f64))
        //     .collect();
        let poisson_lks: Vec<_> = (0..=max_occ).map(|x| max_poisson_lk(x, cov)).collect();
        let likelihood_vectors: Vec<_> = (0..k)
            .map(|cl| LikelihoodVector::new(&contexts, &assignments, cl))
            .collect();
        Self {
            assignments,
            reads: contexts.to_vec(),
            // binom_lks,
            poisson_lks,
            likelihood_vectors,
            cluster_num: k,
        }
    }
    fn lk(&self) -> f64 {
        self.likelihood_vectors
            .iter()
            .map(|lkv| lkv.lk(&self.poisson_lks))
            .sum()
    }
    fn flip<R: Rng>(&mut self, rng: &mut R) -> Option<f64> {
        let idx = rng.gen_range(0..self.reads.len());
        let new = rng.gen_range(0..self.cluster_num);
        let (old, ctx) = (self.assignments[idx], &self.reads[idx]);
        if old == new {
            return None;
        }
        let lk = self.lk();
        // Change the assignment. TODO:Maybe this code is unstable because it includes many sub/add operations.
        self.likelihood_vectors[old].sub(ctx);
        self.likelihood_vectors[new].add(ctx);
        // calc lk
        let new_lk = self.lk();
        let prob = (new_lk - lk).exp().min(1f64);
        if rng.gen_bool(prob) {
            // Success.
            self.assignments[idx] = new;
            Some(new_lk)
        } else {
            // Flip fail. Roll back the the previous state.
            self.likelihood_vectors[old].add(ctx);
            self.likelihood_vectors[new].sub(ctx);
            None
        }
    }
    fn read_gains(&self) -> Vec<(u64, usize, Vec<f64>)> {
        let used_positions: Vec<_> = self
            .likelihood_vectors
            .iter()
            .map(|lkv| lkv.used_positions())
            .collect();
        self.reads
            .iter()
            .map(|ctx| {
                let lks: Vec<f64> = used_positions
                    .iter()
                    .map(|(up, center, down)| {
                        // TODO: How to treamt no-occ units?
                        let up: f64 = ctx
                            .upstream
                            .iter()
                            .map(|(unit, gains)| match up.get(unit) {
                                Some(&idx) => gains[idx],
                                None => 0f64,
                            })
                            .sum();
                        let down: f64 = ctx
                            .downstream
                            .iter()
                            .map(|(unit, gains)| match down.get(unit) {
                                Some(&idx) => gains[idx],
                                None => 0f64,
                            })
                            .sum();
                        let center = ctx.cluster[*center];
                        up + center + down
                    })
                    .collect();
                (ctx.id, ctx.index, lks)
            })
            .collect()
    }
}

#[derive(Debug, Clone)]
struct LikelihoodVector {
    // Sum of the likelihood and its counts.
    center: Vec<f64>,
    count: usize,
    // TODO:They are too slow.
    upstream: HashMap<u64, (usize, Vec<f64>)>,
    downstream: HashMap<u64, (usize, Vec<f64>)>,
}

fn vec2str(xs: &[f64]) -> String {
    xs.iter().fold(String::new(), |x, y| {
        if x.is_empty() {
            format!("{:.3}", y)
        } else {
            x + &format!("\t{:.3}", y)
        }
    })
}

impl std::fmt::Display for LikelihoodVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (unit, (count, lks)) in self.upstream.iter() {
            writeln!(f, "{}\t{}\t{}", unit, count, vec2str(lks))?;
        }
        writeln!(f, "C\t{}\t{}", self.count, vec2str(&self.center))?;
        for (unit, (count, lks)) in self.downstream.iter() {
            writeln!(f, "{}\t{}\t{}", unit, count, vec2str(lks))?;
        }
        Ok(())
    }
}

impl LikelihoodVector {
    fn new(contexts: &[Context], assignments: &[usize], cluster: usize) -> Self {
        let mut center = vec![0f64; contexts[0].cluster.len()];
        let contexts = contexts
            .iter()
            .zip(assignments)
            .filter(|&(_, &asn)| asn == cluster);
        for (ctx, _) in contexts.clone() {
            for (x, y) in center.iter_mut().zip(ctx.cluster.iter()) {
                *x += y;
            }
        }
        let count = assignments.iter().filter(|&&asn| asn == cluster).count();
        let mut upstream: HashMap<u64, (usize, Vec<f64>)> = HashMap::new();
        let mut downstream: HashMap<u64, (usize, Vec<f64>)> = HashMap::new();
        for (ctx, _) in contexts {
            for &(unit, post) in ctx.upstream.iter() {
                let insert = || (0, vec![0f64; post.len()]);
                let (count, slot) = upstream.entry(unit).or_insert_with(insert);
                slot.iter_mut().zip(post).for_each(|(x, y)| *x += y);
                *count += 1;
            }
            for &(unit, post) in ctx.downstream.iter() {
                let insert = || (0, vec![0f64; post.len()]);
                let (count, slot) = downstream.entry(unit).or_insert_with(insert);
                slot.iter_mut().zip(post).for_each(|(x, y)| *x += y);
                *count += 1;
            }
        }
        Self {
            center,
            count,
            downstream,
            upstream,
        }
    }
    fn pos_max(xs: &[f64]) -> f64 {
        xs.iter().fold(0f64, |x, &y| x.max(y))
    }
    // TODO: How to treat zero-count units?
    fn lk(&self, count_lks: &[f64]) -> f64 {
        let center = Self::pos_max(&self.center);
        let upstream: f64 = self
            .upstream
            .values()
            .map(|&(count, ref lks)| Self::pos_max(lks) + count_lks[count])
            .sum();
        let downstream: f64 = self
            .downstream
            .values()
            .map(|&(count, ref lks)| Self::pos_max(lks) + count_lks[count])
            .sum();
        center + upstream + downstream
    }
    fn sub(&mut self, ctx: &Context) {
        self.center
            .iter_mut()
            .zip(ctx.cluster.iter())
            .for_each(|(x, y)| *x -= y);
        self.count -= 1;
        for (unit, post) in ctx.downstream.iter() {
            if let Some((count, slot)) = self.downstream.get_mut(unit) {
                slot.iter_mut().zip(post.iter()).for_each(|(x, y)| *x -= y);
                *count -= 1;
            }
        }
        for (unit, post) in ctx.upstream.iter() {
            if let Some((count, slot)) = self.upstream.get_mut(unit) {
                slot.iter_mut().zip(post.iter()).for_each(|(x, y)| *x -= y);
                *count -= 1;
            }
        }
    }
    fn add(&mut self, ctx: &Context) {
        self.center
            .iter_mut()
            .zip(ctx.cluster.iter())
            .for_each(|(x, y)| *x += y);
        self.count += 1;
        for (unit, post) in ctx.downstream.iter() {
            if let Some((count, slot)) = self.downstream.get_mut(unit) {
                slot.iter_mut().zip(post.iter()).for_each(|(x, y)| *x += y);
                *count += 1;
            }
        }
        for (unit, post) in ctx.upstream.iter() {
            if let Some((count, slot)) = self.upstream.get_mut(unit) {
                slot.iter_mut().zip(post.iter()).for_each(|(x, y)| *x += y);
                *count += 1;
            }
        }
    }
    // Upstream, center, donwstream
    fn used_positions(&self) -> (HashMap<u64, usize>, usize, HashMap<u64, usize>) {
        fn argmax(lks: &[f64]) -> usize {
            lks.iter()
                .enumerate()
                .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
                .unwrap()
                .0
        }
        let upstream: HashMap<_, _> = self
            .upstream
            .iter()
            .filter(|&(_, &(count, _))| 0 < count)
            .map(|(&unit, (_, lks))| (unit, argmax(lks)))
            .collect();
        let center = argmax(&self.center);
        let downstream: HashMap<_, _> = self
            .downstream
            .iter()
            .filter(|&(_, &(count, _))| 0 < count)
            .map(|(&unit, (_, lks))| (unit, argmax(lks)))
            .collect();
        (upstream, center, downstream)
    }
}

#[cfg(test)]
pub mod tests {
    impl<'a> Context<'a> {
        fn with_attrs(
            id: u64,
            index: usize,
            unit: u64,
            cluster: &'a [f64],
            upstream: Vec<(u64, &'a [f64])>,
            downstream: Vec<(u64, &'a [f64])>,
        ) -> Self {
            Self {
                id,
                index,
                unit,
                cluster,
                upstream,
                downstream,
            }
        }
    }
    use super::*;
    #[test]
    fn test_correction() {
        let len = 10;
        let dataset: Vec<_> = (0..len)
            .map(|i| {
                let forward: Vec<(u64, Vec<f64>)> = vec![];
                let (backward, cluster) = if i < len / 2 {
                    (vec![(0, vec![1f64])], vec![1f64, 0f64])
                } else {
                    (vec![(1, vec![1f64])], vec![0f64, 1f64])
                };
                (i, cluster, forward, backward)
            })
            .collect();
        let contexts: Vec<_> = dataset
            .iter()
            .map(|(id, cluster, forward, backward)| {
                let forward: Vec<_> = forward.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                let backward: Vec<_> = backward.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                Context::with_attrs(*id, 0, 2, cluster.as_slice(), forward, backward)
            })
            .collect();
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(942830);
        let cov = 5f64;
        let config = ReadClusteringConfig::default();
        let (asn, _likelihood) = simple_clustering_inner(&contexts, 2, cov, &config, &mut rng);
        let num_correct = asn
            .iter()
            .filter(|&&(id, _, ref post)| {
                if id < len / 2 {
                    post[0] < post[1]
                } else {
                    post[1] < post[0]
                }
            })
            .count();
        let num_correct = num_correct as u64;
        assert!(num_correct.max(len - num_correct) > len * 8 / 10);
    }
    #[test]
    fn test_correction_2() {
        let len = 20;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(942830);
        let dataset: Vec<_> = (0..len)
            .map(|i| {
                let (upstream, center, downstream) = if i < len / 2 {
                    let up = vec![(3, vec![10f64, 0f64])];
                    let center = vec![0f64, rng.gen_range(-5f64..5f64)];
                    let down = vec![(5, vec![0f64, -10f64])];
                    (up, center, down)
                } else {
                    let up = vec![(3, vec![-10f64, 0f64])];
                    let center = vec![0f64, rng.gen_range(-5f64..5f64)];
                    let down = vec![(5, vec![0f64, 10f64])];
                    (up, center, down)
                };
                (i, center, upstream, downstream)
            })
            .collect();
        let contexts: Vec<_> = dataset
            .iter()
            .map(|(id, center, upstream, downstream)| {
                let upstream: Vec<_> = upstream.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                let downstream: Vec<_> =
                    downstream.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                Context::with_attrs(*id, 0, 2, center, upstream, downstream)
            })
            .collect();
        let cov = len as f64 / 2f64;
        let config = ReadClusteringConfig::default();
        let (asn, likelihood) = simple_clustering_inner(&contexts, 2, cov, &config, &mut rng);
        println!("{}", likelihood);
        let num_correct = asn
            .iter()
            .filter(|&&(id, _, ref post)| {
                println!("{}\t{:?}", id, post);
                if id < len / 2 {
                    post[0] < post[1]
                } else {
                    post[1] < post[0]
                }
            })
            .count();
        let num_correct = num_correct as u64;
        assert!(
            num_correct.max(len - num_correct) > len * 8 / 10,
            "{}",
            num_correct
        );
    }
    #[test]
    fn test_correction_del() {
        let len = 20;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(942830);
        let dataset: Vec<_> = (0..len)
            .map(|i| {
                let (upstream, center, downstream) = if i < len / 2 {
                    let up = vec![(3, vec![10f64, 0f64])];
                    let center = vec![0f64, rng.gen_range(-5f64..5f64)];
                    let down = vec![(5, vec![0f64])];
                    (up, center, down)
                } else {
                    let up = vec![(3, vec![-10f64, 0f64])];
                    let center = vec![0f64, rng.gen_range(-5f64..5f64)];
                    let down = vec![(6, vec![0f64])];
                    (up, center, down)
                };
                (i, center, upstream, downstream)
            })
            .collect();
        let contexts: Vec<_> = dataset
            .iter()
            .map(|(id, center, upstream, downstream)| {
                let upstream: Vec<_> = upstream.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                let downstream: Vec<_> =
                    downstream.iter().map(|(i, f)| (*i, f.as_slice())).collect();
                Context::with_attrs(*id, 0, 2, center, upstream, downstream)
            })
            .collect();
        let cov = len as f64 / 2f64;
        let config = ReadClusteringConfig::default();
        let (asn, likelihood) = simple_clustering_inner(&contexts, 2, cov, &config, &mut rng);
        println!("{}", likelihood);
        let num_correct = asn
            .iter()
            .filter(|&&(id, _, ref post)| {
                println!("{}\t{:?}", id, post);
                if id < len / 2 {
                    post[0] < post[1]
                } else {
                    post[1] < post[0]
                }
            })
            .count();
        let num_correct = num_correct as u64;
        assert!(
            num_correct.max(len - num_correct) > len * 8 / 10,
            "{}",
            num_correct
        );
    }
}
