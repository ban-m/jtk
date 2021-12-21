use definitions::DataSet;
use definitions::EncodedRead;
use dirichlet_fit::{GreedyOptimizer, Optimizer};
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{self, Distribution};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
// If the posterior prob is less than this value, it would be removed...
const WEIGHT_FILTER: f64 = 0.00001;
// If the likelihood does not improve by this value, stop iteration.
const THRESHOLD: f64 = 0.0000000001;
// Small  value, to avoid underflow in .ln().
const SMALL_VALUE: f64 = 0.00000000000000001f64;
pub trait DirichletMixtureCorrection {
    fn correct_clustering(&mut self, config: &ClusteringConfig);
}

impl DirichletMixtureCorrection for DataSet {
    fn correct_clustering(&mut self, config: &ClusteringConfig) {
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
    config: &ClusteringConfig,
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
    let (contexts, _up_units, _down_units) = convert_to_contexts(&reads, unit_id, config);
    let (mut new_clustering, _lk, new_k) = (1..=k)
        .map(|k| clustering(&contexts, unit_id, k, config))
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

fn convert_to_contexts(
    reads: &[&EncodedRead],
    unit: u64,
    config: &ClusteringConfig,
) -> (Vec<Context>, HashMap<u64, usize>, HashMap<u64, usize>) {
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    unit_counts.retain(|_, count| config.coverage_thr < *count);
    let mut upstream_units: HashMap<_, usize> = HashMap::new();
    let mut downstream_units: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for (i, _) in read
            .nodes
            .iter()
            .enumerate()
            .filter(|&(_, n)| n.unit == unit)
        {
            let (before, after) = read.nodes.split_at(i);
            // Never panics.
            let (target, after) = after.split_first().unwrap();
            if target.is_forward {
                for n in before
                    .iter()
                    .rev()
                    .filter(|n| unit_counts.contains_key(&n.unit))
                {
                    let len = upstream_units.len();
                    upstream_units.entry(n.unit).or_insert(len);
                }
                for n in after.iter().filter(|n| unit_counts.contains_key(&n.unit)) {
                    let len = downstream_units.len();
                    downstream_units.entry(n.unit).or_insert(len);
                }
            } else {
                for n in after.iter().filter(|n| unit_counts.contains_key(&n.unit)) {
                    let len = upstream_units.len();
                    upstream_units.entry(n.unit).or_insert(len);
                }
                for n in before
                    .iter()
                    .rev()
                    .filter(|n| unit_counts.contains_key(&n.unit))
                {
                    let len = downstream_units.len();
                    downstream_units.entry(n.unit).or_insert(len);
                }
            }
        }
    }
    let contexts: Vec<_> = {
        let mut buffer = vec![];
        for read in reads.iter() {
            buffer.extend(
                (0..read.nodes.len())
                    .filter(|&i| read.nodes[i].unit == unit)
                    .map(|i| Context::new(read, i, &upstream_units, &downstream_units)),
            )
        }
        buffer
    };
    (contexts, upstream_units, downstream_units)
}

fn clustering(
    contexts: &[Context],
    unit: u64,
    cluster_num: usize,
    config: &ClusteringConfig,
) -> (Vec<(u64, usize, Vec<f64>)>, f64, usize) {
    let seed = (unit + 1) * cluster_num as u64;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let take_len = match cluster_num {
        1 => 1,
        _ => config.repeat_num,
    };
    let (asn, likelihood) = (0..take_len)
        .map(|_| clustering_inner(&contexts, cluster_num, &mut rng))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    (asn, likelihood, cluster_num)
}

#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    repeat_num: usize,
    coverage_thr: usize,
}

impl std::default::Default for ClusteringConfig {
    fn default() -> Self {
        Self {
            repeat_num: 20,
            coverage_thr: 5,
        }
    }
}
impl ClusteringConfig {
    pub fn new(repeat_num: usize, coverage_thr: usize) -> Self {
        Self {
            repeat_num,
            coverage_thr,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Context {
    // The ID of the original encoded read.
    id: u64,
    // The original index of this context.
    index: usize,
    #[allow(dead_code)]
    unit: u64,
    // This is log-scaled.
    center: Vec<f64>,
    // Only forward/backward information is preversed, the
    // ordering is arbitrary.
    // LOG posterior probabilities.
    upstream: Vec<(usize, Vec<f64>)>,
    downstream: Vec<(usize, Vec<f64>)>,
}

impl std::fmt::Display for Context {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (unit, prob) in self.upstream.iter().rev() {
            writeln!(f, "{}\t{}", unit, vec2str(prob))?;
        }
        writeln!(f, "C\t{}", vec2str(&self.center))?;
        for (i, (unit, prob)) in self.downstream.iter().enumerate() {
            if i == self.downstream.len() - 1 {
                write!(f, "{}\t{}", unit, vec2str(prob))?;
            } else {
                writeln!(f, "{}\t{}", unit, vec2str(prob))?;
            }
        }
        Ok(())
    }
}

impl Context {
    pub fn with(
        id: u64,
        index: usize,
        unit: u64,
        center: Vec<f64>,
        upstream: Vec<(usize, Vec<f64>)>,
        downstream: Vec<(usize, Vec<f64>)>,
    ) -> Self {
        Self {
            id,
            index,
            unit,
            center,
            upstream,
            downstream,
        }
    }
    fn new(
        read: &EncodedRead,
        index: usize,
        up_map: &HashMap<u64, usize>,
        down_map: &HashMap<u64, usize>,
    ) -> Self {
        let (before, after) = read.nodes.split_at(index);
        let (target, after) = after.split_first().unwrap();
        let (unit, center) = (target.unit, target.posterior.clone());
        let (upstream, downstream) = match read.nodes[index].is_forward {
            true => {
                let upstream: Vec<_> = before
                    .iter()
                    .rev()
                    .filter_map(|n| up_map.get(&n.unit).map(|&i| (i, n.posterior.clone())))
                    .collect();
                let downstream: Vec<_> = after
                    .iter()
                    .filter_map(|n| down_map.get(&n.unit).map(|&i| (i, n.posterior.clone())))
                    .collect();
                (upstream, downstream)
            }
            false => {
                let upstream: Vec<_> = after
                    .iter()
                    .filter_map(|n| up_map.get(&n.unit).map(|&i| (i, n.posterior.clone())))
                    .collect();
                let downstream: Vec<_> = before
                    .iter()
                    .rev()
                    .filter_map(|n| down_map.get(&n.unit).map(|&i| (i, n.posterior.clone())))
                    .collect();
                (upstream, downstream)
            }
        };
        let id = read.id;
        Self {
            id,
            index,
            unit,
            center,
            upstream,
            downstream,
        }
    }
}

pub fn clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64) {
    // ID of this trial.
    let mut weights = if k == 1 {
        vec![vec![1f64]; contexts.len()]
    } else {
        let dir = rand_distr::Dirichlet::new(&vec![1.5f64; k]).unwrap();
        contexts.iter().map(|_| dir.sample(rng)).collect()
    };
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    for (weights, ctx) in weights.iter().zip(contexts.iter()) {
        trace!("{}\t{}", vec2str(&ctx.center), vec2str(weights));
    }
    let mut model = HMMixtureModel::new(contexts, &weights, k);
    let mut lk: f64 = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
    trace!("CORRECT\tLikelihood\t{}\t0\t{}", id, lk);
    for t in 1..20 {
        trace!("CORRECT\tModel\t{}\t{}\n{}", t, id, model);
        model.update(&mut weights, contexts, t);
        let next_lk = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
        trace!("CORRECT\tLikelihood\t{}\t{}\t{}", id, t, next_lk);
        if (next_lk - lk) < THRESHOLD {
            break;
        } else {
            lk = next_lk;
        }
    }
    trace!("CORRECT\tModel\t{}\n{}", id, model);
    let predictions: Vec<_> = contexts
        .iter()
        .zip(weights.into_iter())
        .map(|(ctx, weight)| (ctx.id, ctx.index, weight))
        .collect();
    trace!("CORRECT\tFinal\t{}\t{}\t{:.4}", id, k, lk);
    (predictions, lk)
}

struct HMMixtureModel {
    k: usize,
    // Number of upstream units.
    upstream_units: usize,
    // Number of dowstream units.
    downstream_units: usize,
    fractions: Vec<f64>,
    models: Vec<HMModel>,
}

impl std::fmt::Display for HMMixtureModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}\t{}\t{}",
            self.k, self.upstream_units, self.downstream_units
        )?;
        for (i, (fr, m)) in self.fractions.iter().zip(self.models.iter()).enumerate() {
            if i != self.k - 1 {
                writeln!(f, "Cluster({})\t{:.3}\n{}", i, fr, m)?;
            } else {
                write!(f, "Cluster({})\t{:.3}\n{}", i, fr, m)?;
            }
        }
        Ok(())
    }
}

impl HMMixtureModel {
    fn new(contexts: &[Context], weights: &[Vec<f64>], cluster_num: usize) -> Self {
        let fractions = sum_and_normalize(weights);
        let up_units = contexts
            .iter()
            .filter_map(|ctx| ctx.upstream.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap_or(0);
        let down_units = contexts
            .iter()
            .filter_map(|ctx| ctx.downstream.iter().map(|x| x.0 + 1).max())
            .max()
            .unwrap_or(0);
        let mut up_dims: Vec<_> = vec![0; up_units];
        let mut down_dims: Vec<_> = vec![0; down_units];
        for ctx in contexts.iter() {
            for up in ctx.upstream.iter() {
                up_dims[up.0] = up_dims[up.0].max(up.1.len());
            }
            for down in ctx.downstream.iter() {
                down_dims[down.0] = down_dims[down.0].max(down.1.len());
            }
        }
        let models: Vec<_> = (0..cluster_num)
            .map(|cl| HMModel::new(contexts, weights, cl, &up_dims, &down_dims))
            .collect();
        Self {
            k: cluster_num,
            upstream_units: up_units,
            downstream_units: down_units,
            fractions,
            models,
        }
    }
    fn get_likelihood(&self, context: &Context) -> f64 {
        let lks: Vec<_> = self
            .fractions
            .iter()
            .zip(self.models.iter())
            .map(|(f, m)| f.ln() + m.lk(context))
            .collect();
        logsumexp(&lks)
    }
    fn update(&mut self, weights: &mut [Vec<f64>], contexts: &[Context], iteration: usize) {
        // E-step
        // Posterior prob and alignment expectation.
        let (posterior_prob, align_expt) = self.e_step(contexts);
        for (i, post) in posterior_prob.iter().enumerate() {
            let post: Vec<_> = post.iter().map(|x| x.log10()).collect();
            trace!("{}\t{}", i, vec2str(&post));
        }
        // M-step
        self.fractions = sum_and_normalize(&posterior_prob);
        for (cl, model) in self.models.iter_mut().enumerate() {
            let (mut ws, mut ctces, mut alns) = (vec![], vec![], vec![]);
            for ((weight, ctx), align) in weights
                .iter()
                .zip(contexts.iter())
                .zip(align_expt[cl].iter())
            // .filter(|&((ws, _), _)| WEIGHT_FILTER < ws[cl])
            {
                ws.push(weight[cl]);
                ctces.push(ctx);
                alns.push(align);
            }
            assert_eq!(ws.len(), ctces.len());
            assert_eq!(ws.len(), alns.len());
            model.update(&ws, &ctces, &alns, iteration);
        }
        for (ws, posterior) in weights.iter_mut().zip(posterior_prob) {
            *ws = posterior;
        }
    }
    // 2nd return value: the expectation of aligning information of each cluster k.
    fn e_step(&self, contexts: &[Context]) -> (Vec<Vec<f64>>, Vec<Vec<AlignInfo>>) {
        let alignments: Vec<Vec<_>> = self
            .models
            .iter()
            .map(|m| contexts.iter().map(|ctx| m.align(ctx)).collect())
            .collect();
        let posterior: Vec<_> = (0..contexts.len())
            .map(|idx| {
                let mut lks: Vec<_> = self
                    .fractions
                    .iter()
                    .zip(alignments.iter())
                    .map(|(f, alns)| f.ln() + alns[idx].lk)
                    .collect();
                let total = logsumexp(&lks);
                lks.iter_mut().for_each(|x| *x = (*x - total).exp());
                lks
            })
            .collect();
        (posterior, alignments)
    }
}

fn vec2str(xs: &[f64]) -> String {
    let xs: Vec<_> = xs.iter().map(|x| format!("{:.2}", x)).collect();
    xs.join(",")
}

struct AlignInfo {
    // Likelihood. Sum of center, upstream and downstream.
    lk: f64,
    // For each j-th position of the reference and the i-th position of the query,
    // upstream[j][i] return the probability that the j-th pos of the refr emit the i-th pos of the query.
    // Note that the last element of each row(upstream[j].last()) would be the prob of deletion.
    upstream: Vec<Vec<f64>>,
    downstream: Vec<Vec<f64>>,
}

// Sum xs and normalize the result so that it summing up to 1.
fn sum_and_normalize(xss: &[Vec<f64>]) -> Vec<f64> {
    let mut sumed = vec![SMALL_VALUE; xss[0].len()];
    let mut total = 0f64;
    for xs in xss {
        sumed
            .iter_mut()
            .zip(xs.iter())
            .for_each(|(acc, x)| *acc += x);
        total += xs.iter().sum::<f64>();
    }
    sumed.iter_mut().for_each(|x| *x /= total);
    sumed
}

impl std::fmt::Display for HMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for up in self.upstream.iter() {
            writeln!(f, "{}", up)?;
        }
        writeln!(f, "Center\t{}", self.center)?;
        for (i, down) in self.downstream.iter().enumerate() {
            if i != self.downstream_len - 1 {
                writeln!(f, "{}", down)?;
            } else {
                write!(f, "{}", down)?;
            }
        }
        Ok(())
    }
}

struct HMModel {
    upstream_len: usize,
    downstream_len: usize,
    center: Dirichlet,
    // the f64 parts summing up to 1 on each position.
    upstream: Vec<DirichletMixture>,
    // the f64 parts summing up to 1 on each position.
    downstream: Vec<DirichletMixture>,
}

const NEG_LARGE: f64 = -10000000000000f64;
type DP = Vec<Vec<f64>>;
impl HMModel {
    fn new(
        contexts: &[Context],
        weights: &[Vec<f64>],
        cluster: usize,
        up: &[usize],
        down: &[usize],
    ) -> Self {
        let upstream_len: usize = contexts
            .iter()
            .map(|ctx| ctx.upstream.len())
            .max()
            .unwrap_or(0);
        let downstream_len: usize = contexts
            .iter()
            .map(|ctx| ctx.downstream.len())
            .max()
            .unwrap_or(0);
        let weights: Vec<_> = weights.iter().map(|ws| ws[cluster]).collect();
        let center: Vec<_> = contexts.iter().map(|ctx| ctx.center.as_slice()).collect();
        let center = Dirichlet::new(&center, &weights, center[0].len());
        let mut upstream_slots = vec![vec![]; upstream_len];
        for (ctx, &w) in contexts.iter().zip(weights.iter()) {
            for (&(unit, ref prob), slot) in ctx.upstream.iter().zip(upstream_slots.iter_mut()) {
                slot.push((unit, w, prob.as_slice()));
            }
        }
        let upstream: Vec<_> = upstream_slots
            .iter()
            .map(|bucket| DirichletMixture::new(bucket, up))
            .collect();
        let mut downstream_slots = vec![vec![]; downstream_len];
        for (ctx, &w) in contexts.iter().zip(weights.iter()) {
            for (&(unit, ref prob), slot) in ctx.downstream.iter().zip(downstream_slots.iter_mut())
            {
                slot.push((unit, w, prob.as_slice()));
            }
        }
        let downstream: Vec<_> = downstream_slots
            .iter()
            .map(|bucket| DirichletMixture::new(bucket, down))
            .collect();
        Self {
            upstream_len,
            downstream_len,
            center,
            upstream,
            downstream,
        }
    }
    fn lk(&self, ctx: &Context) -> f64 {
        let up_forward = Self::forward(&self.upstream, &ctx.upstream);
        let up_lk = up_forward[ctx.upstream.len()][self.upstream_len];
        let down_forward = Self::forward(&self.downstream, &ctx.downstream);
        let down_lk = down_forward[ctx.downstream.len()][self.downstream_len];
        let center = self.center.lk(&ctx.center);
        up_lk + center + down_lk
    }
    fn align(&self, ctx: &Context) -> AlignInfo {
        let up_forward = Self::forward(&self.upstream, &ctx.upstream);
        let up_backward = Self::backward(&self.upstream, &ctx.upstream);
        let down_forward = Self::forward(&self.downstream, &ctx.downstream);
        let down_backward = Self::backward(&self.downstream, &ctx.downstream);
        let up_lk = up_forward[ctx.upstream.len()][self.upstream_len];
        let down_lk = down_forward[ctx.downstream.len()][self.downstream_len];
        let center = self.center.lk(&ctx.center);
        let upstream = Self::to_align(&self.upstream, &ctx.upstream, &up_forward, &up_backward);
        let downstream = Self::to_align(
            &self.downstream,
            &ctx.downstream,
            &down_forward,
            &down_backward,
        );
        AlignInfo {
            lk: up_lk + center + down_lk,
            upstream,
            downstream,
        }
    }
    fn forward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)]) -> DP {
        // Init
        let mut dp = vec![vec![0f64; refr.len() + 1]; query.len() + 1];
        for i in 1..query.len() + 1 {
            dp[i][0] = NEG_LARGE;
        }
        for j in 1..refr.len() + 1 {
            dp[0][j] = dp[0][j - 1] + refr[j - 1].del();
        }
        for (i, obs) in query.iter().enumerate().map(|(i, x)| (i + 1, x)) {
            for (j, dir) in refr.iter().enumerate().map(|(j, x)| (j + 1, x)) {
                let del_trans = dp[i][j - 1] + dir.del();
                let match_trans = dp[i - 1][j - 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans);
            }
        }
        dp
    }
    fn backward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)]) -> DP {
        // Init
        let mut dp = vec![vec![0f64; refr.len() + 1]; query.len() + 1];
        for i in 0..query.len() + 1 {
            dp[i][refr.len()] = NEG_LARGE;
        }
        for j in 0..refr.len() + 1 {
            dp[query.len()][j] = 0f64;
        }
        for (i, obs) in query.iter().enumerate().rev() {
            for (j, dir) in refr.iter().enumerate().rev() {
                let del_trans = dp[i][j + 1] + dir.del();
                let match_trans = dp[i + 1][j + 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans);
            }
        }
        dp
    }
    fn to_align(
        refr: &[DirichletMixture],
        query: &[(usize, Vec<f64>)],
        forward: &DP,
        backward: &DP,
    ) -> Vec<Vec<f64>> {
        refr.iter()
            .enumerate()
            .map(|(j, dir)| {
                let mut lks: Vec<_> = query
                    .iter()
                    .enumerate()
                    .map(|(i, obs)| forward[i][j] + dir.mat(obs) + backward[i + 1][j + 1])
                    .collect();
                // Push deletion
                {
                    let del_probs: Vec<_> = (0..query.len())
                        .map(|i| forward[i][j] + dir.del() + backward[i][j + 1])
                        .collect();
                    lks.push(logsumexp(&del_probs));
                }
                let total = logsumexp(&lks);
                lks.iter_mut().for_each(|x| *x = (*x - total).exp());
                lks
            })
            .collect()
    }
    fn update(
        &mut self,
        weights: &[f64],
        contexts: &[&Context],
        align_expt: &[&AlignInfo],
        iteration: usize,
    ) {
        let center: Vec<_> = contexts.iter().map(|ctx| ctx.center.as_slice()).collect();
        self.center.update(&center, weights, iteration);
        // Update upstream
        let up_align_expectation: Vec<_> = align_expt
            .iter()
            .map(|aln| aln.upstream.as_slice())
            .collect();
        let up_contexts: Vec<_> = contexts.iter().map(|ctx| ctx.upstream.as_slice()).collect();
        Self::update_oneside(
            &mut self.upstream,
            weights,
            &up_contexts,
            &up_align_expectation,
            iteration,
        );
        let down_align_expectation: Vec<&[Vec<f64>]> = align_expt
            .iter()
            .map(|aln| aln.downstream.as_slice())
            .collect();
        let down_contexts: Vec<_> = contexts
            .iter()
            .map(|ctx| ctx.downstream.as_slice())
            .collect();
        Self::update_oneside(
            &mut self.downstream,
            weights,
            &down_contexts,
            &down_align_expectation,
            iteration,
        );
    }
    fn update_oneside(
        parameters: &mut [DirichletMixture],
        weights: &[f64],
        contexts: &[&[(usize, Vec<f64>)]],
        alignments: &[&[Vec<f64>]],
        _iteration: usize,
    ) {
        for (j, dir) in parameters.iter_mut().enumerate() {
            let alns = alignments.iter().map(|aln| aln[j].iter());
            let total_ws: f64 = alns
                .clone()
                .zip(weights.iter())
                .map(|(alns, w)| w * alns.sum::<f64>())
                .sum();
            // It never panics.
            let del_ws: f64 = alns
                .clone()
                .zip(weights.iter())
                .map(|(alns, w)| w * alns.last().unwrap())
                .sum();
            dir.del_prob = (del_ws / total_ws).max(SMALL_VALUE);
            // First re-estimate parameters on categorical distributions.
            let mut fraction = vec![0f64; dir.dirichlets.len()];
            for ((w, ctx), alns) in weights.iter().zip(contexts.iter()).zip(alns.clone()) {
                for (&(unit, _), weight) in ctx.iter().zip(alns) {
                    fraction[unit] += w * weight;
                }
            }
            let total: f64 = fraction.iter().sum();
            dir.dirichlets
                .iter_mut()
                .zip(fraction)
                .for_each(|(x, f)| x.0 = f / total);
            // Then, re-estimate parameters on dirichlet distirbution on high-frequent unit.
            for (target, (_, dir)) in dir.dirichlets.iter_mut().enumerate()
            //.filter(|(_, (fr, dir))| WEIGHT_FILTER < *fr && 1 < dir.dim)
            {
                if dir.dim <= 1 {
                    continue;
                }
                let (weights, dataset): (Vec<_>, Vec<_>) = weights
                    .iter()
                    .zip(contexts.iter())
                    .zip(alns.clone())
                    .filter_map(|((&w, ctx), alns)| {
                        let (ws, obss): (Vec<f64>, Vec<_>) = ctx
                            .iter()
                            .zip(alns)
                            .filter_map(|(&(unit, ref prob), &weight)| {
                                (unit == target).then(|| (weight, prob.as_slice()))
                            })
                            .unzip();
                        (!ws.is_empty()).then(|| ((w, ws), obss))
                    })
                    .unzip();
                assert!(dir.dim > 1);
                trace!("MSTEP\t{}\t{}\t{}", j, target, weights.len());
                if !weights.is_empty() {
                    let mut opt = dirichlet_fit::GreedyOptimizer::new(dir.dim).set_norm(DIR_NORM);
                    dir.param =
                        dirichlet_fit::fit_multiple_with(&dataset, &weights, &mut opt, &dir.param);
                }
            }
        }
    }
}

struct DirichletMixture {
    del_prob: f64,
    dirichlets: Vec<(f64, Dirichlet)>,
}

impl std::fmt::Display for DirichletMixture {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dump: Vec<_> = self
            .dirichlets
            .iter()
            .enumerate()
            .filter(|&(_, &(fr, _))| fr > 0.1)
            .map(|(unit, (fr, dir))| format!("{}:{:.3}:{}", unit, fr, dir))
            .collect();
        write!(f, "{:.2}\t{}", self.del_prob, dump.join("\t"))
    }
}

impl DirichletMixture {
    // Return log del_prob
    fn del(&self) -> f64 {
        self.del_prob.ln()
    }
    // Return log (1-e)Dir(prob|param).
    fn mat(&self, &(unit, ref prob): &(usize, Vec<f64>)) -> f64 {
        let &(frac, ref dir) = &self.dirichlets[unit];
        (1f64 - self.del_prob).ln() + frac.ln() + dir.lk(prob)
    }
    // Create new instance from the given data.
    fn new(observed: &[(usize, f64, &[f64])], dims: &[usize]) -> Self {
        let del_prob = 0.1;
        let mut post_probs = vec![vec![]; dims.len()];
        let mut prob_weights = vec![vec![]; dims.len()];
        // TODO: Maybe 1 would be better?
        let mut fractions = vec![0f64; dims.len()];
        for &(unit, w, prob) in observed {
            fractions[unit] += w;
            post_probs[unit].push(prob);
            prob_weights[unit].push(w);
        }
        let sum: f64 = fractions.iter().sum();
        let dirichlets: Vec<_> = post_probs
            .iter()
            .zip(prob_weights.iter())
            .zip(fractions.iter())
            .zip(dims.iter())
            .map(|(((probs, ws), fr), &dim)| {
                assert_eq!(probs.len(), ws.len());
                (fr / sum, Dirichlet::new(probs, ws, dim))
            })
            .collect();
        Self {
            del_prob,
            dirichlets,
        }
    }
}

struct Dirichlet {
    dim: usize,
    param: Vec<f64>,
}

const DIR_NORM: f64 = 2f64;
impl Dirichlet {
    fn lk(&self, log_prob: &[f64]) -> f64 {
        assert_eq!(self.dim, log_prob.len());
        match self.dim {
            1 => 0f64,
            _ => dirichlet_fit::dirichlet_log(log_prob, &self.param),
        }
    }
    fn update<T: std::borrow::Borrow<[f64]>>(&mut self, obs: &[T], weights: &[f64], _t: usize) {
        if 1 < self.dim {
            let (data, weights) = (&[obs], &[(1f64, weights)]);
            // let mut optimizer = dirichlet_fit::AdamOptimizer::new(self.dim)
            //     .norm(DIR_NORM)
            //     .learning_rate(0.01 / t as f64);
            let mut optimizer = dirichlet_fit::GreedyOptimizer::new(self.dim).set_norm(DIR_NORM);
            self.param =
                dirichlet_fit::fit_multiple_with(data, weights, &mut optimizer, &self.param);
        }
    }
    fn new<T: std::borrow::Borrow<[f64]>>(obs: &[T], weights: &[f64], dim: usize) -> Self {
        // If there is no observation, return uniformal distribution.
        // let param = vec![DIR_NORM * (dim as f64).recip().sqrt(); dim];
        // let mut optimizer = AdamOptimizer::new(dim).norm(DIR_NORM);
        let mut param = vec![DIR_NORM * (dim as f64).recip(); dim];
        if !obs.is_empty() && !weights.is_empty() {
            let mut optimizer = GreedyOptimizer::new(dim).set_norm(DIR_NORM);
            let weights = &[(1f64, weights)];
            param = match dim {
                1 => vec![DIR_NORM],
                _ => dirichlet_fit::fit_multiple_with(&[obs], weights, &mut optimizer, &param),
            };
        }
        Self { dim, param }
    }
}

impl std::fmt::Display for Dirichlet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[")?;
        for (i, p) in self.param.iter().enumerate() {
            if i != self.dim - 1 {
                write!(f, "{:.2},", p)?;
            } else {
                write!(f, "{:.2}", p)?;
            }
        }
        write!(f, "]")
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        0.
    } else if xs.len() == 1 {
        xs[0]
    } else {
        let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
        assert!(sum >= 0., "{:?}->{}", xs, sum);
        max + sum
    }
}
// Log(exp(x) + exp(y)) = x + Log(1+exp(y-x)).
fn logsumexp2(x: f64, y: f64) -> f64 {
    let (x, y) = (x.max(y), x.min(y));
    x + (1f64 + (y - x).exp()).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn lse() {
        let x = 10.3f64;
        let y = 3.1f64;
        let test = logsumexp2(x, y);
        let answer = (x.exp() + y.exp()).ln();
        assert!((test - answer).abs() < 0.0001);
    }
    // (id, [(unit, posterior distr)]
    type Haplotype = (Vec<(usize, usize)>, usize, Vec<(usize, usize)>);
    type ClusterNum = (Vec<usize>, usize, Vec<usize>);
    type IsError = (Vec<bool>, bool, Vec<bool>);
    type UnitArray = Vec<(usize, Vec<f64>)>;
    type SimRead = (UnitArray, Vec<f64>, UnitArray);
    fn sim_read<R: Rng>(
        (hap_up, hap, hap_down): &Haplotype,
        (up_cl_num, cl_num, down_cl_num): &ClusterNum,
        (up_err, err, down_err): &IsError,
        drop: f64,
        rng: &mut R,
    ) -> SimRead {
        let mut upstream = vec![];
        for &(unit, cl) in hap_up.iter() {
            if rng.gen_bool(drop) {
                break;
            }
            let (dim, err) = (up_cl_num[unit], up_err[unit]);
            upstream.push((unit, gen_post_dist(cl, dim, err, rng)));
        }
        let center = gen_post_dist(*hap, *cl_num, *err, rng);
        let mut downstream = vec![];
        for &(unit, cl) in hap_down.iter() {
            if rng.gen_bool(drop) {
                break;
            }
            let (dim, err) = (down_cl_num[unit], down_err[unit]);
            downstream.push((unit, gen_post_dist(cl, dim, err, rng)));
        }
        (upstream, center, downstream)
    }
    fn gen_post_dist<R: Rng>(cl: usize, dim: usize, err: bool, rng: &mut R) -> Vec<f64> {
        if dim == 1 {
            return vec![0f64];
        }
        let mut param = vec![0.5f64; dim];
        if !err {
            param[cl] = 20.0;
        }
        let dirichlet = rand_distr::Dirichlet::new(&param).unwrap();
        dirichlet.sample(rng).into_iter().map(|x| x.ln()).collect()
    }

    fn gen_reads<R: Rng>(
        hap1: &Haplotype,
        hap2: &Haplotype,
        clnum: &ClusterNum,
        errors: &IsError,
        rng: &mut R,
        num: usize,
    ) -> Vec<Context> {
        (0..2 * num)
            .map(|i| match i / num {
                0 => (i, sim_read(hap1, clnum, errors, 0f64, rng)),
                1 => (i, sim_read(hap2, clnum, errors, 0f64, rng)),
                _ => panic!(),
            })
            .map(|(id, (upstream, center, downstream))| Context {
                id: id as u64,
                index: 0,
                unit: 0,
                center,
                upstream,
                downstream,
            })
            .collect()
    }
    fn gen_reads_drop<R: Rng>(
        haps: &[Haplotype],
        clnum: &ClusterNum,
        errors: &IsError,
        rng: &mut R,
        drop: f64,
        num: usize,
    ) -> Vec<Context> {
        (0..haps.len() * num)
            .map(|i| {
                let idx = i / num;
                (i, sim_read(&haps[idx], clnum, errors, drop, rng))
            })
            .map(|(id, (upstream, center, downstream))| Context {
                id: id as u64,
                index: 0,
                unit: 0,
                center,
                upstream,
                downstream,
            })
            .collect()
    }
    #[test]
    fn easy() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = (vec![], 0, vec![]);
        let hap2s = (vec![], 1, vec![]);
        let cl_num = { (vec![], 2, vec![]) };
        let is_uninformative = (vec![], false, vec![]);
        let num = 10;
        let reads = gen_reads(&hap1s, &hap2s, &cl_num, &is_uninformative, &mut rng, num);
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
    #[test]
    fn easy2() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = {
            let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
            let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
            (hap1_up, 0, hap1_down)
        };
        let hap2s = {
            let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
            let hap2_down = vec![(0, 1), (1, 1), (2, 1)];
            (hap2_up, 0, hap2_down)
        };
        let cl_num = { (vec![2, 2, 2], 2, vec![2, 2, 2]) };
        let is_uninformative = { (vec![false, false, false], false, vec![false, false, false]) };
        let num = 10;
        let reads = gen_reads(&hap1s, &hap2s, &cl_num, &is_uninformative, &mut rng, num);
        for (i, read) in reads.iter().enumerate() {
            trace!("{}\t{}\t{}", i, read.upstream.len(), read.downstream.len());
        }
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
    #[test]
    fn easy3() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = {
            let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
            let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
            (hap1_up, 0, hap1_down)
        };
        let hap2s = {
            let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
            let hap2_down = vec![(0, 1), (1, 1), (2, 1)];
            (hap2_up, 0, hap2_down)
        };
        let haps = vec![hap1s, hap2s];
        let drop = 0.15;
        let cl_num = { (vec![2, 2, 2], 2, vec![2, 2, 2]) };
        let is_uninformative = { (vec![false, false, false], false, vec![false, false, false]) };
        let num = 10;

        let reads = gen_reads_drop(&haps, &cl_num, &is_uninformative, &mut rng, drop, num);
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
    #[test]
    fn easy4() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = {
            let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
            let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
            (hap1_up, 0, hap1_down)
        };
        let hap2s = {
            let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
            let hap2_down = vec![(0, 1), (1, 1), (2, 1)];
            (hap2_up, 0, hap2_down)
        };
        let haps = vec![hap1s, hap2s];
        let drop = 0.15;
        let cl_num = { (vec![2, 2, 2], 2, vec![2, 2, 2]) };
        let is_uninformative = { (vec![false, false, false], true, vec![false, false, false]) };
        let num = 10;

        let reads = gen_reads_drop(&haps, &cl_num, &is_uninformative, &mut rng, drop, num);
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
    #[test]
    fn easy5() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = {
            let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
            let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
            (hap1_up, 0, hap1_down)
        };
        let hap2s = {
            let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
            let hap2_down = vec![(3, 1), (4, 1), (5, 1)];
            (hap2_up, 0, hap2_down)
        };
        let haps = vec![hap1s, hap2s];
        let drop = 0.15;
        let cl_num = { (vec![2, 2, 2], 2, vec![2, 2, 2, 2, 2, 2]) };
        let is_uninformative = { (vec![false; 3], false, vec![false; 6]) };
        let num = 10;
        let reads = gen_reads_drop(&haps, &cl_num, &is_uninformative, &mut rng, drop, num);
        for ctx in reads.iter() {
            println!("{}\n", ctx);
        }
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        for (_, _, res) in result.iter() {
            println!("{}", vec2str(res));
        }
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
    #[test]
    fn hard1() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = {
            let hap1_up = vec![(0, 0), (1, 0), (2, 0)];
            let hap1_down = vec![(0, 0), (1, 0), (2, 0)];
            (hap1_up, 0, hap1_down)
        };
        let hap2s = {
            let hap2_up = vec![(0, 1), (1, 1), (2, 1)];
            let hap2_down = vec![(3, 1), (4, 1), (5, 1)];
            (hap2_up, 0, hap2_down)
        };
        let haps = vec![hap1s, hap2s];
        let drop = 0.20;
        let cl_num = { (vec![2, 2, 2], 3, vec![2, 1, 2, 1, 2, 2]) };
        let is_uninformative = { (vec![false; 3], true, vec![false; 6]) };
        let num = 10;
        let reads = gen_reads_drop(&haps, &cl_num, &is_uninformative, &mut rng, drop, num);
        for ctx in reads.iter() {
            println!("{}\n", ctx);
        }
        let (result, _) = clustering_inner(&reads, 2, &mut rng);
        for (_, _, res) in result.iter() {
            println!("{}", vec2str(res));
        }
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
}
