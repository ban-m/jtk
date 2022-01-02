use definitions::DataSet;
use definitions::EncodedRead;
use dirichlet_fit::Dirichlet;
const DIR_NORM: Option<f64> = Some(2f64);
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{self, Distribution};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
const APPROX: bool = true;
const LOG_WEIGHT_FILTER: f64 = -30f64;
// If the likelihood does not improve by this value, stop iteration.
const THRESHOLD: f64 = 0.0001;
// If the clustering with different seeds does not different at least this value w.r.t likelihood,
// it is the same
const THRESHOLD_CL: f64 = 0.1f64;
// use std::println as trace;
// Small  value, to avoid underflow in .ln().
// shoule be smaller than LOG_WEIGHT_FILTER afte taking log.
const SMALL_VALUE: f64 = 0.00000000000000000000000000000000000000000000001f64;
use std::collections::HashSet;
pub trait DirichletMixtureCorrection {
    fn correct_clustering(&mut self, config: &ClusteringConfig);
    fn correct_clustering_on_selected(
        &mut self,
        config: &ClusteringConfig,
        selection: &HashSet<u64>,
    );
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
        // let (argmax, max) = cl_num.iter().max_by_key(|x| x.1).unwrap();
        let selections: Vec<_> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        let posterior_distributions: Vec<_> = selections
            .par_iter()
            .map(|&(id, cluster_num)| correct_unit(self, id, cluster_num, config))
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
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = result.remove(&read.id) {
                for (pos, post) in corrected {
                    read.nodes[pos].cluster = argmax(&post) as u64;
                }
            }
        }
    }

    fn correct_clustering_on_selected(
        &mut self,
        config: &ClusteringConfig,
        selection: &HashSet<u64>,
    ) {
        let selections: Vec<_> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .filter(|(id, _)| selection.contains(id))
            .collect();
        let posterior_distributions: Vec<_> = selections
            .par_iter()
            .map(|&(id, cluster_num)| correct_unit(self, id, cluster_num, config))
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
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = result.remove(&read.id) {
                for (pos, post) in corrected {
                    read.nodes[pos].cluster = argmax(&post) as u64;
                }
            }
        }
    }
}

fn argmax(xs: &[f64]) -> usize {
    xs.iter()
        .enumerate()
        .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
        .unwrap()
        .0
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
    let (contexts, up_units, down_units) = convert_to_contexts(&reads, unit_id, config);
    let tot = {
        let (mut up_clnum, mut down_clnum) = (HashMap::new(), HashMap::new());
        for ctx in contexts.iter() {
            for up in ctx.upstream.iter() {
                up_clnum.entry(up.0).or_insert(up.1.len());
            }
            for down in ctx.downstream.iter() {
                down_clnum.entry(down.0).or_insert(down.1.len());
            }
        }
        up_clnum.values().sum::<usize>()
            + down_clnum.values().sum::<usize>()
            + contexts[0].center.len()
    };
    debug!(
        "ReadClustering\t{}\tBEGIN\t{}\t{}\t{}\t{}",
        unit_id,
        contexts.len(),
        up_units.len(),
        down_units.len(),
        tot,
    );
    let init = match k {
        0..=3 => 1,
        _ => k - 3,
    };
    let (mut new_clustering, _lk, new_k) = (init..=k)
        .map(|k| clustering(&contexts, unit_id, k, config))
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    // debug!("ReadClustering\tFinal\t{}\t{}\t{:.4}", unit_id, new_k, lk);
    let pad_len = k.saturating_sub(new_k);
    for (_, _, prob) in new_clustering.iter_mut() {
        prob.extend(std::iter::repeat(0f64).take(pad_len));
    }
    let end = std::time::Instant::now();
    let duration = (end - start).as_secs();
    debug!("ReadClustering\t{}\tEND\t{}", unit_id, duration,);
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

type PosteriorAtLocation = (u64, usize, Vec<f64>);
fn clustering(
    contexts: &[Context],
    unit: u64,
    cluster_num: usize,
    config: &ClusteringConfig,
) -> (Vec<PosteriorAtLocation>, f64, usize) {
    let seed = (unit + 1) * cluster_num as u64;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let retry = config.retry_num;
    if cluster_num == 1 {
        let (asn, lk) = clustering_inner(contexts, cluster_num, retry, &mut rng);
        return (asn, lk, cluster_num);
    }
    let (mut asn, mut likelihood) = clustering_inner(contexts, cluster_num, retry, &mut rng);
    let mut last_upd = 0;
    let start = std::time::Instant::now();
    for loop_count in 0.. {
        let (asn_n, lk_n) = clustering_inner(contexts, cluster_num, retry, &mut rng);
        if THRESHOLD_CL < lk_n - likelihood {
            likelihood = lk_n;
            asn = asn_n;
            last_upd = loop_count;
        }
        if config.follow_through < loop_count - last_upd {
            break;
        }
    }
    let end = std::time::Instant::now();
    let duration = (end - start).as_secs();
    trace!(
        "CORRECTION\tMAX\t{}\t{}\t{}\t{}",
        unit,
        last_upd,
        cluster_num,
        duration
    );
    (asn, likelihood, cluster_num)
}

#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    // How many times we re-estimate parameters in EM-algorithm after the
    // log-likelihood reached a plateau.
    follow_through: usize,
    // How many times we re-execute whole EM-algorithm after reaching local maximum.
    retry_num: usize,
    coverage_thr: usize,
}

impl std::default::Default for ClusteringConfig {
    fn default() -> Self {
        Self {
            retry_num: 5,
            follow_through: 50,
            coverage_thr: 5,
        }
    }
}
impl ClusteringConfig {
    pub fn new(follow_through: usize, retry_num: usize, coverage_thr: usize) -> Self {
        Self {
            follow_through,
            retry_num,
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
    retry: usize,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64) {
    // ID of this trial.
    let start = std::time::Instant::now();
    let mut weights = if k == 1 {
        vec![vec![1f64]; contexts.len()]
    } else {
        let dir = rand_distr::Dirichlet::new(&vec![1.5f64; k]).unwrap();
        contexts.iter().map(|_| dir.sample(rng)).collect()
    };
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut model = HMMixtureModel::new(contexts, &weights, k);
    trace!("CORRECT\tModel\t{}\n{}", id, model);
    let mut lk: f64 = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
    trace!("CORRECT\tLikelihood\t{}\t0\t{}", id, lk);
    let mut count = 0;
    for t in 1.. {
        for _l in 0..retry {
            count += 1;
            model.update(&mut weights, contexts, t);
        }
        let next_lk = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
        trace!("CORRECT\tLikelihood\t{}\t{}\t{}", id, t, next_lk,);
        if next_lk - lk < THRESHOLD {
            break;
        }
        lk = next_lk;
    }
    let end = std::time::Instant::now();
    let duration = (end - start).as_millis();
    trace!("CORRECT\tTIME\t{}\t{}\t{}\t{}", id, k, duration, count);
    // trace!("CORRECT\tModel\t{}\n{}", id, model);
    let predictions: Vec<_> = contexts
        .iter()
        .zip(weights.into_iter())
        .map(|(ctx, weight)| (ctx.id, ctx.index, weight))
        .collect();
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
            .map(|cl| {
                let ws: Vec<_> = weights.iter().map(|ws| ws[cl]).collect();
                HMModel::new(contexts, &ws, &up_dims, &down_dims)
            })
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
    // #[allow(dead_code)]
    // fn q_value(&self, weights: &[f64], context: &Context, align_expt: &[&AlignInfo]) -> f64 {
    //     assert_eq!(weights.len(), align_expt.len());
    //     self.fractions
    //         .iter()
    //         .zip(self.models.iter())
    //         .zip(weights.iter())
    //         .zip(align_expt.iter())
    //         .map(|(((f, model), w), alns)| w * (f.ln() + model.q_value(context, alns)))
    //         .sum()
    // }
    // #[allow(dead_code)]
    // fn calc_q_value(
    //     &self,
    //     weights: &[Vec<f64>],
    //     contexts: &[Context],
    //     align_expt: &[Vec<AlignInfo>],
    // ) -> f64 {
    //     contexts
    //         .iter()
    //         .zip(weights.iter())
    //         .enumerate()
    //         .map(|(i, (ctx, ws))| {
    //             let alns: Vec<_> = align_expt.iter().map(|alns| &alns[i]).collect();
    //             self.q_value(ws, ctx, &alns)
    //         })
    //         .sum()
    // }
    fn update(&mut self, weights: &mut Vec<Vec<f64>>, contexts: &[Context], iteration: usize) {
        // E-step
        let align_expt = self.e_step(contexts, weights);
        // M-step
        self.fractions = sum_and_normalize(weights);
        let align_expt = align_expt.chunks(contexts.len());
        for (cl, (model, alns)) in self.models.iter_mut().zip(align_expt).enumerate() {
            let ws: Vec<_> = weights.iter().map(|ws| ws[cl]).collect();
            model.update(&ws, contexts, alns, iteration);
        }
    }
    // 2nd return value: the expectation of aligning information of each cluster k.
    fn e_step(&self, contexts: &[Context], weights: &mut Vec<Vec<f64>>) -> Vec<AlignInfo> {
        let mut alignments = Vec::with_capacity(contexts.len() * self.k + 1);
        let mut posterior = vec![0f64; contexts.len() * self.k];
        for (cl, (m, f)) in self.models.iter().zip(self.fractions.iter()).enumerate() {
            let begin = contexts.len() * cl;
            let f = f.ln();
            alignments.extend(contexts.iter().map(|ctx| m.align(ctx)));
            for (idx, aln) in alignments[begin..].iter().enumerate() {
                posterior[self.k * idx + cl] = aln.lk + f;
            }
        }
        for (ws, post) in weights.iter_mut().zip(posterior.chunks(self.k)) {
            assert!(post.iter().all(|x| !x.is_nan()), "{}", self);
            let total = logsumexp(post);
            ws.iter_mut()
                .zip(post)
                .for_each(|(w, &p)| *w = (p - total).exp());
        }
        alignments
    }
}

fn vec2str(xs: &[f64]) -> String {
    let xs: Vec<_> = xs
        .iter()
        .map(|&x| {
            if x < NEG_LARGE + 100f64 {
                "  !   ".to_string()
            } else {
                format!("{:6.1}", x)
            }
        })
        .collect();
    xs.join(",")
}

struct AlignInfo {
    // Likelihood. Sum of center, upstream and downstream.
    lk: f64,
    upstream: AlnToArm,
    downstream: AlnToArm,
}

struct AlnToArm {
    // reference index -> query index -> prob of match at that position,
    // reference index -> last -> prob of deletion at that position.
    // For each reference index, the sum would be one.
    match_expt: Vec<Vec<f64>>,
    // reference index + 1 -> the probability to move to the dead position at that location.
    // at the 0-th position it has the probability to drop.
    drop_expt: Vec<f64>,
    // expectation arriving at the i-th state.
    arrive_expt: Vec<f64>,
}

impl AlnToArm {
    fn merge_drop_prob<'a, I: std::iter::Iterator<Item = &'a Self>>(
        alns: I,
        weights: &[f64],
    ) -> Vec<f64> {
        let (mut total, mut sum) = (vec![], vec![]);
        for (aln, w) in alns.zip(weights.iter()) {
            if sum.is_empty() {
                sum = aln
                    .arrive_expt
                    .iter()
                    .zip(aln.drop_expt.iter())
                    .map(|(x, y)| w * x * y)
                    .collect();
                total = aln.arrive_expt.iter().map(|x| w * x).collect();
            } else {
                for (((sum, total), arrive), drop) in sum
                    .iter_mut()
                    .zip(total.iter_mut())
                    .zip(aln.arrive_expt.iter())
                    .zip(aln.drop_expt.iter())
                {
                    *sum += w * arrive * drop;
                    *total += w * arrive;
                }
            }
        }
        sum.iter_mut()
            .zip(total)
            .for_each(|(x, t)| *x = (*x + SMALL_VALUE) / (t + SMALL_VALUE));
        sum
    }
}

// Sum xs and normalize the result so that it summing up to 1.
fn sum_and_normalize(xss: &[Vec<f64>]) -> Vec<f64> {
    let mut sumed = vec![SMALL_VALUE; xss[0].len()];
    let mut total = SMALL_VALUE;
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

struct HMModel {
    #[allow(dead_code)]
    upstream_len: usize,
    downstream_len: usize,
    center: Dirichlet,
    // Dropping/continue probability at first position. Logged.
    up_drop: (f64, f64),
    upstream: Vec<DirichletMixture>,
    // Dropping/continuing probability at first position. Logged.
    down_drop: (f64, f64),
    downstream: Vec<DirichletMixture>,
}

impl std::fmt::Display for HMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for up in self.upstream.iter().rev() {
            writeln!(f, "{}", up)?;
        }
        writeln!(f, "{:.3}", self.up_drop.0.exp())?;
        writeln!(f, "Center\t{}", self.center)?;
        writeln!(f, "{:.3}", self.down_drop.0.exp())?;
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

const NEG_LARGE: f64 = -10000000000000f64;
type DP = Vec<Vec<f64>>;
impl HMModel {
    fn new(contexts: &[Context], weights: &[f64], up: &[usize], down: &[usize]) -> Self {
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
        let center = {
            let total: f64 = weights.iter().sum();
            let mut xs = weights
                .iter()
                .zip(contexts.iter().map(|ctx| &ctx.center))
                .fold(Vec::new(), |mut acc, (w, center)| {
                    if acc.is_empty() {
                        acc = center.iter().map(|x| x * w).collect();
                    } else {
                        acc.iter_mut()
                            .zip(center.iter())
                            .for_each(|(a, c)| *a += w * c);
                    }
                    acc
                });
            xs.iter_mut().for_each(|x| *x /= total);
            Dirichlet::fit(&xs, DIR_NORM)
        };
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
        let mut model = Self {
            upstream_len,
            downstream_len,
            center,
            up_drop: (0.05f64.ln(), 0.95f64.ln()),
            upstream,
            down_drop: (0.05f64.ln(), 0.95f64.ln()),
            downstream,
        };
        // Tune some-loop
        let mut lk = std::f64::NEG_INFINITY;
        for _t in 1.. {
            let alignments: Vec<_> = contexts.iter().map(|ctx| model.align(ctx)).collect();
            let new_lk = alignments.iter().map(|x| x.lk).sum::<f64>();
            model.update(weights, contexts, &alignments, 0);
            if new_lk < lk + THRESHOLD {
                break;
            }
            lk = new_lk;
        }
        model
    }
    fn up_drop(&self) -> f64 {
        self.up_drop.0
    }
    fn up_cont(&self) -> f64 {
        self.up_drop.1
    }
    fn down_drop(&self) -> f64 {
        self.down_drop.0
    }
    fn down_cont(&self) -> f64 {
        self.down_drop.1
    }
    fn lk(&self, ctx: &Context) -> f64 {
        let up_lk = {
            let forward = Self::forward(&self.upstream, &ctx.upstream, self.up_cont());
            let lk = self
                .upstream
                .iter()
                .zip(forward.iter().filter_map(|obs_q| obs_q.last()).skip(1))
                .map(|(dir, obs)| obs + dir.drop_ln());
            let lk = logsumexp_str(lk);
            logsumexp2(lk, forward[0].last().unwrap() + self.up_drop())
        };
        let down_lk = {
            let forward = Self::forward(&self.downstream, &ctx.downstream, self.down_cont());
            let lk = self
                .downstream
                .iter()
                .zip(forward.iter().filter_map(|obs_q| obs_q.last()).skip(1))
                .map(|(dir, obs)| obs + dir.drop_ln());
            let lk = logsumexp_str(lk);
            logsumexp2(lk, self.down_drop() + forward[0].last().unwrap())
        };
        let center = self.center.lk(&ctx.center);
        assert!(!(up_lk + center + down_lk).is_nan(), "\n{}\n{}", self, ctx);
        up_lk + center + down_lk
    }
    fn align(&self, ctx: &Context) -> AlignInfo {
        let (up_aln, up_lk) = {
            let forward = Self::forward(&self.upstream, &ctx.upstream, self.up_cont());
            let lk = self
                .upstream
                .iter()
                .zip(forward.iter().filter_map(|obs_q| obs_q.last()).skip(1))
                .map(|(dir, obs)| obs + dir.drop_ln())
                .chain(std::iter::once(forward[0].last().unwrap() + self.up_drop()));
            let lk = logsumexp_str(lk);
            let backward = Self::backward(&self.upstream, &ctx.upstream, self.up_cont());
            let drop = self.up_drop;
            let aln = Self::to_align(&self.upstream, &ctx.upstream, forward, &backward, drop);
            (aln, lk)
        };
        let (down_aln, down_lk) = {
            let forward = Self::forward(&self.downstream, &ctx.downstream, self.down_cont());
            let lk = self
                .downstream
                .iter()
                .zip(forward.iter().filter_map(|obs_q| obs_q.last()).skip(1))
                .map(|(dir, obs)| obs + dir.drop_ln())
                .chain(std::iter::once(
                    forward[0].last().unwrap() + self.down_drop(),
                ));
            let lk = logsumexp_str(lk);
            let backward = Self::backward(&self.downstream, &ctx.downstream, self.down_cont());
            let drop = self.down_drop;
            let aln = Self::to_align(&self.downstream, &ctx.downstream, forward, &backward, drop);
            (aln, lk)
        };
        let center = self.center.lk(&ctx.center);
        AlignInfo {
            lk: up_lk + center + down_lk,
            upstream: up_aln,
            downstream: down_aln,
        }
    }
    fn forward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)], cont: f64) -> DP {
        let mut dp = vec![vec![0f64; query.len() + 1]; refr.len() + 1];
        for j in 1..query.len() + 1 {
            dp[0][j] = NEG_LARGE;
        }
        for i in 1..refr.len() + 1 {
            let cont_prev = if 1 < i { refr[i - 2].cont_ln() } else { cont };
            dp[i][0] = dp[i - 1][0] + cont_prev + refr[i - 1].del();
        }
        for (i, dir) in refr.iter().enumerate().map(|(i, x)| (i + 1, x)) {
            let cont_prev = if 1 < i { refr[i - 2].cont_ln() } else { cont };
            for (j, obs) in query.iter().enumerate().map(|(j, x)| (j + 1, x)) {
                let del_trans = dp[i - 1][j] + dir.del();
                let match_trans = dp[i - 1][j - 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans) + cont_prev;
            }
        }
        dp
    }
    // i->j-> the probability to see j:N from the the i-1 th node.
    fn backward(refr: &[DirichletMixture], query: &[(usize, Vec<f64>)], cont: f64) -> DP {
        let mut dp = vec![vec![0f64; query.len() + 1]; refr.len() + 1];
        for j in 0..query.len() + 1 {
            dp[refr.len()][j] = NEG_LARGE;
        }
        for i in 0..refr.len() + 1 {
            dp[i][query.len()] = if 0 < i { refr[i - 1].cont_ln() } else { cont };
        }
        for (i, dir) in refr.iter().enumerate().rev() {
            let cont_prev = if 0 < i { refr[i - 1].cont_ln() } else { cont };
            for (j, obs) in query.iter().enumerate().rev() {
                let del_trans = dp[i + 1][j] + dir.del();
                let match_trans = dp[i + 1][j + 1] + dir.mat(obs);
                dp[i][j] = logsumexp2(del_trans, match_trans) + cont_prev;
            }
        }
        dp
    }
    fn to_align(
        refr: &[DirichletMixture],
        query: &[(usize, Vec<f64>)],
        mut forward: DP,
        backward: &DP,
        (init_drop, init_cont): (f64, f64),
    ) -> AlnToArm {
        // the probability to drop at the i-th position is
        // forward[i+1][0:N], then drop by drop[i].
        // the probability not to drop at the i-th position is
        // forward[i+1][0:j], suceed by (1-drop[i]),
        // then
        // 1. see the j-th obserbation at the i+1-th node and rest by backward[i+2][j+1]
        // 2. see the deltion at the i+1-th node and res by backward[i+2][j]
        // summing up all j.
        let mut drop_probs = vec![];
        if !refr.is_empty() {
            // Immediate drop.
            let drop = forward[0].last().unwrap() + init_drop;
            let suc_match = forward[0]
                .iter()
                .zip(backward[1].iter().skip(1))
                .zip(query.iter())
                .map(|((f, b), obs)| f + init_cont + refr[0].mat(obs) + b);
            let suc_del = forward[0]
                .iter()
                .zip(backward[1].iter())
                .map(|(f, b)| f + init_cont + refr[0].del() + b);
            let total = logsumexp_str(suc_match.chain(suc_del).chain(std::iter::once(drop)));
            drop_probs.push((drop - total).exp());
            for i in 0..refr.len().max(1) - 1 {
                let drop = forward[i + 1].last().unwrap() + refr[i].drop_ln();
                let suc_match = forward[i + 1]
                    .iter()
                    .zip(backward[i + 2].iter().skip(1))
                    .zip(query.iter())
                    .map(|((f, b), obs)| f + refr[i].cont_ln() + refr[i + 1].mat(obs) + b);
                let suc_del = forward[i + 1]
                    .iter()
                    .zip(backward[i + 2].iter())
                    .map(|(f, b)| f + refr[i].cont_ln() + refr[i + 1].del() + b);
                let total = logsumexp_str(suc_match.chain(suc_del).chain(std::iter::once(drop)));
                drop_probs.push((drop - total).exp());
            }
        }
        // The drop prob of the last state is almost 1 always.
        drop_probs.push(1f64 - SMALL_VALUE);
        // Arriving probability.
        // i -> probability to arrive the i-1 th node. The first element is always 1,
        // as we start with that node.
        // This is sum_j(F[i][j] * B[i][j]), which is proportional to the
        // prob to arrive at the i-1 th node.
        let step_at_prob: Vec<f64> = forward
            .iter()
            .zip(backward.iter())
            .map(|(fs, bs)| logsumexp_str(fs.iter().zip(bs.iter()).map(|(f, b)| (f + b))))
            .collect();
        // This is sum_j(F[i][N] * B[i][N]), which is proportional to the
        // prob. to *not* seeing the i th node. So, the skip(1) is needed.
        let drop_at_prob: Vec<f64> = forward
            .iter()
            .zip(backward.iter())
            .map(|(fs, bs)| fs.last().unwrap() + bs.last().unwrap())
            .collect();
        let mut arrive_prob = vec![1f64];
        arrive_prob.extend(
            step_at_prob
                .iter()
                .skip(1)
                .zip(drop_at_prob.iter())
                .map(|(&step, &drop)| (step - logsumexp2(step, drop)).exp()),
        );
        // The matching probability.
        // The probability to match i the position to the j-th obaservation is
        // forward[i][j], then suceed by refr[i-1].cont_ln(),
        // match by refr[i].match(obs[j]), then see the rest by backward[i+1][j+1].
        // The del probs is forward[i][j], then suceed by refr[i-1].cont_ln(),
        // del at refr[i], then see the rest by backward[i+1][j],
        // summing over all j.
        // Also, it shoule be multiplied by the probability to enter into the i-th state.
        // it is forward[i][j], then see the rest by backward[i][j]. summing up all over i.
        let suc_probs = std::iter::once(init_cont).chain(refr.iter().map(|dir| dir.cont_ln()));
        forward.pop();
        forward
            .iter_mut()
            .zip(refr.iter())
            .zip(backward.iter().skip(1))
            .zip(suc_probs)
            .for_each(|(((forward, dir), backward), cont_ln)| {
                let del_prob = forward
                    .iter()
                    .zip(backward.iter())
                    .map(|(f, b)| f + cont_ln + dir.del() + b);
                let del_prob = logsumexp_str(del_prob);
                forward
                    .iter_mut()
                    .zip(backward.iter().skip(1))
                    .zip(query.iter())
                    .for_each(|((f, b), obs)| *f += cont_ln + dir.mat(obs) + b);
                *forward.last_mut().unwrap() = del_prob;
                if forward.iter().any(|x| x.is_nan()) {
                    eprintln!("{:?}", forward);
                    eprintln!("{}", dir);
                    eprintln!("{}", vec2str(forward));
                    eprintln!("{}", vec2str(backward));
                    eprintln!("{}", cont_ln);
                    panic!();
                }

                let total = logsumexp(forward);
                for x in forward.iter_mut() {
                    *x = (*x - total).exp();
                }
            });
        AlnToArm {
            match_expt: forward,
            drop_expt: drop_probs,
            arrive_expt: arrive_prob,
        }
    }
    fn update<A: std::borrow::Borrow<AlignInfo>, C: std::borrow::Borrow<Context>>(
        &mut self,
        weights: &[f64],
        contexts: &[C],
        align_expt: &[A],
        _iteration: usize,
    ) {
        let contexts = contexts.iter().map(|ctx| ctx.borrow());
        let align_expt = align_expt.iter().map(|aln| aln.borrow());
        let up_alns = align_expt.clone().map(|aln| &aln.upstream);
        let down_alns = align_expt.clone().map(|aln| &aln.downstream);
        let up_drop_prob = AlnToArm::merge_drop_prob(up_alns.clone(), weights);
        let down_drop_prob = AlnToArm::merge_drop_prob(down_alns.clone(), weights);
        {
            // Update the center.
            let (mut total, mut sums) = (SMALL_VALUE, Vec::new());
            for (w, center) in weights.iter().zip(contexts.clone().map(|ctx| &ctx.center)) {
                total += w;
                if sums.is_empty() {
                    sums.extend(center.iter().map(|x| w * x));
                } else {
                    sums.iter_mut()
                        .zip(center.iter())
                        .for_each(|(a, x)| *a += x * w);
                };
            }
            sums.iter_mut().for_each(|x| *x /= total);
            self.center.update(&sums, DIR_NORM);
            let up_drop = up_drop_prob[0].max(SMALL_VALUE);
            self.up_drop = (up_drop.ln(), (1f64 - up_drop).max(SMALL_VALUE).ln());
            let down_drop = down_drop_prob[0].max(SMALL_VALUE);
            self.down_drop = (down_drop.ln(), (1f64 - down_drop).max(SMALL_VALUE).ln());
        }
        if !self.upstream.is_empty() {
            // Update the model in the upstream region.
            let contexts = contexts.clone().map(|ctx| ctx.upstream.as_slice());
            self.upstream
                .iter_mut()
                .zip(up_drop_prob.iter().skip(1))
                .for_each(|(dir, &drop)| dir.set_drop_prob(drop));
            Self::update_oneside(&mut self.upstream, weights, contexts, up_alns);
        }
        if !self.downstream.is_empty() {
            let contexts = contexts.clone().map(|ctx| ctx.downstream.as_slice());
            self.downstream
                .iter_mut()
                .zip(down_drop_prob.iter().skip(1))
                .for_each(|(dir, &drop)| dir.set_drop_prob(drop));
            Self::update_oneside(&mut self.downstream, weights, contexts, down_alns);
        }
    }
    fn update_oneside<
        'a,
        'b,
        I: std::iter::Iterator<Item = &'a AlnToArm>,
        J: std::iter::Iterator<Item = &'b [(usize, Vec<f64>)]>,
    >(
        parameters: &mut [DirichletMixture],
        weights: &[f64],
        contexts: J,
        alignments: I,
    ) {
        // Allocate all the elements at once, it is inefficient?
        // j -> (total weight of the j, del prob at j, c -> Dirichlet)
        type SummaryOnLocation = (f64, f64, Vec<(f64, Vec<f64>)>);
        let mut sum_ups: Vec<SummaryOnLocation> = parameters
            .iter()
            .map(|dir| {
                let state: Vec<_> = dir
                    .dirichlets
                    .iter()
                    .map(|(_, dir)| (SMALL_VALUE, vec![0f64; dir.dim()]))
                    .collect();
                (SMALL_VALUE, SMALL_VALUE, state)
            })
            .collect();
        // Alns: j -> i -> prob to align j to the i-th.
        for ((&w_data, ctx), aln) in weights.iter().zip(contexts).zip(alignments) {
            let (matches, arrive) = (aln.match_expt.iter(), aln.arrive_expt.iter().skip(1));
            for (((weight_sum, del_sum, state), aln), arrive) in
                sum_ups.iter_mut().zip(matches).zip(arrive)
            {
                assert!(*aln.last().unwrap() <= 1f64);
                *del_sum += w_data * arrive * aln.last().unwrap();
                *weight_sum += w_data * arrive;
                for (&(unit, ref obs), aln_prob) in ctx.iter().zip(aln.iter()) {
                    state[unit].0 += w_data * arrive * aln_prob;
                    state[unit]
                        .1
                        .iter_mut()
                        .zip(obs.iter())
                        .for_each(|(s, x)| *s += w_data * arrive * aln_prob * x);
                }
            }
        }
        parameters
            .iter_mut()
            .zip(sum_ups)
            .for_each(|(dir, (weight_sum, del_sum, state))| {
                dir.set_del_prob(del_sum / weight_sum);
                let match_sum: f64 = state.iter().map(|x| x.0).sum();
                dir.dirichlets.iter_mut().zip(state).for_each(
                    |((fr, di), (mat_prob, mut total))| {
                        *fr = (mat_prob / match_sum).ln();
                        if !APPROX || LOG_WEIGHT_FILTER < *fr {
                            total.iter_mut().for_each(|x| *x /= mat_prob);
                            di.update(&total, DIR_NORM);
                        }
                    },
                );
            });
    }
}

struct DirichletMixture {
    // **log dropping probability**
    drop_prob: f64,
    cont_prob: f64,
    // **log deletion prob**
    del_prob: f64,
    mat_prob: f64,
    // **(ln cluster composition, dirichlet) **.
    dirichlets: Vec<(f64, Dirichlet)>,
}

impl std::fmt::Display for DirichletMixture {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let dump: Vec<_> = self
            .dirichlets
            .iter()
            .enumerate()
            .filter(|&(_, &(fr, _))| fr.exp() > 0.1)
            .map(|(unit, (fr, dir))| format!("{}:{:.3}:{}", unit, fr.exp(), dir))
            .collect();
        write!(
            f,
            "{:.2}\t{}\t{:.2}",
            self.del_prob.exp(),
            dump.join("\t"),
            self.drop_prob.exp()
        )
    }
}

impl DirichletMixture {
    fn drop_ln(&self) -> f64 {
        self.drop_prob
    }
    fn cont_ln(&self) -> f64 {
        self.cont_prob
    }
    // Return log del_prob
    fn del(&self) -> f64 {
        self.del_prob
    }
    // Return log (1-e)Dir(prob|param).
    fn mat(&self, &(unit, ref prob): &(usize, Vec<f64>)) -> f64 {
        self.dirichlets
            .get(unit)
            .map(|(frac, dir)| {
                if APPROX && frac + self.mat_prob < LOG_WEIGHT_FILTER {
                    self.mat_prob + frac
                } else {
                    self.mat_prob + frac + dir.lk(prob)
                }
            })
            .unwrap()
    }
    fn set_del_prob(&mut self, del_prob: f64) {
        self.del_prob = del_prob.max(SMALL_VALUE).ln();
        self.mat_prob = (1f64 - del_prob).max(SMALL_VALUE).ln();
    }
    fn set_drop_prob(&mut self, drop_prob: f64) {
        self.drop_prob = drop_prob.max(SMALL_VALUE).ln();
        self.cont_prob = (1f64 - drop_prob).max(SMALL_VALUE).ln();
    }
    // Create new instance from the given data.
    fn new(observed: &[(usize, f64, &[f64])], dims: &[usize]) -> Self {
        let del_prob = 0.1;
        let drop_prob = 0.1;
        let mut post_probs = vec![vec![]; dims.len()];
        let mut prob_weights = vec![vec![]; dims.len()];
        let mut fractions = vec![SMALL_VALUE; dims.len()];
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
                let fr = (fr / sum).ln();
                let total: f64 = ws.iter().sum::<f64>() + SMALL_VALUE;
                let mut xs = vec![0f64; dim];
                for (w, prob) in ws.iter().zip(probs.iter()) {
                    xs.iter_mut()
                        .zip(prob.iter())
                        .for_each(|(a, x)| *a += w * x);
                }
                xs.iter_mut().for_each(|x| *x /= total);
                (fr, Dirichlet::fit(&xs, DIR_NORM))
            })
            .collect();
        Self {
            mat_prob: (1f64 - del_prob).ln(),
            del_prob: del_prob.ln(),
            cont_prob: (1f64 - drop_prob).ln(),
            drop_prob: drop_prob.ln(),
            dirichlets,
        }
    }
}

// TODO:Maybe we can make this faster?
// TODO: Maybe we should return -Inf for the empty iterator.
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
        0 => 0f64,
        1 => max,
        _ => accum.ln() + max,
    }
}

// LogSumExp(xs). Return if there's nan or inf.
// If the vector is empty, it returns zero.
fn logsumexp(xs: &[f64]) -> f64 {
    match xs.len() {
        0 => 0f64,
        1 => xs[0],
        2 => logsumexp2(xs[0], xs[1]),
        _ => {
            let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
            let sum = xs
                .iter()
                .filter(|&x| !(APPROX && 20f64 < max - x))
                .map(|x| (x - max).exp())
                .sum::<f64>()
                .ln();
            max + sum
        }
    }
}

// Log(exp(x) + exp(y)) = x + Log(1+exp(y-x)).
fn logsumexp2(x: f64, y: f64) -> f64 {
    let (x, y) = (x.max(y), x.min(y));
    if APPROX && 20f64 < x - y {
        x
    } else {
        x + (1f64 + (y - x).exp()).ln()
    }
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
    fn easy1() {
        let seed = 423304982094;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let hap1s = (vec![], 0, vec![]);
        let hap2s = (vec![], 1, vec![]);
        let cl_num = { (vec![], 2, vec![]) };
        let is_uninformative = (vec![], false, vec![]);
        let num = 20;
        let reads = gen_reads(&hap1s, &hap2s, &cl_num, &is_uninformative, &mut rng, num);
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
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
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
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
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
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
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
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
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
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
        let (result, _) = clustering_inner(&reads, 2, 20, &mut rng);
        for (_, _, res) in result.iter() {
            println!("{}", vec2str(res));
        }
        let answer = vec![vec![0; num], vec![1; num]].concat();
        let pred: Vec<_> = result.iter().map(|(_, _, x)| (x[0] < x[1]) as u8).collect();
        let idx = crate::local_clustering::rand_index(&answer, &pred);
        assert!(idx > 0.8);
    }
}
