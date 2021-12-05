//! This module contains codes to correct (or polish) local clusterings.
//! Usually, it is sufficient to call `ds.correct_clustering(&config)` with an appropriate configuration `config`.

// Hyper parameters. Maybe we need to tune this, but
// `tuning` is something we should postpone as far as possible.
const PRIOR_FRACTION: f64 = 1f64;
const CLUSTER_PRIOR: f64 = 1f64;

const SMALL: f64 = 0.0000000000000000000001;
const LOGSMALL: f64 = -100f64;
use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::HashMap;

pub trait DirichletCorrection {
    /// Correct the posterior probability on each node in each read.
    fn correct_clustering(&mut self, config: &Config);
}

/// Top-level configuration. It handles the parameters shared by the dirichlet correction procedure.
#[derive(Debug, Clone)]
pub struct Config {
    repeat_num: usize,
    coverage_thr: usize,
}

impl Config {
    pub fn new(repeat_num: usize, coverage_thr: usize) -> Self {
        Self {
            repeat_num,
            coverage_thr,
        }
    }
}

impl std::default::Default for Config {
    fn default() -> Self {
        Self {
            repeat_num: 20,
            coverage_thr: 4,
        }
    }
}

impl DirichletCorrection for DataSet {
    fn correct_clustering(&mut self, config: &Config) {
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
                .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
                .unwrap()
                .0
        };
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = result.remove(&read.id) {
                for (pos, post) in corrected {
                    read.nodes[pos].cluster = argmax(&post) as u64;
                    read.nodes[pos].posterior = post;
                }
            }
        }
    }
}

pub fn correct_unit(
    ds: &DataSet,
    unit_id: u64,
    k: usize,
    config: &Config,
) -> Vec<(u64, usize, Vec<f64>)> {
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
        .collect();
    if reads.is_empty() {
        return vec![];
    }
    let (mut new_clustering, lk, new_k) = (1..=k)
        .flat_map(|k| match k {
            1 => std::iter::repeat(k).take(1),
            _ => std::iter::repeat(k).take(config.repeat_num),
        })
        .enumerate()
        .map(|(i, k)| {
            let seed = unit_id * (i * k) as u64;
            clustering(&reads, unit_id, seed, k, &config)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    trace!("CORRECT\tPickedLK\t{}", lk);
    let pad_len = k.saturating_sub(new_k);
    for (_, _, prob) in new_clustering.iter_mut() {
        prob.extend(std::iter::repeat(0f64).take(pad_len));
    }
    new_clustering
}

use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

fn clustering(
    reads: &[&EncodedRead],
    unit: u64,
    seed: u64,
    cluster_num: usize,
    config: &Config,
) -> (Vec<(u64, usize, Vec<f64>)>, f64, usize) {
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    unit_counts.retain(|_, count| config.coverage_thr < *count);
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
    trace!("CORRECT\t{}\t{}", cluster_num, contexts.len());
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let (asn, likelihood, _) = simple_clustering_inner(&contexts, cluster_num, config, &mut rng);
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
    // Only forward/backward information is preversed, the
    // ordering is arbitrary.
    forward: Vec<(u64, &'a [f64])>,
    backward: Vec<(u64, &'a [f64])>,
}

impl std::fmt::Display for Context<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let forward: Vec<_> = self
            .forward
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .collect();
        let backward: Vec<_> = self
            .backward
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .collect();
        write!(
            f,
            "{}\t({},{:?})\t{}",
            forward.join("-"),
            self.unit,
            self.cluster,
            backward.join("-")
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
        let forward: Vec<_> = nodes
            .clone()
            .skip(index + 1)
            .map(|n| (n.unit, n.posterior.as_slice()))
            .filter(|n| unit_counts.contains_key(&n.0))
            .collect();
        let backward: Vec<_> = nodes
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
                forward,
                backward,
            },
            false => Self {
                id,
                index,
                unit,
                cluster,
                forward: backward,
                backward: forward,
            },
        }
    }
}

// Return vector of (read id, position, posterior probability),
// the final likelihood, and the number of parameters.
fn simple_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    _config: &Config,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64, f64) {
    // ID of this trial.
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    let gen_weight = |_| -> Vec<f64> {
        let mut bufs = vec![0f64; k];
        for _ in 0..100 {
            bufs[rng.gen_range(0..k)] += 1f64;
        }
        let sum: f64 = bufs.iter().sum();
        bufs.iter_mut().for_each(|x| *x /= sum);
        let sum: f64 = bufs.iter().sum();
        assert!((1f64 - sum).abs() < 0.00001);
        bufs
    };
    let mut weights: Vec<_> = contexts.iter().map(gen_weight).collect();
    let mut model = DirichletModel::new(contexts, &weights, k);
    trace!("CORRECT\tModel\t{}\n{}", id, model);
    let mut lk: f64 = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
    trace!("CORRECT\tLikelihood\t{}\t{}", id, lk);
    for _ in 0..20 {
        trace!("CORRECT\tModel\t{}\n{}", id, model);
        model.update(&mut weights, contexts);
        model = DirichletModel::new(contexts, &weights, k);
        let next_lk = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
        trace!("CORRECT\tLikelihood\t{}\t{}", id, next_lk);
        if (next_lk - lk) < 0.001 {
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
    let num_param = model.num_parameters() as f64;
    trace!("CORRECT\tFinal\t{}\t{}\t{:.4}\t{}", id, k, lk, num_param);
    (predictions, lk, num_param)
}

#[derive(Debug, Clone)]
struct DirichletModel {
    unit: u64,
    fraction: Vec<f64>,
    // Consensus of each cluster. Actually, it is a bug of dirichlet distributions on chunks.
    consensus: Vec<Consensus>,
}

impl std::fmt::Display for DirichletModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (i, (fr, cons)) in self.fraction.iter().zip(self.consensus.iter()).enumerate() {
            writeln!(f, "{}\t{:.3}\n{}", i, fr, cons)?;
        }
        Ok(())
    }
}

impl DirichletModel {
    // Create new model.
    fn new(contexts: &[Context<'_>], weights: &[Vec<f64>], cluster_num: usize) -> Self {
        assert!(!weights.is_empty());
        let unit = contexts[0].unit;
        let fraction = {
            let mut sum_of_weights = vec![PRIOR_FRACTION; cluster_num];
            for ws in weights.iter() {
                sum_of_weights.iter_mut().zip(ws).for_each(|(x, y)| *x += y);
            }
            let sum: f64 = sum_of_weights.iter().sum();
            sum_of_weights.iter_mut().for_each(|x| *x /= sum);
            sum_of_weights
        };
        let consensus: Vec<_> = (0..cluster_num)
            .map(|cl| Consensus::new(&contexts, weights, cl))
            .collect();
        Self {
            unit,
            fraction,
            consensus,
        }
    }
    fn num_parameters(&self) -> usize {
        use std::collections::HashSet;
        let mut used_units: HashSet<_> = HashSet::new();
        used_units.insert(self.unit);
        for cons in self.consensus.iter() {
            used_units.extend(cons.forward.keys().copied());
            used_units.extend(cons.backward.keys().copied());
        }
        let model_param = used_units.len();
        model_param + self.fraction.len() - 1
    }
    fn update(&mut self, weights: &mut [Vec<f64>], contexts: &[Context]) {
        for (weight, context) in weights.iter_mut().zip(contexts.iter()) {
            // let prev = weight.clone();
            self.update_weight(weight, context);
            // trace!("[{}]->[{}]", vec2str(&prev), vec2str(weight));
        }
    }
    fn update_weight(&self, weight: &mut [f64], context: &Context<'_>) {
        // Unnormalized lielihood.
        weight
            .iter_mut()
            .zip(&self.fraction)
            .zip(&self.consensus)
            .for_each(|((w, f), c)| *w = f.ln() + c.lk(context));
        let total = logsumexp(weight);
        weight.iter_mut().for_each(|w| *w = (*w - total).exp());
    }
    fn get_likelihood(&self, context: &Context) -> f64 {
        let lks: Vec<_> = self
            .fraction
            .iter()
            .zip(&self.consensus)
            .map(|(f, c)| f.ln() + c.lk(context))
            .collect();
        logsumexp(&lks)
    }
}

fn vec2str(xs: &[f64]) -> String {
    let dump: Vec<_> = xs.iter().map(|x| format!("{:.1}", x)).collect();
    dump.join(",")
}

#[derive(Debug, Clone)]
struct Consensus {
    // Dirichlet parameters.
    center: Vec<f64>,
    // The first element is normalized across sum, the second element is not.
    forward: HashMap<u64, (f64, Vec<f64>)>,
    backward: HashMap<u64, (f64, Vec<f64>)>,
}

impl std::fmt::Display for Consensus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut backward: Vec<_> = self.backward.iter().collect();
        backward.sort_by_key(|x| x.0);
        for (u, c) in backward.iter() {
            writeln!(f, "B\t{}\t{:.2}\t[{}]", u, c.0, vec2str(&c.1))?;
        }
        writeln!(f, "C\t1.0\t{}", vec2str(&self.center))?;
        let mut forward: Vec<_> = self.forward.iter().collect();
        forward.sort_by_key(|x| x.0);
        for (u, c) in forward.iter() {
            writeln!(f, "F\t{}\t{:.2}\t[{}]", u, c.0, vec2str(&c.1))?;
        }
        Ok(())
    }
}

impl Consensus {
    fn new(xs: &[Context<'_>], weights: &[Vec<f64>], cl: usize) -> Self {
        // Column sum, take the max value.
        let center = {
            let mut center = vec![CLUSTER_PRIOR; xs[0].cluster.len()];
            for (ctx, w) in xs.iter().zip(weights.iter()).map(|(x, ws)| (x, ws[cl])) {
                center
                    .iter_mut()
                    .zip(ctx.cluster.iter())
                    .for_each(|(x, y)| *x += w * y);
            }
            let sum: f64 = center.iter().sum();
            center.iter_mut().for_each(|x| *x /= sum);
            center
        };
        let mut forward: HashMap<u64, Vec<f64>> = HashMap::new();
        let mut backward: HashMap<u64, Vec<f64>> = HashMap::new();
        for (ctx, w) in xs.iter().zip(weights.iter()).map(|(ctx, ws)| (ctx, ws[cl])) {
            for &(unit, post) in ctx.forward.iter() {
                let insert = || vec![CLUSTER_PRIOR; post.len()];
                let post = post.iter().map(|x| x * w);
                let dirichlet = forward.entry(unit).or_insert_with(insert);
                dirichlet.iter_mut().zip(post).for_each(|(x, y)| *x += y)
            }
            for &(unit, post) in ctx.backward.iter() {
                let insert = || vec![CLUSTER_PRIOR; post.len()];
                let post = post.iter().map(|x| x * w);
                let dirichlet = backward.entry(unit).or_insert_with(insert);
                dirichlet.iter_mut().zip(post).for_each(|(x, y)| *x += y);
            }
        }
        // Normalizing.
        fn normalize(xs: &mut [f64]) -> f64 {
            let sum: f64 = xs.iter().sum();
            xs.iter_mut().for_each(|x| *x /= sum);
            sum
        }
        let forward_sum = forward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let forward: HashMap<u64, (f64, Vec<f64>)> = forward
            .into_iter()
            .map(|(unit, mut prob)| {
                let sum = normalize(&mut prob);
                (unit, (sum / forward_sum, prob))
            })
            .collect();
        let backward_sum = backward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let backward: HashMap<u64, (f64, Vec<f64>)> = backward
            .into_iter()
            .map(|(unit, mut prob)| {
                let sum = normalize(&mut prob);
                (unit, (sum / backward_sum, prob))
            })
            .collect();
        Self {
            center,
            forward,
            backward,
        }
    }
    fn lk(&self, context: &Context) -> f64 {
        let center_lk = dirichlet(context.cluster, &self.center);
        assert!(!center_lk.is_nan());
        let forward: f64 = context
            .forward
            .iter()
            .map(|(unit, prob)| {
                self.forward
                    .get(unit)
                    .map(|(frac, param)| frac.max(SMALL).ln() + dirichlet(prob, param))
                    .unwrap_or(LOGSMALL)
            })
            .sum();
        assert!(!forward.is_nan());
        let backward: f64 = context
            .backward
            .iter()
            .map(|(unit, prob)| {
                self.backward
                    .get(unit)
                    .map(|(frac, param)| frac.max(SMALL).ln() + dirichlet(prob, param))
                    .unwrap_or(LOGSMALL)
            })
            .sum();
        assert!(!backward.is_nan(), "{}\n{}\n{}", context, self, backward);
        center_lk + forward + backward
    }
}

// Log Dir(p|q) = log ( G(sum(q_i))/prod(G(q_i)) prod(p_i^{q_i-1}))
// = log G(sum(q_i)) - sum(log(G(q_i))) + sum((q_i-1)log(p_i)).
// here, G is the gamma function.
fn dirichlet(prob: &[f64], param: &[f64]) -> f64 {
    assert_eq!(prob.len(), param.len());
    let without_scale: f64 = prob
        .iter()
        .zip(param.iter())
        .map(|(p, q)| (q - 1f64) * p.max(SMALL).ln())
        .sum();
    let sum: f64 = param.iter().sum();
    let coef = unsafe {
        let denom: f64 = param.iter().map(|&q| lgamma(q)).sum::<f64>();
        lgamma(sum) - denom
    };
    coef + without_scale
}

#[link(name = "m")]
extern "C" {
    fn lgamma(x: f64) -> f64;
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

#[cfg(test)]
mod tests {
    impl<'a> Context<'a> {
        fn with_attrs(
            id: u64,
            index: usize,
            unit: u64,
            cluster: &'a [f64],
            forward: Vec<(u64, &'a [f64])>,
            backward: Vec<(u64, &'a [f64])>,
        ) -> Self {
            Self {
                id,
                index,
                unit,
                cluster,
                forward,
                backward,
            }
        }
        // fn dump_with(&self, fdigit: usize, bdigit: usize) -> String {
        //     let mut f_slots = vec!["---".to_string(); fdigit];
        //     let mut b_slots = vec!["---".to_string(); bdigit];
        //     for (i, (u, c)) in self.forward.iter().enumerate().take(fdigit) {
        //         f_slots[i] = format!("{}-{:?}", u, c);
        //     }
        //     for (i, (u, c)) in self.backward.iter().enumerate().take(bdigit) {
        //         b_slots[bdigit - i - 1] = format!("{}-{:?}", u, c);
        //     }
        //     format!(
        //         "{}\t{}-{:?}\t{}",
        //         b_slots.join("\t"),
        //         self.unit,
        //         self.cluster,
        //         f_slots.join("\t")
        //     )
        // }
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
        let config = Config::default();
        let (asn, _likelihood, _) = simple_clustering_inner(&contexts, 2, &config, &mut rng);
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
    fn dirichlet_test() {
        let pdf = dirichlet(&[0.1, 0.9], &[1f64, 100f64]);
        let answer = -5.8255;
        eprintln!("{}\n{}", pdf, answer);
        assert!((pdf - answer).abs() < 0.0001);

        let pdf = dirichlet(&[0.5, 0.5], &[1f64, 1f64]);
        let answer = 0f64;
        eprintln!("{}\n{}", pdf, answer);
        assert!((pdf - answer).abs() < 0.0001);

        let pdf = dirichlet(&[0.2, 0.6, 0.2], &[1f64, 1f64, 1f64]);
        let answer = 0.693147;
        eprintln!("{}\n{}", pdf, answer);
        assert!((pdf - answer).abs() < 0.0001);

        let pdf = dirichlet(&[0.2, 0.6, 0.2], &[2f64, 6f64, 2f64]);
        let answer = 2.2413317;
        let diff = (pdf - answer).abs();
        eprintln!("{}\n{}", pdf, answer);
        assert!(diff < 0.001, "{},{},{}", pdf, answer, diff);
    }
}
