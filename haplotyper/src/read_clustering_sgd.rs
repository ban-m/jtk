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
    _cov: f64,
    config: &ReadClusteringConfig,
) -> (Vec<(u64, usize, Vec<f64>)>, f64, usize) {
    let seed = (unit + 1) * cluster_num as u64;
    let (contexts, _up_units, _down_units) = convert_to_contexts(&reads, unit, config);
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let take_len = match cluster_num {
        1 => 1,
        _ => config.repeat_num,
    };
    // let (up_len, down_len) = (up_units.len(), down_units.len());
    let (asn, likelihood) = (0..take_len)
        .map(|_| simple_clustering_inner(&contexts, cluster_num, config, &mut rng))
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    (asn, likelihood, cluster_num)
}

// Convert reads into contexts.
// return (contexts, # of units in upstream, and # of units in downstream.
fn convert_to_contexts<'a>(
    reads: &[&'a EncodedRead],
    unit: u64,
    config: &ReadClusteringConfig,
) -> (Vec<Context<'a>>, HashMap<u64, usize>, HashMap<u64, usize>) {
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
            let (up_slot, down_slot) = match read.nodes[i].is_forward {
                true => (&mut upstream_units, &mut downstream_units),
                false => (&mut downstream_units, &mut upstream_units),
            };
            for n in read.nodes.iter().take(i) {
                let len = up_slot.len();
                if unit_counts.contains_key(&n.unit) {
                    up_slot.entry(n.unit).or_insert(len);
                }
            }
            for n in read.nodes.iter().skip(i + 1) {
                let len = down_slot.len();
                if unit_counts.contains_key(&n.unit) {
                    down_slot.entry(n.unit).or_insert(len);
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

#[derive(Debug, Clone)]
struct Context<'a> {
    // The ID of the original encoded read.
    id: u64,
    // The original index of this context.
    index: usize,
    #[allow(dead_code)]
    unit: u64,
    cluster: &'a [f64],
    // Only upstream/downstream information is preversed, the
    // ordering is arbitrary.
    upstream: Vec<(usize, &'a [f64])>,
    downstream: Vec<(usize, &'a [f64])>,
}

impl std::fmt::Display for Context<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let upstream: Vec<_> = self
            .upstream
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .rev()
            .collect();
        let downstream: Vec<_> = self
            .downstream
            .iter()
            .map(|(u, c)| format!("({},{:?})", u, c))
            .collect();
        write!(
            f,
            "{}\t({:?})\t{}",
            upstream.join("-"),
            self.cluster,
            downstream.join("\t")
        )
    }
}

impl<'a> Context<'a> {
    fn new(
        read: &'a EncodedRead,
        index: usize,
        up_map: &HashMap<u64, usize>,
        down_map: &HashMap<u64, usize>,
    ) -> Self {
        let (unit, cluster) = (
            read.nodes[index].unit,
            read.nodes[index].posterior.as_slice(),
        );
        let nodes = read.nodes.iter();
        let (upstream, downstream) = match read.nodes[index].is_forward {
            true => {
                let upstream: Vec<_> = nodes
                    .clone()
                    .take(index)
                    .rev()
                    .filter_map(|n| up_map.get(&n.unit).map(|&i| (i, n.posterior.as_slice())))
                    .collect();
                let downstream: Vec<_> = nodes
                    .clone()
                    .skip(index + 1)
                    .filter_map(|n| down_map.get(&n.unit).map(|&i| (i, n.posterior.as_slice())))
                    .collect();
                (upstream, downstream)
            }
            false => {
                let upstream: Vec<_> = nodes
                    .clone()
                    .skip(index + 1)
                    .filter_map(|n| up_map.get(&n.unit).map(|&i| (i, n.posterior.as_slice())))
                    .collect();
                let downstream: Vec<_> = nodes
                    .clone()
                    .take(index)
                    .rev()
                    .filter_map(|n| down_map.get(&n.unit).map(|&i| (i, n.posterior.as_slice())))
                    .collect();
                (upstream, downstream)
            }
        };
        let id = read.id;
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

fn simple_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    _config: &ReadClusteringConfig,
    rng: &mut R,
) -> (Vec<(u64, usize, Vec<f64>)>, f64) {
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    let mut model = ReadClusteringModel::new(contexts.len());
    let mut parameters = model.new_params(&contexts, k, rng);
    let mut loss = model.loss(&parameters, &contexts);
    let mut no_improve = 0;
    // MAP estimation.
    for t in 0.. {
        for i in rand::seq::index::sample(rng, contexts.len(), contexts.len()) {
            let data = &contexts[i];
            model.update(&mut parameters, data, i);
        }
        model.update_learn_rate();
        let new_loss = model.loss(&parameters, &contexts);
        use println as debug;
        debug!("ReadClustering\tLoss\t{}\t{}\t{}", id, t, new_loss);
        if new_loss < loss {
            no_improve = 0;
            loss = new_loss;
        } else {
            no_improve += 1;
        }
        if 100 < no_improve {
            break;
        }
    }
    (model.predict(&parameters, contexts), loss)
}

// Maybe it is more appropriate call this struct as ``optimizer'' and
// integrate the learning rate or other parameters into this struct. How about it?
#[derive(Debug, Clone)]
pub struct ReadClusteringModel {
    learning_rate: f64,
    reg_term: f64,
    data_size: usize,
    clock: usize,
}

impl ReadClusteringModel {
    fn new(data_size: usize) -> Self {
        Self {
            reg_term: 0f64,
            data_size,
            learning_rate: 0.01,
            clock: 0,
        }
    }
    fn loss(&self, param: &ReadClusteringParam, contexts: &[Context]) -> f64 {
        let fractions: Vec<_> = param
            .consensus
            .iter()
            .map(|cons| cons.to_fraction())
            .collect();
        let data_loss: f64 = contexts
            .iter()
            .enumerate()
            .map(|(i, ctx)| self.loss_data(param, ctx, i, &fractions))
            .sum();
        let reg_loss: f64 = self.reg_loss(&param) * self.reg_term;
        data_loss + reg_loss
    }
    fn reg_loss(&self, param: &ReadClusteringParam) -> f64 {
        param.consensus.iter().map(|cons| cons.norm()).sum()
    }
    fn loss_data(
        &self,
        param: &ReadClusteringParam,
        context: &Context,
        i: usize,
        fractions: &[Fraction],
    ) -> f64 {
        assert_eq!(param.weights[i].len(), param.consensus.len());
        let total = logsumexp(&param.weights[i]);
        param.weights[i]
            .iter()
            .zip(param.consensus.iter())
            .zip(fractions.iter())
            .map(|((w, cons), fs)| (w - total).exp() * cons.loss(fs, context))
            .sum()
    }
    fn update_learn_rate(&mut self) {
        self.learning_rate *= (self.clock as f64) / (self.clock + 1) as f64;
        self.clock += 1;
    }
    // TODO: Numerical -> hand-written differential
    fn update(&self, param: &mut ReadClusteringParam, data: &Context, i: usize) {
        for (k, cons) in param.consensus.iter().enumerate() {
            println!("{}\t{}", k, cons);
        }
        println!();
        // Calc gradient.
        let mut grad_weight = self.get_weight_grad(param, data, i);
        grad_weight
            .iter_mut()
            .for_each(|x| *x *= self.learning_rate);
        let mut grad_cons = self.get_consensus_grad(param, data, i);
        grad_cons.iter_mut().for_each(|cons| {
            cons.add_scalar(2f64 * self.reg_term / self.data_size as f64);
            cons.mul_scalar(self.learning_rate);
        });
        // Update.
        for (cons, grad) in param.consensus.iter_mut().zip(grad_cons.iter()) {
            cons.sub(grad);
        }
        param.weights[i]
            .iter_mut()
            .zip(grad_weight.iter())
            .for_each(|(x, y)| *x -= y);
    }
    fn get_weight_grad(
        &self,
        param: &mut ReadClusteringParam,
        data: &Context,
        i: usize,
    ) -> Vec<f64> {
        let fractions: Vec<_> = param
            .consensus
            .iter()
            .map(|cons| cons.to_fraction())
            .collect();
        let total = logsumexp(&param.weights[i]);
        let weights: Vec<_> = param.weights[i].iter().map(|x| (x - total).exp()).collect();
        let loss: Vec<_> = param
            .consensus
            .iter()
            .zip(fractions.iter())
            .map(|(cons, fs)| cons.loss(fs, data))
            .collect();
        let mean_loss: f64 = weights.iter().zip(loss.iter()).map(|(x, y)| x * y).sum();
        weights
            .iter()
            .zip(loss.iter())
            .map(|(z_k, loss_k)| z_k * (loss_k - mean_loss))
            .collect()
    }
    // Gradient of consensus.
    fn get_consensus_grad(
        &self,
        param: &mut ReadClusteringParam,
        data: &Context,
        i: usize,
    ) -> Vec<Consensus> {
        let total = logsumexp(&param.weights[i]);
        let ws: Vec<_> = param.weights[i].clone();
        param
            .consensus
            .iter_mut()
            .zip(ws)
            .map(|(cons, w)| {
                let z_ik = (w - total).exp();
                assert!(0f64 <= z_ik && z_ik <= 1f64);
                let mut grad = cons.get_grad(data);
                println!("Data\t{:.3}\t{}", z_ik, data);
                println!("Cons\t{}", cons);
                println!("Grad\t{}\n", grad);
                grad.mul_scalar(z_ik);
                grad
            })
            .collect()
    }
    fn predict(
        &self,
        param: &ReadClusteringParam,
        contexts: &[Context],
    ) -> Vec<(u64, usize, Vec<f64>)> {
        contexts
            .iter()
            .enumerate()
            .map(|(i, ctx)| {
                let total = logsumexp(&param.weights[i]);
                let posterior: Vec<_> =
                    param.weights[i].iter().map(|x| (x - total).exp()).collect();
                (ctx.id, ctx.index, posterior)
            })
            .collect()
    }
    // Generate new parameters.
    fn new_params<R: Rng>(
        &self,
        contexts: &[Context],
        k: usize,
        rng: &mut R,
    ) -> ReadClusteringParam {
        let mut num_cluster_in_up_units = vec![];
        let mut num_cluster_in_down_units = vec![];
        let num_cluster = contexts[0].cluster.len();
        for ctx in contexts.iter() {
            for &(unit, prob) in ctx.upstream.iter() {
                let len = prob.len();
                let ext_len = (unit + 1).saturating_sub(num_cluster_in_up_units.len());
                num_cluster_in_up_units.extend(vec![0; ext_len]);
                num_cluster_in_up_units[unit] = len;
            }
            for &(unit, prob) in ctx.downstream.iter() {
                let len = prob.len();
                let ext_len = (unit + 1).saturating_sub(num_cluster_in_down_units.len());
                num_cluster_in_down_units.extend(vec![0; ext_len]);
                num_cluster_in_down_units[unit] = len;
            }
        }
        ReadClusteringParam::new(
            k,
            contexts.len(),
            num_cluster,
            &num_cluster_in_up_units,
            &num_cluster_in_down_units,
            rng,
        )
    }
}

#[derive(Debug, Clone)]
pub struct ReadClusteringParam {
    // The weight of the i-th read on the k-th cluster.
    weights: Vec<Vec<f64>>,
    // k->j->m
    consensus: Vec<Consensus>,
}

impl ReadClusteringParam {
    fn new<R: Rng>(
        k: usize,
        num_data: usize,
        num_cluster: usize,
        in_up_units: &[usize],
        in_down_units: &[usize],
        rng: &mut R,
    ) -> Self {
        let weights: Vec<_> = (0..num_data).map(|_| Consensus::gen_prob(k, rng)).collect();
        let consensus: Vec<_> = (0..k)
            .map(|_| Consensus::new(num_cluster, in_up_units, in_down_units, rng))
            .collect();
        Self { weights, consensus }
    }
}

// It can be regarded as "log-ed" version of the posterior probabilites.
#[derive(Debug, Clone)]
pub struct Consensus {
    upstream: Vec<Vec<f64>>,
    downstream: Vec<Vec<f64>>,
    center: Vec<f64>,
}

impl std::fmt::Display for Consensus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, probs) in self.upstream.iter().enumerate().rev() {
            let probs: Vec<_> = probs.iter().map(|x| format!("{:.2}", x)).collect();
            write!(f, "({},[{}])", i, probs.join(","))?;
        }
        let probs: Vec<_> = self.center.iter().map(|x| format!("{:.2}", x)).collect();
        write!(f, "\t([{}])\t", probs.join(","))?;
        for (i, probs) in self.downstream.iter().enumerate() {
            let probs: Vec<_> = probs.iter().map(|x| format!("{:.2}", x)).collect();
            write!(f, "({},[{}])", i, probs.join(","))?;
        }
        Ok(())
    }
}

impl Consensus {
    fn with(upstream: Vec<Vec<f64>>, center: Vec<f64>, downstream: Vec<Vec<f64>>) -> Self {
        Self {
            upstream,
            center,
            downstream,
        }
    }
    fn to_fraction(&self) -> Fraction {
        let mut center = logsumexp(&self.center);
        let mut upstream: Vec<_> = self.upstream.iter().map(|logs| logsumexp(logs)).collect();
        let mut downstream: Vec<_> = self.downstream.iter().map(|logs| logsumexp(logs)).collect();
        let all_units: Vec<_> = std::iter::once(center)
            .chain(upstream.iter().copied())
            .chain(downstream.iter().copied())
            .collect();
        let total = logsumexp(&all_units);
        center = (center - total).exp();
        upstream.iter_mut().for_each(|x| *x = (*x - total).exp());
        downstream.iter_mut().for_each(|x| *x = (*x - total).exp());
        Fraction {
            center,
            upstream,
            downstream,
        }
    }
    fn new<R: Rng>(
        num_cluster: usize,
        in_up_units: &[usize],
        in_down_units: &[usize],
        rng: &mut R,
    ) -> Self {
        // Return *log probabilities*.
        let center = Self::gen_prob(num_cluster, rng);
        let upstream: Vec<_> = in_up_units
            .iter()
            .map(|&dim| Self::gen_prob(dim, rng))
            .collect();
        let downstream: Vec<_> = in_down_units
            .iter()
            .map(|&dim| Self::gen_prob(dim, rng))
            .collect();
        Self::with(upstream, center, downstream)
    }
    // Log-prob of dim-dimentional probability.
    fn gen_prob<R: Rng>(dim: usize, rng: &mut R) -> Vec<f64> {
        let mut slots = vec![0f64; dim];
        for _ in 0..500 {
            slots[rng.gen_range(0..dim)] += 1f64;
        }
        let sum: f64 = slots.iter().sum();
        slots.iter_mut().for_each(|x| *x = (*x / sum).ln());
        slots
    }
    fn norm(&self) -> f64 {
        std::iter::once(&self.center)
            .chain(self.upstream.iter())
            .chain(self.downstream.iter())
            .map(|probs| probs.iter().map(|x| x * x).sum::<f64>())
            .sum()
    }
    fn mul_scalar(&mut self, factor: f64) {
        let center = self.center.iter_mut();
        let upstream = self.upstream.iter_mut().flat_map(|xs| xs.iter_mut());
        let downstream = self.downstream.iter_mut().flat_map(|xs| xs.iter_mut());
        center
            .chain(upstream)
            .chain(downstream)
            .for_each(|x| *x *= factor);
    }
    fn loss(&self, fs: &Fraction, context: &Context) -> f64 {
        let center = fs.center * sq_diff(&self.center, &context.cluster);
        let upstream: f64 = context
            .upstream
            .iter()
            .map(|&(unit, ref prob)| fs.upstream[unit] * sq_diff(&self.upstream[unit], prob))
            .sum();
        let downstream: f64 = context
            .downstream
            .iter()
            .map(|&(unit, ref prob)| fs.downstream[unit] * sq_diff(&self.downstream[unit], prob))
            .sum();
        center + upstream + downstream
    }
    // TODO: This is numerical diff.( f(x+h) - f(x-h) ) / 2h
    fn get_grad(&mut self, data: &Context) -> Self {
        let diff = 0.00001;
        let center: Vec<_> = (0..self.center.len())
            .map(|i| {
                let loss_b = {
                    self.center[i] -= diff;
                    let frac = self.to_fraction();
                    self.loss(&frac, data)
                };
                let loss_f = {
                    self.center[i] += 2f64 * diff;
                    let frac = self.to_fraction();
                    self.loss(&frac, data)
                };
                println!("C\t{}\t{}-{}", i, loss_f, loss_b);
                self.center[i] -= diff;
                (loss_f - loss_b) / 2f64 / diff
            })
            .collect();
        let upstream: Vec<Vec<_>> = (0..self.upstream.len())
            .map(|j| {
                (0..self.upstream[j].len())
                    .map(|m| {
                        let loss_b = {
                            self.upstream[j][m] -= diff;
                            let frac = self.to_fraction();
                            self.loss(&frac, data)
                        };
                        let loss_f = {
                            self.upstream[j][m] += 2f64 * diff;
                            let frac = self.to_fraction();
                            self.loss(&frac, data)
                        };
                        self.upstream[j][m] += diff;
                        (loss_f - loss_b) / 2f64 / diff
                    })
                    .collect()
            })
            .collect();
        let downstream: Vec<_> = (0..self.downstream.len())
            .map(|j| {
                (0..self.downstream[j].len())
                    .map(|m| {
                        let loss_b = {
                            self.downstream[j][m] -= diff;
                            let frac = self.to_fraction();
                            self.loss(&frac, data)
                        };
                        let loss_f = {
                            self.downstream[j][m] += 2f64 * diff;
                            let frac = self.to_fraction();
                            self.loss(&frac, data)
                        };
                        self.downstream[j][m] += diff;
                        (loss_f - loss_b) / 2f64 / diff
                    })
                    .collect()
            })
            .collect();
        Self::with(upstream, center, downstream)
    }
    fn add_scalar(&mut self, grad: f64) {
        let center = self.center.iter_mut();
        let upstream = self.upstream.iter_mut().flat_map(|xs| xs.iter_mut());
        let downstream = self.downstream.iter_mut().flat_map(|xs| xs.iter_mut());
        center
            .chain(upstream)
            .chain(downstream)
            .for_each(|x| *x += grad);
    }
    fn sub(&mut self, other: &Self) {
        self.center
            .iter_mut()
            .zip(other.center.iter())
            .for_each(|(x, y)| *x -= y);
        for (xs, ys) in self.upstream.iter_mut().zip(other.upstream.iter()) {
            xs.iter_mut().zip(ys).for_each(|(x, y)| *x -= y);
        }
        for (xs, ys) in self.downstream.iter_mut().zip(other.downstream.iter()) {
            xs.iter_mut().zip(ys).for_each(|(x, y)| *x -= y);
        }
    }
}

#[derive(Debug, Clone)]
pub struct Fraction {
    upstream: Vec<f64>,
    downstream: Vec<f64>,
    center: f64,
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

fn sq_diff(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(xs.len(), ys.len());
    xs.iter().zip(ys.iter()).map(|(x, y)| (x - y).powi(2)).sum()
}

#[cfg(test)]
pub mod tests {
    impl<'a> Context<'a> {
        fn with_attrs(
            id: u64,
            index: usize,
            unit: u64,
            cluster: &'a [f64],
            upstream: Vec<(usize, &'a [f64])>,
            downstream: Vec<(usize, &'a [f64])>,
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
                let forward: Vec<(usize, Vec<f64>)> = vec![];
                let (backward, cluster) = if i < len / 2 {
                    (vec![(0, vec![0f64])], vec![0f64, -1000f64])
                } else {
                    (vec![(1, vec![0f64])], vec![-1000f64, 0f64])
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
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(94);
        let config = ReadClusteringConfig::default();
        let (asn, _likelihood) = simple_clustering_inner(&contexts, 2, &config, &mut rng);
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
    fn test_sq_diff_test() {
        let test = sq_diff(&[0.1, 0.9], &[0.5, 0.5]);
        let answer = 0.4 * 0.4 * 2f64;
        assert!((test - answer).abs() < 0.0001);
    }
    #[test]
    fn test_cons_sub() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(5945);
        let mut cons = Consensus::new(2, &[2, 2, 2], &[3, 3, 3], &mut rng);
        let cons_2 = cons.clone();
        cons.sub(&cons_2);
        assert!(cons.norm() < 0.0001);
    }
    #[test]
    fn test_cons_add_scalar() {
        let center = vec![0f64, 0f64];
        let upstream = vec![vec![0f64; 3], vec![0f64; 3], vec![0f64; 3]];
        let downstream = vec![vec![0f64; 2]; 3];
        let mut cons = Consensus::with(upstream, center, downstream);
        cons.add_scalar(1f64);
        let answer = (2 + 3 * 3 + 2 * 3) as f64;
        assert!((answer - cons.norm()).abs() < 0.0001);
    }
    #[test]
    fn test_cons_loss() {
        let cluster = vec![0.9f64.ln(), 0.1f64.ln()];
        let upstream: Vec<_> = vec![(0, vec![0.1f64.ln(), 0.9f64.ln()])];
        let downstream: Vec<_> = vec![(0, vec![0.1f64.ln(), 0.9f64.ln()])];
        let upstream: Vec<_> = upstream.iter().map(|(x, y)| (*x, y.as_slice())).collect();
        let cons = {
            let upstream: Vec<_> = upstream.iter().map(|x| x.1.to_vec()).collect();
            let downstream: Vec<_> = downstream.iter().map(|x| x.1.to_vec()).collect();
            let center = cluster.clone();
            Consensus::with(upstream, center, downstream)
        };
        let downstream: Vec<_> = downstream.iter().map(|(x, y)| (*x, y.as_slice())).collect();
        let context = Context::with_attrs(0, 0, 0, &cluster, upstream, downstream);
        let fraction = cons.to_fraction();
        let loss = cons.loss(&fraction, &context);
        let answer = 0f64;
        assert!((loss - answer).abs() < 0.0001);
    }
}
