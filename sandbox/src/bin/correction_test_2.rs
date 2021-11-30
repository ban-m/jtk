const PRIOR_FRACTION: f64 = 1f64;
const CLUSTER_PRIOR: f64 = 0f64;
use definitions::*;
use std::collections::HashMap;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    correct(&mut ds);
    println!("{}", serde_json::ser::to_string_pretty(&ds).unwrap());
    Ok(())
}

#[derive(Debug, Clone)]
pub struct Config {
    repeat_num: usize,
    seed: u64,
    cluster_num: usize,
    coverage_thr: usize,
    pub to_use_offset: bool,
    focal: u64,
}

impl Config {
    pub fn repeat_num(&self) -> usize {
        self.repeat_num
    }
    pub fn new(
        repeat_num: usize,
        seed: u64,
        cluster_num: usize,
        focal: u64,
        to_use_offset: bool,
        coverage: usize,
    ) -> Self {
        Self {
            repeat_num,
            seed,
            cluster_num,
            focal,
            to_use_offset,
            coverage_thr: coverage,
        }
    }
}

// use rayon::prelude::*;

fn correct(ds: &mut DataSet) {
    let posterior_distributions: Vec<_> = ds
        .selected_chunks
        // .par_iter()
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .filter(|&(unit, _)| unit == 99)
        .map(|(id, cluster_num)| correct_unit(ds, id, cluster_num))
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
    for read in ds.encoded_reads.iter_mut() {
        if let Some(corrected) = result.remove(&read.id) {
            for (pos, post) in corrected {
                read.nodes[pos].cluster = argmax(&post) as u64;
                read.nodes[pos].posterior = post;
            }
        }
    }
}

fn correct_unit(ds: &DataSet, unit_id: u64, k: usize) -> Vec<(u64, usize, Vec<f64>)> {
    // let repeat_num = 20;
    let repeat_num = 1;
    let coverage_thr = 5;
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
        .collect();
    if reads.is_empty() {
        return vec![];
    }
    let (new_clustering, lk, _new_k) = (1..=k)
        .flat_map(|k| match k {
            1 => std::iter::repeat(k).take(1),
            _ => std::iter::repeat(k).take(repeat_num),
        })
        .enumerate()
        .map(|(i, k)| {
            let seed = unit_id * (i * k) as u64;
            let config = Config::new(repeat_num, seed, k, unit_id, true, coverage_thr);
            clustering(&reads, &config)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    trace!("CORRECT\tPickedLK\t{}", lk);
    // let id_is_hapa: HashMap<_, _> = ds
    //     .raw_reads
    //     .iter()
    //     .map(|r| (r.id, r.name.contains("hapA")))
    //     .collect();
    // for (id, _, cl) in new_clustering.iter() {
    //     trace!("CORRECT\tPredict\t{}\t{}\t{}", id, id_is_hapa[id], cl);
    // }
    new_clustering
}

use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

fn clustering(
    reads: &[&EncodedRead],
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
                if read.nodes[index].unit == config.focal {
                    buffer.push(Context::new(read, index, &unit_counts));
                }
            }
        }
        buffer
    };
    let cluster_num = config.cluster_num;
    trace!("CORRECT\t{}\t{}", cluster_num, contexts.len());
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(config.seed);
    let (asn, likelihood, _) = simple_clustering_inner(&contexts, cluster_num, &mut rng);
    (asn, likelihood, cluster_num)
}

#[derive(Debug, Clone)]
pub struct Context<'a> {
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
    pub fn dump_with(&self, fdigit: usize, bdigit: usize) -> String {
        let mut f_slots = vec!["---".to_string(); fdigit];
        let mut b_slots = vec!["---".to_string(); bdigit];
        for (i, (u, c)) in self.forward.iter().enumerate().take(fdigit) {
            f_slots[i] = format!("{}-{:?}", u, c);
        }
        for (i, (u, c)) in self.backward.iter().enumerate().take(bdigit) {
            b_slots[bdigit - i - 1] = format!("{}-{:?}", u, c);
        }
        format!(
            "{}\t{}-{:?}\t{}",
            b_slots.join("\t"),
            self.unit,
            self.cluster,
            f_slots.join("\t")
        )
    }
    pub fn with_attrs(
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

use log::*;

pub fn simple_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
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
    let mut lk: f64 = contexts.iter().map(|ctx| model.get_likelihood(ctx)).sum();
    trace!("CORRECT\tLikelihood\t{}\t{}", id, lk);
    trace!("CORRECT\tModel\t{}\n{}", id, model);
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
        let dump: Vec<_> = self
            .fraction
            .iter()
            .zip(self.consensus.iter())
            .map(|(frac, m)| format!("{:.3}-({})", frac, m))
            .collect();
        write!(f, "{}", dump.join("\n"))
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
            self.update_weight(weight, context);
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
    dump.join(":")
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
        let backward = backward
            .iter()
            .map(|(u, c)| format!("{}-{:.2}[{}]", u, c.0, vec2str(&c.1)));
        let mut forward: Vec<_> = self.forward.iter().collect();
        forward.sort_by_key(|x| x.0);
        let forward = forward
            .iter()
            .map(|(u, c)| format!("{}-{:.2}[{}]", u, c.0, vec2str(&c.1)));
        let center = std::iter::once(format!("C-1.0[{}]", vec2str(&self.center)));
        let dump: Vec<_> = backward.chain(center).chain(forward).collect();
        write!(f, "{}", dump.join("\t"))
    }
}

impl Consensus {
    fn new(xs: &[Context<'_>], weights: &[Vec<f64>], cl: usize) -> Self {
        // Column sum, take the max value.
        let center = {
            let mut center = vec![CLUSTER_PRIOR; xs[0].cluster.len()];
            for ctx in xs.iter() {
                center
                    .iter_mut()
                    .zip(ctx.cluster.iter())
                    .for_each(|(x, y)| *x += y);
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
        let forward_sum = forward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let forward: HashMap<u64, (f64, Vec<f64>)> = forward
            .into_iter()
            .map(|(unit, prob)| {
                let sum: f64 = prob.iter().sum();
                (unit, (sum / forward_sum, prob))
            })
            .collect();
        let backward_sum = backward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let backward: HashMap<u64, (f64, Vec<f64>)> = backward
            .into_iter()
            .map(|(unit, prob)| {
                let sum: f64 = prob.iter().sum();
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
                    .map(|(frac, param)| frac.ln() + dirichlet(prob, param))
                    .unwrap_or((0.000000001f64).ln())
            })
            .sum();
        assert!(!forward.is_nan());
        let backward: f64 = context
            .backward
            .iter()
            .map(|(unit, prob)| {
                self.backward
                    .get(unit)
                    .map(|(frac, param)| frac.ln() + dirichlet(prob, param))
                    .unwrap_or((0.0000001f64).ln())
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
        .map(|(p, q)| (q - 1f64) * p.ln())
        .sum();
    let sum: f64 = param.iter().sum();
    let lk = unsafe { lgamma(sum) - param.iter().map(|&q| lgamma(q)).sum::<f64>() + without_scale };
    assert!(!lk.is_nan(), "{:?}\t{:?}", prob, param);
    lk
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
