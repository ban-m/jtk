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

use rayon::prelude::*;

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

fn correct(ds: &mut DataSet) {
    let key_values: Vec<_> = ds
        .selected_chunks
        .iter()
        .map(|ref_unit| (ref_unit.id, ref_unit.cluster_num))
        // .filter(|x| x.0 == 99)
        .collect();
    let result: Vec<_> = key_values
        .par_iter()
        .map(|&(id, cluster_num)| correct_unit(ds, id, cluster_num))
        .collect();
    let result: HashMap<u64, Vec<(usize, u64)>> =
        result.iter().fold(HashMap::new(), |mut acc, results| {
            for &(id, pos, cluster) in results {
                acc.entry(id).or_default().push((pos, cluster));
            }
            acc
        });
    for read in ds.encoded_reads.iter_mut() {
        if let Some(corrected) = result.get(&read.id) {
            for &(pos, cluster) in corrected {
                read.nodes[pos].cluster = cluster;
            }
        }
    }
}

fn correct_unit(ds: &DataSet, unit_id: u64, k: usize) -> Vec<(u64, usize, u64)> {
    let repeat_num = 20;
    let coverage_thr = 5;
    let reads: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
        .collect();
    if reads.is_empty() {
        return vec![];
    }
    let (new_clustering, dist, _new_k) = (1..=k)
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
    trace!("CORRECT\tPickedDist\t{}", -dist);
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

fn clustering(reads: &[&EncodedRead], config: &Config) -> (Vec<(u64, usize, u64)>, f64, usize) {
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
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(config.seed);
    let cluster_num = config.cluster_num;
    let (asn, dist, _) = simple_clustering_inner(&contexts, cluster_num, &mut rng);
    (asn, -dist, cluster_num)
}

// KL(p||q) = sum(p_i ln p/q).
fn kl_divergence(ps: &[f64], qs: &[f64]) -> f64 {
    assert_eq!(ps.len(), qs.len());
    ps.iter()
        .zip(qs.iter())
        .filter(|&(&p, _)| 0.0001 < p)
        .map(|(p, q)| p * (p / q.max(0.00001)).ln())
        .sum()
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
) -> (Vec<(u64, usize, u64)>, f64, f64) {
    // ID of this trial.
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    let mut assignments: Vec<_> = contexts.iter().map(|_| rng.gen_range(0..k)).collect();
    let mut model = SimpleModel::new(contexts, &assignments, k);
    let mut dist: f64 = contexts.iter().map(|ctx| model.get_dist(ctx)).sum();
    trace!("CORRECT\tDist\t{}\t{}", id, dist);
    for _ in 0..20 {
        model.update(&mut assignments, contexts);
        model = SimpleModel::new(contexts, &assignments, k);
        let next_dist = contexts.iter().map(|ctx| model.get_dist(ctx)).sum();
        trace!("CORRECT\tDist\t{}\t{}", id, next_dist);
        if (dist - next_dist) < 0.001 {
            break;
        } else {
            dist = next_dist;
        }
    }
    trace!("CORRECT\tModel\t{}\n{}", id, model);
    let predictions: Vec<_> = contexts
        .iter()
        .zip(assignments.iter())
        .map(|(ctx, &cluster)| (ctx.id, ctx.index, cluster as u64))
        .collect();
    let num_param = model.num_parameters() as f64;
    trace!("CORRECT\tFinal\t{}\t{}\t{:.4}\t{}", id, k, dist, num_param);
    (predictions, dist, num_param)
}

#[derive(Debug, Clone)]
struct SimpleModel {
    unit: u64,
    // Fraction of each consensus. Currently not used.
    fraction: Vec<f64>,
    // Consensus of each cluster. Actually, it is a bug of probability distribution on chunks.
    consensus: Vec<Consensus>,
}

impl std::fmt::Display for SimpleModel {
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

impl SimpleModel {
    // Create new model.
    fn new(contexts: &[Context<'_>], assignments: &[usize], cluster_num: usize) -> Self {
        let unit = contexts[0].unit;
        let fraction = {
            let mut fraction = vec![0f64; cluster_num];
            for &asn in assignments {
                fraction[asn] += 1f64;
            }
            let sum: f64 = fraction.iter().sum();
            fraction.iter_mut().for_each(|x| *x /= sum);
            fraction
        };
        let consensus: Vec<_> = (0..cluster_num)
            .map(|cl| {
                let contexts: Vec<_> = contexts
                    .iter()
                    .zip(assignments)
                    .filter_map(|(c, &asn)| (asn == cl).then(|| c))
                    .collect();
                Consensus::new(&contexts)
            })
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
    fn update(&mut self, assignments: &mut [usize], contexts: &[Context]) {
        for (asn, context) in assignments.iter_mut().zip(contexts.iter()) {
            *asn = self.get_nearest_consensus(context).0;
        }
    }
    fn get_nearest_consensus(&self, context: &Context<'_>) -> (usize, f64) {
        self.get_dists(context)
            .into_iter()
            .enumerate()
            .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
            .unwrap()
    }
    fn get_dists(&self, context: &Context) -> Vec<f64> {
        self.consensus.iter().map(|m| m.get_dist(context)).collect()
    }
    fn get_dist(&self, context: &Context) -> f64 {
        self.get_nearest_consensus(context).1
    }
}

fn vec2str(xs: &[f64]) -> String {
    let dump: Vec<_> = xs.iter().map(|x| format!("{:.3}", x)).collect();
    dump.join(":")
}

#[derive(Debug, Clone)]
struct Consensus {
    center: Vec<f64>,
    forward: HashMap<u64, (f64, Vec<f64>)>,
    backward: HashMap<u64, (f64, Vec<f64>)>,
}

impl std::fmt::Display for Consensus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut backward: Vec<_> = self.backward.iter().collect();
        backward.sort_by_key(|x| x.0);
        let backward = backward
            .iter()
            .map(|(u, c)| format!("{}-{}", u, vec2str(&c.1)));
        let mut forward: Vec<_> = self.forward.iter().collect();
        forward.sort_by_key(|x| x.0);
        let forward = forward
            .iter()
            .map(|(u, c)| format!("{}-{}", u, vec2str(&c.1)));
        let center = std::iter::once(format!("C-{}", vec2str(&self.center)));
        let dump: Vec<_> = backward.chain(center).chain(forward).collect();
        write!(f, "{}", dump.join("\t"))
    }
}

impl Consensus {
    fn new(xs: &[&Context<'_>]) -> Self {
        // Column sum, take the max value.
        let mut forward: HashMap<u64, Vec<f64>> = HashMap::new();
        let mut backward: HashMap<u64, Vec<f64>> = HashMap::new();
        let mut center = xs.iter().fold(Vec::new(), |mut sum, ctx| {
            if sum.is_empty() {
                sum = vec![0f64; ctx.cluster.len()];
            }
            sum.iter_mut()
                .zip(ctx.cluster.iter())
                .for_each(|(x, y)| *x += y);
            sum
        });
        for ctx in xs.iter() {
            for &(unit, post) in ctx.forward.iter() {
                forward
                    .entry(unit)
                    .and_modify(|sum| {
                        sum.iter_mut().zip(post.iter()).for_each(|(x, y)| *x += y);
                    })
                    .or_insert_with(|| post.to_vec());
            }
            for &(unit, post) in ctx.backward.iter() {
                backward
                    .entry(unit)
                    .and_modify(|sum| {
                        sum.iter_mut().zip(post.iter()).for_each(|(x, y)| *x += y);
                    })
                    .or_insert_with(|| post.to_vec());
            }
        }
        // Normalizing.
        let forward_sum = forward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let backward_sum = forward.values().flat_map(|x| x).sum::<f64>().max(1f64);
        let center_sum: f64 = center.iter().sum();
        center.iter_mut().for_each(|x| *x /= center_sum);
        let forward: HashMap<u64, (f64, Vec<f64>)> = forward
            .into_iter()
            .map(|(unit, mut prob)| {
                let sum: f64 = prob.iter().sum();
                prob.iter_mut().for_each(|x| *x /= sum);
                (unit, (sum / forward_sum, prob))
            })
            .collect();
        let backward: HashMap<u64, (f64, Vec<f64>)> = backward
            .into_iter()
            .map(|(unit, mut prob)| {
                let sum: f64 = prob.iter().sum();
                prob.iter_mut().for_each(|x| *x /= sum);
                (unit, (sum / backward_sum, prob))
            })
            .collect();
        Self {
            center,
            forward,
            backward,
        }
    }
    fn get_dist(&self, context: &Context) -> f64 {
        let kl_center = match self.center.is_empty() {
            false => kl_divergence(&self.center, context.cluster),
            true => -(0.00001f64).ln(),
        };
        Self::kl_divergences(&self.forward, &context.forward)
            + Self::kl_divergences(&self.backward, &context.backward)
            + kl_center
    }
    fn kl_divergences(cons: &HashMap<u64, (f64, Vec<f64>)>, query: &[(u64, &[f64])]) -> f64 {
        query
            .iter()
            .map(|(unit, distribution)| match cons.get(unit) {
                Some((f, model)) => kl_divergence(model, distribution) - f.ln(),
                None => {
                    let cluster_num = distribution.len() as f64;
                    let entropy: f64 = distribution.iter().map(|x| x.max(0.00001).ln()).sum();
                    -entropy / cluster_num - cluster_num.ln() - (0.0001f64).ln()
                }
            })
            .sum()
    }
}
