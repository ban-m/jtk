use definitions::*;
use rayon::prelude::*;
const SEED: u64 = 1221;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::HashMap;
use std::collections::HashSet;
#[derive(Debug, Clone)]
pub struct Config {
    repeat_num: usize,
    seed: u64,
    cluster_num: usize,
    coverage_thr: usize,
    focal: u64,
}

impl Config {
    pub fn new(
        repeat_num: usize,
        seed: u64,
        cluster_num: usize,
        focal: u64,
        coverage: usize,
    ) -> Self {
        Self {
            repeat_num,
            seed,
            cluster_num,
            focal,
            coverage_thr: coverage,
        }
    }
}

pub trait ClusteringCorrection {
    fn correct_clustering(self, repeat_num: usize, coverage_thr: usize) -> Self;
    fn correct_clustering_em(self, repeat_num: usize, coverage_thr: usize) -> Self;
}

impl ClusteringCorrection for DataSet {
    fn correct_clustering(mut self, repeat_num: usize, coverage_thr: usize) -> Self {
        let id_to_name: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
        let id_to_desc: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
        let result: Vec<_> = self
            .selected_chunks
            .par_iter()
            .flat_map(|ref_unit| {
                let unit_id = ref_unit.id;
                let (read_indices, reads): (Vec<usize>, Vec<_>) = self
                    .encoded_reads
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| r.nodes.iter().any(|n| n.unit == unit_id))
                    .unzip();
                let k = ref_unit.cluster_num;
                let config = Config::new(repeat_num, SEED, k, unit_id, coverage_thr);
                if reads.is_empty() {
                    debug!("Unit {} does not appear in any read.", unit_id);
                    return vec![];
                }
                let new_clustering = clustering(&self.selected_chunks, &reads, &config);
                if log_enabled!(log::Level::Debug) {
                    for cl in 0..k {
                        for (read, &cluster) in reads.iter().zip(new_clustering.iter()) {
                            let id = read.id;
                            if cluster == cl {
                                let name = id_to_name[&id];
                                let desc = match id_to_desc.get(&id) {
                                    Some(res) => res.as_str(),
                                    None => "",
                                };
                                debug!("IMP\t{}\t{}\t{}\t{}\t{}", unit_id, cl, id, name, desc);
                            }
                        }
                    }
                }
                let mut result = vec![];
                for ((read, read_idx), cluster) in
                    reads.into_iter().zip(read_indices).zip(new_clustering)
                {
                    for (idx, node) in read.nodes.iter().enumerate() {
                        if node.unit == unit_id {
                            result.push((read_idx, read.id, idx, node.unit, cluster));
                        }
                    }
                }
                result
            })
            .collect();
        for (read_idx, read_id, position, unit_id, cluster) in result {
            assert_eq!(self.encoded_reads[read_idx].id, read_id);
            assert_eq!(self.encoded_reads[read_idx].nodes[position].unit, unit_id);
            self.encoded_reads[read_idx].nodes[position].cluster = cluster as u64;
        }
        self
    }
    fn correct_clustering_em(mut self, repeat_num: usize, coverage_thr: usize) -> Self {
        // let id_to_name: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.name)).collect();
        // let id_to_desc: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, &r.desc)).collect();
        let result: Vec<_> = self
            .selected_chunks
            .par_iter()
            .map(|ref_unit| {
                let unit_id = ref_unit.id;
                let reads: Vec<_> = self
                    .encoded_reads
                    .iter()
                    .filter(|r| r.nodes.iter().any(|n| n.unit == unit_id))
                    .collect();
                let k = ref_unit.cluster_num;
                if reads.is_empty() {
                    debug!("Unit {} does not appear in any read.", unit_id);
                    return vec![];
                }
                let (new_clustering, lk) = (0..repeat_num as u64)
                    .map(|s| {
                        let config = Config::new(repeat_num, unit_id * s, k, unit_id, coverage_thr);
                        em_clustering(&reads, &config)
                    })
                    .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
                    .unwrap();
                trace!("MaxLK:{}", lk);
                new_clustering
            })
            .collect();
        let result: HashMap<u64, Vec<(usize, u64)>> =
            result.iter().fold(HashMap::new(), |mut acc, results| {
                for &(id, pos, cluster) in results {
                    acc.entry(id).or_default().push((pos, cluster));
                }
                acc
            });
        for read in self.encoded_reads.iter_mut() {
            if let Some(corrected) = result.get(&read.id) {
                for &(pos, cluster) in corrected {
                    read.nodes[pos].cluster = cluster;
                }
            }
        }
        self
    }
}

pub fn clustering(units: &[Unit], reads: &[&EncodedRead], config: &Config) -> Vec<usize> {
    let weights: Vec<_> = reads
        .iter()
        .map(|r| r.nodes.iter().find(|n| n.unit == config.focal).unwrap())
        .map(|node| {
            let mut weight = vec![0.; config.cluster_num];
            weight[node.cluster as usize] += 1.;
            weight
        })
        .collect();
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(config.focal * SEED);
    let (asn, lk) = em_clustering_old(units, reads, weights, config);
    let (asn, lk) = (0..config.repeat_num)
        .map(|_| {
            let weights: Vec<_> = reads
                .iter()
                .map(|_| {
                    let mut weight = vec![0.; config.cluster_num];
                    weight[rng.gen_range(0..config.cluster_num)] += 1.;
                    weight
                })
                .collect();
            em_clustering_old(units, reads, weights, config)
        })
        .fold((asn, lk), |(argmax, max_lk), (arg, lk)| {
            if lk < max_lk {
                (argmax, max_lk)
            } else {
                (arg, lk)
            }
        });
    trace!("OPTLK\t{}", lk);
    asn
}

fn em_clustering_old(
    units: &[Unit],
    reads: &[&EncodedRead],
    mut weights: Vec<Vec<f64>>,
    config: &Config,
) -> (Vec<usize>, f64) {
    trace!("==================");
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    let use_units: HashSet<_> = unit_counts
        .iter()
        .filter(|&(_, &c)| c > config.coverage_thr)
        .map(|(&x, _)| x)
        .collect();
    let k = config.cluster_num;
    let mut model = Model::new(units, reads, &weights, k, &use_units);
    let mut diff = 10.;
    let mut lk = reads.iter().map(|read| model.lk(read)).sum::<f64>();
    trace!("LK:{}", lk);
    while diff > 0.00000001 {
        weights = reads.iter().map(|read| model.weight(read)).collect();
        model = Model::new(units, reads, &weights, k, &use_units);
        let new_lk = reads.iter().map(|read| model.lk(read)).sum::<f64>();
        trace!("LK:{}", lk);
        diff = new_lk - lk;
        assert!(diff > -0.01, "{}", diff);
        lk = new_lk;
    }
    for (idx, w) in weights.iter().enumerate() {
        let w: Vec<_> = w.iter().map(|x| format!("{:.1}", x)).collect();
        let lk = model.lk(&reads[idx]);
        trace!("WEIGHT\t{}\t{}\t{}", idx, w.join("\t"), lk);
    }
    let predictions: Vec<_> = weights
        .iter()
        .map(|ws| {
            ws.iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|x| x.0)
                .unwrap()
        })
        .collect();
    (predictions, lk)
}

struct Model {
    // Fraction for each cluster.
    fraction: Vec<f64>,
    // Cluster -> Unit -> Category
    // For all i in 0..cluster_num, j in units,
    // category_fraction[i][j].iter().sum::<f64>() == 1 holds.
    category_fraction: Vec<HashMap<u64, Vec<f64>>>,
}
const SMALL_WEIGHT: f64 = 0.00000001;

impl Model {
    fn new(
        units: &[Unit],
        reads: &[&EncodedRead],
        weights: &[Vec<f64>],
        k: usize,
        use_units: &HashSet<u64>,
    ) -> Self {
        let category_num: HashMap<_, _> = units
            .iter()
            .filter(|u| use_units.contains(&u.id))
            .map(|u| (u.id, u.cluster_num))
            .collect();
        assert!(weights.iter().all(|ws| ws.len() == k));
        let total_weights: Vec<_> = (0..k)
            .map(|cluster| weights.iter().map(|ws| ws[cluster]).sum::<f64>())
            .collect();
        assert!((reads.len() as f64 - total_weights.iter().sum::<f64>()).abs() < 0.001);
        let fraction: Vec<_> = total_weights
            .iter()
            .map(|f| f / weights.len() as f64 + SMALL_WEIGHT)
            .collect();
        let sums_to_one = (1. - fraction.iter().sum::<f64>()).abs() < 0.0001;
        assert!(
            sums_to_one,
            "{:?}\t{}\t{}",
            fraction,
            weights.len(),
            reads.len(),
        );
        let mut category_fraction: Vec<HashMap<_, Vec<f64>>> = (0..k)
            .map(|_| {
                category_num
                    .iter()
                    .map(|(&unit_id, &category_num)| {
                        let slots = vec![SMALL_WEIGHT; category_num + 1];
                        (unit_id, slots)
                    })
                    .collect()
            })
            .collect();
        for (read, weight) in reads.iter().zip(weights.iter()) {
            for (cluster, w) in weight.iter().enumerate() {
                for node in read.nodes.iter() {
                    let (d, m) = (node.unit, node.cluster as usize);
                    category_fraction[cluster]
                        .entry(d)
                        .and_modify(|slots| slots[m] += w);
                }
            }
        }
        for clusterwise_weight in category_fraction.iter_mut() {
            for position_fraction in clusterwise_weight.values_mut() {
                let sum = position_fraction.iter().sum::<f64>();
                position_fraction.iter_mut().for_each(|x| *x /= sum);
            }
        }
        for xs in category_fraction.iter() {
            for x in xs.values() {
                assert!((1. - x.iter().sum::<f64>()).abs() < 0.001);
            }
        }
        Self {
            category_fraction,
            fraction,
        }
    }
    fn weight(&self, read: &EncodedRead) -> Vec<f64> {
        let log_weight: Vec<_> = self
            .fraction
            .iter()
            .zip(self.category_fraction.iter())
            .map(|(f, cluster)| {
                let read_lk = read
                    .nodes
                    .iter()
                    .filter_map(|node| cluster.get(&node.unit).map(|xs| xs[node.cluster as usize]))
                    .map(|x| x.ln())
                    .sum::<f64>();
                read_lk + f.ln()
            })
            .collect();
        let log_total_weight = logsumexp(&log_weight);
        let weight: Vec<_> = log_weight
            .iter()
            .map(|x| (x - log_total_weight).exp())
            .collect();
        assert!((1. - weight.iter().sum::<f64>()).abs() < 0.001);
        weight
    }
    fn lk(&self, read: &EncodedRead) -> f64 {
        let log_lks: Vec<_> = self
            .fraction
            .iter()
            .zip(self.category_fraction.iter())
            .map(|(f, cluster)| {
                let read_lk = read
                    .nodes
                    .iter()
                    .filter_map(|node| cluster.get(&node.unit).map(|xs| xs[node.cluster as usize]))
                    .map(|x| x.ln())
                    .sum::<f64>();
                read_lk + f.ln()
            })
            .collect();
        logsumexp(&log_lks)
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    } else if xs.len() == 1 {
        xs[0]
    } else {
        let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
        let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
        assert!(sum >= 0., "{:?}->{}", xs, sum);
        max + sum
    }
}

/// Return the id of the read, the position at that read, and the clusters predicted.
pub fn em_clustering(reads: &[&EncodedRead], config: &Config) -> (Vec<(u64, usize, u64)>, f64) {
    trace!("==================");
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    let use_units: HashSet<_> = unit_counts
        .iter()
        .filter(|&(_, &c)| c > config.coverage_thr)
        .map(|(&x, _)| x)
        .collect();
    let k = config.cluster_num;
    let contexts: Vec<_> = {
        let mut buffer = vec![];
        for read in reads.iter() {
            for index in 0..read.nodes.len() {
                if read.nodes[index].unit == config.focal {
                    buffer.push(Context::new(read, index, &use_units));
                }
            }
        }
        buffer
    };
    // for (i, (context, read)) in contexts.iter().zip(reads.iter()).enumerate() {
    //     let focal = read.nodes.iter().find(|n| n.unit == config.focal).unwrap();
    //     let backward = context
    //         .backward
    //         .iter()
    //         .rev()
    //         .map(|(c, u)| format!("{}-{}", c, u));
    //     let forward = context.forward.iter().map(|(c, u)| format!("{}-{}", c, u));
    //     let center = std::iter::once(format!("{}-{}", context.unit, context.cluster));
    //     let dump: Vec<_> = backward.chain(center).chain(forward).collect();
    //     debug!("Context\t{}\t{}", i, dump.join("\t"));
    // }
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(config.seed);
    let mut weights: Vec<_> = (0..contexts.len())
        .map(|_| {
            let mut weight = vec![0f64; k];
            const NUM: usize = 100;
            for _ in 0..NUM {
                weight[rng.gen_range(0..k)] += 1f64;
            }
            weight.iter_mut().for_each(|x| *x /= NUM as f64);
            weight
        })
        .collect();
    let mut model = EMModel::new(&contexts, &weights, k);
    let mut lk = -1000000f64 * contexts.len() as f64;
    loop {
        let (prev_model, prev_weights) = (model.clone(), weights.clone());
        let prev_lk = model.update(&mut weights, &contexts);
        trace!("LK:{}", prev_lk);
        if (prev_lk - lk) < 0.001 {
            // Revert the update.
            model = prev_model;
            weights = prev_weights;
            break;
        } else {
            lk = prev_lk;
        }
    }
    // for (idx, w) in weights.iter().enumerate() {
    //     let w: Vec<_> = w.iter().map(|x| format!("{:.1}", x)).collect();
    //     let lk = model.lk(&contexts[idx]);
    //     trace!("WEIGHT\t{}\t{}\t{}", idx, w.join("\t"), lk);
    // }
    let read_ids: Vec<_> = reads
        .iter()
        .flat_map(|read| {
            read.nodes
                .iter()
                .filter(|n| n.unit == config.focal)
                .map(|_| read.id)
                .collect::<Vec<_>>()
        })
        .collect();
    let predictions: Vec<_> = contexts
        .iter()
        .zip(weights.iter())
        .zip(read_ids)
        .map(|((ctx, ws), id)| {
            let cluster: usize = ws
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|x| x.0)
                .unwrap();
            (id, ctx.index, cluster as u64)
        })
        .collect();
    trace!("MODEL:{}", model);
    for (context, &(_, _, cl)) in contexts.iter().zip(predictions.iter()) {
        trace!("{}", cl);
        let (forward, _) = model.models[cl as usize].forward_mapping(&context);
        trace!("F:{:?}\t{:?}", forward, context.forward);
        let (backward, _) = model.models[cl as usize].backward_mapping(&context);
        trace!("B:{:?}\t{:?}", backward, context.backward);
    }
    (predictions, lk)
}

// TODO: Add direction.
#[derive(Debug, Clone)]
pub struct Context {
    // The original index of this context.
    index: usize,
    unit: u64,
    cluster: u64,
    forward: Vec<(u64, u64)>,
    backward: Vec<(u64, u64)>,
}

impl Context {
    fn new(read: &EncodedRead, index: usize, use_unit: &HashSet<u64>) -> Self {
        let (unit, cluster) = (read.nodes[index].unit, read.nodes[index].cluster);
        let nodes = read.nodes.iter();
        let forward: Vec<_> = nodes
            .clone()
            .skip(index + 1)
            .map(|n| (n.unit, n.cluster))
            .filter(|n| use_unit.contains(&n.0))
            .collect();
        let backward: Vec<_> = nodes
            .clone()
            .take(index)
            .rev()
            .map(|n| (n.unit, n.cluster))
            .filter(|n| use_unit.contains(&n.0))
            .collect();
        if read.nodes[index].is_forward {
            Self {
                index,
                unit,
                cluster,
                forward,
                backward,
            }
        } else {
            Self {
                index,
                unit,
                cluster,
                forward: backward,
                backward: forward,
            }
        }
    }
}

#[derive(Debug, Clone)]
struct EMModel {
    cluster_num: usize,
    fraction: Vec<f64>,
    models: Vec<RawModel>,
}

impl std::fmt::Display for EMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{}", self.cluster_num)?;
        for (frac, m) in self.fraction.iter().zip(self.models.iter()) {
            writeln!(f, "{:.3}\n{}", frac, m)?;
        }
        Ok(())
    }
}

impl EMModel {
    // Create new model.
    fn new(contexts: &[Context], weights: &[Vec<f64>], cluster_num: usize) -> Self {
        let mut fraction = vec![0f64; cluster_num];
        for weight in weights.iter() {
            fraction.iter_mut().zip(weight).for_each(|(x, y)| *x += y);
        }
        let sum: f64 = fraction.iter().sum();
        fraction.iter_mut().for_each(|x| *x /= sum);
        let mut models = vec![RawModel::new(); cluster_num];
        for (weight, context) in weights.iter().zip(contexts) {
            models.iter_mut().zip(weight).for_each(|(model, &w)| {
                if SMALL < w {
                    model.add(w, context);
                }
            });
        }
        models.iter_mut().for_each(|model| model.normalize());
        Self {
            cluster_num,
            fraction,
            models,
        }
    }
    fn clear(&mut self) {
        self.fraction.iter_mut().for_each(|x| *x = 0f64);
        self.models.iter_mut().for_each(|m| m.clear());
    }
    // Return the **previous** likelihood!!!!!!!!!!
    fn update(&mut self, weights: &mut [Vec<f64>], contexts: &[Context]) -> f64 {
        let mut total_lk = 0f64;
        // Calculate the next weights, mapping position, and likelihood.
        let mapping_positions: Vec<_> = weights
            .iter_mut()
            .zip(contexts.iter())
            .map(|(weight, context)| {
                let (mapping_positions, lk) = self.get_weight(context, weight);
                total_lk += lk;
                mapping_positions
            })
            .collect();
        // Clear the previous model.
        self.clear();
        // Update fraction.
        for weight in weights.iter() {
            for (i, w) in weight.iter().enumerate() {
                self.fraction[i] += w;
            }
        }
        let sum: f64 = self.fraction.iter().sum();
        self.fraction.iter_mut().for_each(|x| *x /= sum);
        // Update paramters.
        let summaries = weights.iter().zip(contexts).zip(mapping_positions);
        for ((weight, context), mapping_positions) in summaries {
            let each_model = weight
                .iter()
                .zip(self.models.iter_mut())
                .zip(mapping_positions);
            for ((&w, model), (forward_map, backward_map)) in each_model {
                if SMALL < w {
                    *model.center.entry(context.cluster).or_default() += w;
                    for (&elm, pos) in context.forward.iter().zip(forward_map) {
                        *model.forward[pos].entry(elm).or_default() += w;
                    }
                    for (&elm, pos) in context.backward.iter().zip(backward_map) {
                        *model.backward[pos].entry(elm).or_default() += w;
                    }
                }
            }
        }
        self.models.iter_mut().for_each(|model| model.normalize());
        total_lk
    }
    // Mapping position is the location of the j-th element to the i-th model.[i][j]
    fn get_weight(
        &self,
        context: &Context,
        weight: &mut [f64],
    ) -> (Vec<(Vec<usize>, Vec<usize>)>, f64) {
        let (lks, mapping_positions): (Vec<_>, Vec<_>) = self
            .models
            .iter()
            .zip(self.fraction.iter())
            .map(|(m, f)| {
                let center_lk = m.center[&context.cluster].max(SMALL).ln();
                let (forward_mapping, lkf) = m.forward_mapping(context);
                let (backward_mapping, lkb) = m.backward_mapping(context);
                let lk = f.max(SMALL).ln() + center_lk + lkf + lkb;
                (lk, (forward_mapping, backward_mapping))
            })
            .unzip();
        let total_lk = logsumexp(&lks);
        weight
            .iter_mut()
            .zip(lks)
            .for_each(|(w, lk)| *w = (lk - total_lk).exp());
        (mapping_positions, total_lk)
    }
    // fn lk(&self, context: &Context) -> f64 {
    //     let lks: Vec<_> = self
    //         .fraction
    //         .iter()
    //         .zip(self.models.iter())
    //         .map(|(f, m)| m.lk(context) + f.max(SMALL).ln())
    //         .collect();
    //     logsumexp(&lks)
    // }
}

const DEL_PROB: f64 = 0.05;
const MATCH_PROB: f64 = 1f64 - DEL_PROB;
const SMALL: f64 = 0.00000000000001;

#[derive(Debug, Clone)]
struct RawModel {
    center: HashMap<u64, f64>,
    forward: Vec<HashMap<(u64, u64), f64>>,
    backward: Vec<HashMap<(u64, u64), f64>>,
}

impl std::fmt::Display for RawModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Center:{:?}", self.center)?;
        for (i, fr) in self.forward.iter().enumerate() {
            let dump: Vec<_> = fr
                .iter()
                .map(|((x, y), val)| format!("({}-{}):{:.3}", x, y, val))
                .collect();
            writeln!(f, "F{}:{}", i, dump.join("\t"))?;
        }
        for (i, b) in self.backward.iter().enumerate() {
            let dump: Vec<_> = b
                .iter()
                .map(|((x, y), val)| format!("({}-{}):{:.3}", x, y, val))
                .collect();
            writeln!(f, "B{}:{}", i, dump.join("\t"))?;
        }
        Ok(())
    }
}

impl RawModel {
    fn new() -> Self {
        Self {
            center: HashMap::new(),
            forward: vec![],
            backward: vec![],
        }
    }
    // Add `context` into self. The position is the same as context.
    fn add(&mut self, w: f64, context: &Context) {
        // Enlarge model if needed.
        let forward_len = context.forward.len().saturating_sub(self.forward.len());
        self.forward
            .extend(std::iter::repeat(HashMap::new()).take(forward_len));
        let backward_len = context.backward.len().saturating_sub(self.backward.len());
        self.backward
            .extend(std::iter::repeat(HashMap::new()).take(backward_len));
        // Add center.
        *self.center.entry(context.cluster).or_default() += w;
        // Add context.
        self.forward
            .iter_mut()
            .zip(context.forward.iter())
            .for_each(|(slot, &elm)| *slot.entry(elm).or_default() += w);
        self.backward
            .iter_mut()
            .zip(context.backward.iter())
            .for_each(|(slot, &elm)| *slot.entry(elm).or_default() += w);
    }
    // fn lk(&self, context: &Context) -> f64 {
    //     let center_lk = self.center[&context.cluster].max(SMALL).ln();
    //     self.forward_mapping(context).1 + center_lk + self.backward_mapping(context).1
    // }
    // Normalize this model to stat model.
    fn normalize(&mut self) {
        let sum: f64 = self.center.values().sum();
        self.center.values_mut().for_each(|x| *x /= sum);
        fn normalize(slots: &mut HashMap<(u64, u64), f64>) {
            let sum: f64 = slots.values().sum();
            slots.values_mut().for_each(|x| *x /= sum);
            slots.retain(|_, prob| SMALL < *prob);
        }
        self.forward.iter_mut().for_each(normalize);
        self.backward.iter_mut().for_each(normalize);
    }
    // Clear all the weights.
    fn clear(&mut self) {
        self.center.values_mut().for_each(|x| *x = 0f64);
        fn clear(slots: &mut HashMap<(u64, u64), f64>) {
            slots.values_mut().for_each(|x| *x = 0f64);
        }
        self.forward.iter_mut().for_each(clear);
        self.backward.iter_mut().for_each(clear);
    }
    fn align_profile(
        elements: &[(u64, u64)],
        profiles: &[HashMap<(u64, u64), f64>],
    ) -> (Vec<usize>, f64) {
        assert!(
            elements.len() <= profiles.len(),
            "{},{}",
            elements.len(),
            profiles.len()
        );
        let minimum = SMALL.ln() * (elements.len() * profiles.len()) as f64;
        let mut dp = vec![vec![minimum; profiles.len() + 1]; elements.len() + 1];
        // True for match, false for deletion.
        let mut is_matched = vec![vec![false; profiles.len() + 1]; elements.len() + 1];
        for j in 0..profiles.len() {
            dp[0][j] = j as f64 * DEL_PROB.ln();
        }
        for (i, elm) in elements.iter().enumerate().map(|(i, p)| (i + 1, p)) {
            for (j, slot) in profiles.iter().enumerate().map(|(j, p)| (j + 1, p)) {
                let match_prob = MATCH_PROB.ln() + slot.get(elm).copied().unwrap_or(0f64);
                let match_score = dp[i - 1][j - 1] + match_prob.max(SMALL).ln();
                let del_score = dp[i][j - 1] + DEL_PROB.ln();
                if del_score < match_score {
                    dp[i][j] = match_score;
                    is_matched[i][j] = true;
                } else {
                    dp[i][j] = del_score;
                    is_matched[i][j] = false;
                }
            }
        }
        // Start from the corner.
        let mut mapping_position = vec![0; elements.len()];
        let mut i = elements.len();
        let (mut j, lk) = dp[i]
            .iter()
            .enumerate()
            .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
            .unwrap();
        while 0 < i || 0 < j {
            if is_matched[i][j] {
                assert!(0 < j && 0 < i, "{},{}", i, j);
                i -= 1;
                j -= 1;
                mapping_position[i] = j;
            } else {
                assert!(0 < j, "{},{}", i, j);
                j -= 1;
            }
        }
        assert_eq!(i, 0);
        assert_eq!(j, 0);
        (mapping_position, *lk)
    }
    // Return the mapping position of each unit in the forward direction of context.
    // In other words, the i-th element would be the index of the self.forward such that
    // the location would give the maximum likelihood.
    fn forward_mapping(&self, context: &Context) -> (Vec<usize>, f64) {
        Self::align_profile(&context.forward, &self.forward)
    }
    // Return the mapping position of each unit in the backward direction of context.
    // In other words, the i-th element would be the index of the self.backward such that
    // the location would give the maximum likelihood.
    fn backward_mapping(&self, context: &Context) -> (Vec<usize>, f64) {
        Self::align_profile(&context.backward, &self.backward)
    }
}

#[cfg(test)]
mod test {}
