use definitions::*;
use rayon::prelude::*;
// const SEED: u64 = 1221;
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
    to_use_offset: bool,
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

pub trait ClusteringCorrection {
    fn correct_clustering_em(
        &mut self,
        repeat_num: usize,
        coverage_thr: usize,
        to_regularize: bool,
    );
    fn correct_clustering_em_on_selected(
        &mut self,
        repeat_num: usize,
        coverage_thr: usize,
        to_regularize: bool,
        selection: &HashSet<u64>,
    );
}

impl ClusteringCorrection for DataSet {
    fn correct_clustering_em(
        &mut self,
        repeat_num: usize,
        coverage_thr: usize,
        to_regularize: bool,
    ) {
        // First, try to "squish" all the units with no variants.
        let to_squish = {
            let mut pvalues = crate::unit_correlation::calc_p_values(self, coverage_thr as u32);
            pvalues.retain(|_, pvalue| 0.05 < *pvalue);
            pvalues
        };
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .filter(|n| to_squish.contains_key(&n.unit))
            .for_each(|n| n.cluster = 0);
        // Then, usual correction.
        // Maybe we can cluster with compressed .. no.
        // This is because the consensus would be very dirty.
        let targets: HashSet<_> = self
            .selected_chunks
            .iter()
            .filter_map(|c| (0.1 < c.score).then(|| c.id))
            .collect();
        debug!("EM\tSquish\tTarget\tTotal");
        debug!(
            "EM\t{}\t\t{}\t{}",
            to_squish.len(),
            targets.len(),
            self.selected_chunks.len(),
        );
        self.correct_clustering_em_on_selected(repeat_num, coverage_thr, to_regularize, &targets);
    }
    fn correct_clustering_em_on_selected(
        &mut self,
        repeat_num: usize,
        coverage_thr: usize,
        to_regularize: bool,
        selection: &HashSet<u64>,
    ) {
        let result: Vec<_> = self
            .selected_chunks
            .par_iter()
            .filter(|c| selection.contains(&c.id))
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
                let (new_clustering, score, new_k) = (1..=k)
                    .flat_map(|k| match k {
                        1 => std::iter::repeat(k).take(1),
                        _ => std::iter::repeat(k).take(repeat_num),
                    })
                    .enumerate()
                    .map(|(i, k)| {
                        let seed = unit_id * (i * k) as u64;
                        let config =
                            Config::new(repeat_num, seed, k, unit_id, to_regularize, coverage_thr);
                        em_clustering(&reads, &config)
                    })
                    .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
                    .unwrap();
                trace!("CLUSTERED\t{}\t{}\t{:.3}", unit_id, new_k, score);
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

/// Return the id of the read, the position at that read, and the clusters predicted.
pub fn em_clustering(
    reads: &[&EncodedRead],
    config: &Config,
) -> (Vec<(u64, usize, u64)>, f64, usize) {
    let mut unit_counts: HashMap<_, usize> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *unit_counts.entry(node.unit).or_default() += 1;
        }
    }
    unit_counts.retain(|_, c| config.coverage_thr < *c);
    // let use_units: HashSet<_> = unit_counts
    //     .iter()
    //     .filter(|&(_, &c)| c > config.coverage_thr)
    //     .map(|(&x, _)| x)
    //     .collect();
    let contexts: Vec<_> = {
        let mut buffer = vec![];
        for read in reads.iter() {
            for index in 0..read.nodes.len() {
                if read.nodes[index].unit == config.focal {
                    // buffer.push(Context::new(read, index, &use_units));
                    buffer.push(Context::new(read, index, &unit_counts));
                }
            }
        }
        buffer
    };
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(config.seed);
    let cluster_num = config.cluster_num;
    //let (asn, lk, offset) = em_clustering_inner(&contexts, cluster_num, &mut rng);
    let (asn, lk, offset) = simple_clustering_inner(&contexts, cluster_num, &mut rng);
    let lk = if config.to_use_offset {
        lk - offset
    } else {
        lk
    };
    (asn, lk, cluster_num)
}

pub fn initialize_weights<R: Rng>(contexts: &[Context], k: usize, rng: &mut R) -> Vec<Vec<f64>> {
    fn norm(c1: &Context, c2: &Context) -> f64 {
        let center_match = (c1.cluster != c2.cluster) as u32;
        let (f, flen) = c1
            .forward
            .iter()
            .zip(c2.forward.iter())
            .fold((0, 0), |(acc, len), (c1, c2)| {
                (acc + (c1 != c2) as u32, len + 1)
            });
        let (b, blen) = c1
            .backward
            .iter()
            .zip(c2.backward.iter())
            .fold((0, 0), |(acc, len), (c1, c2)| {
                (acc + (c1 != c2) as u32, len + 1)
            });
        (center_match + b + f) as f64 / (1 + flen + blen) as f64
    }
    let mut centers: Vec<&Context> = vec![];
    let indices: Vec<_> = (0..contexts.len()).collect();
    // Choosing centers.
    use rand::seq::SliceRandom;
    while centers.len() < k as usize {
        // calculate distance to the most nearest centers.
        let mut dists: Vec<_> = contexts
            .iter()
            .map(|ctx| {
                centers
                    .iter()
                    .map(|c| norm(ctx, c))
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap_or(1f64)
                    .powi(2)
            })
            .collect();
        let total: f64 = dists.iter().sum();
        if total > 0.00001 {
            dists.iter_mut().for_each(|x| *x /= total);
        } else {
            dists.iter_mut().for_each(|x| *x = 1f64);
        }
        let idx = match indices.choose_weighted(rng, |&idx| dists[idx]) {
            Ok(res) => *res,
            Err(why) => {
                debug!("{},{}", total, centers.len());
                for ctx in contexts {
                    debug!("{}", ctx);
                }
                panic!("{:?},{:?}", dists, why)
            }
        };
        centers.push(&contexts[idx]);
    }
    contexts
        .iter()
        .map(|ctx| {
            let idx = centers
                .iter()
                .enumerate()
                .map(|(i, c)| (i, norm(c, ctx)))
                .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
                .unwrap()
                .0;
            let mut weights = vec![100f64; k];
            weights[idx] += 100f64;
            let sum: f64 = weights.iter().sum();
            weights.iter_mut().for_each(|x| *x /= sum);
            weights
        })
        .collect()
}
pub fn em_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    rng: &mut R,
) -> (Vec<(u64, usize, u64)>, f64, f64) {
    let mut weights: Vec<_> = match rng.gen_bool(0.5) {
        true => contexts
            .iter()
            .map(|_| {
                let mut ws = vec![0f64; k];
                for _ in 0..100 {
                    ws[rng.gen_range(0..k)] += 1f64;
                }
                ws.iter_mut().for_each(|x| *x /= 100f64);
                ws
            })
            .collect(),
        false => initialize_weights(contexts, k, rng),
    };
    let mut model = EMModel::new(contexts, &weights, k);
    let mut lk: f64 = contexts
        .iter()
        .map(|ctx| model.get_weight(ctx, &mut []).1)
        .sum();
    trace!("LK:{}", lk);
    loop {
        let (prev_model, prev_weights) = (model.clone(), weights.clone());
        model.update(&mut weights, contexts);
        let next_lk = contexts
            .iter()
            .map(|ctx| model.get_weight(ctx, &mut []).1)
            .sum();
        trace!("LK:{}", next_lk);
        if (next_lk - lk) < 0.001 {
            model = prev_model;
            weights = prev_weights;
            break;
        } else {
            lk = next_lk;
        }
    }
    let predictions = predict_from_weights(contexts, &weights, k, rng);
    trace!("MODEL:{}", model);
    let (flen, blen) = contexts
        .iter()
        .map(|c| (c.forward.len(), c.backward.len()))
        .fold((0, 0), |(f, b), (x, y)| (f.max(x), b.max(y)));
    for (weight, context) in weights.iter().zip(contexts.iter()) {
        let (cluster, _) = max_and_position(weight);
        let unit = context.unit;
        trace!(
            "DUMP\t{}\t{}\t{}\t{}",
            unit,
            context.id,
            cluster,
            context.dump_with(flen, blen)
        );
    }
    // Num of parameters.
    let offset = {
        let mut units: HashMap<_, HashSet<_>> = HashMap::new();
        for ctx in contexts.iter() {
            units.entry(ctx.unit).or_default().insert(ctx.cluster);
            for &(u, c) in ctx.forward.iter().chain(ctx.backward.iter()) {
                units.entry(u).or_default().insert(c);
            }
        }
        let model_param: usize = units.values().map(|x| x.len() - 1).sum();
        (model_param * k + k - 1) as f64
    };
    trace!("FinalLK\t{:.4}", lk);
    (predictions, lk, offset)
}

#[derive(Debug, Clone)]
pub struct Context {
    // The ID of the original encoded read.
    id: u64,
    // The original index of this context.
    index: usize,
    unit: u64,
    cluster: u64,
    forward: Vec<(u64, u64)>,
    backward: Vec<(u64, u64)>,
}

impl std::fmt::Display for Context {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let forward: Vec<_> = self
            .forward
            .iter()
            .map(|(u, c)| format!("({},{})", u, c))
            .collect();
        let backward: Vec<_> = self
            .backward
            .iter()
            .map(|(u, c)| format!("({},{})", u, c))
            .collect();
        write!(
            f,
            "{}\t({},{})\t{}",
            forward.join("-"),
            self.unit,
            self.cluster,
            backward.join("-")
        )
    }
}

impl Context {
    pub fn dump_with(&self, fdigit: usize, bdigit: usize) -> String {
        let mut f_slots = vec!["---".to_string(); fdigit];
        let mut b_slots = vec!["---".to_string(); bdigit];
        for (i, (u, c)) in self.forward.iter().enumerate().take(fdigit) {
            f_slots[i] = format!("{}-{}", u, c);
        }
        for (i, (u, c)) in self.backward.iter().enumerate().take(bdigit) {
            b_slots[bdigit - i - 1] = format!("{}-{}", u, c);
        }
        format!(
            "{}\t{}-{}\t{}",
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
        cluster: u64,
        forward: Vec<(u64, u64)>,
        backward: Vec<(u64, u64)>,
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
    // fn new(read: &EncodedRead, index: usize, use_unit: &HashSet<u64>) -> Self {
    fn new(read: &EncodedRead, index: usize, unit_counts: &HashMap<u64, usize>) -> Self {
        let (unit, cluster) = (read.nodes[index].unit, read.nodes[index].cluster);
        let nodes = read.nodes.iter();
        let forward: Vec<_> = nodes
            .clone()
            .skip(index + 1)
            .map(|n| (n.unit, n.cluster))
            //.filter(|n| use_unit.contains(&n.0))
            .filter(|n| unit_counts.contains_key(&n.0))
            .collect();
        let backward: Vec<_> = nodes
            .clone()
            .take(index)
            .rev()
            .map(|n| (n.unit, n.cluster))
            //.filter(|n| use_unit.contains(&n.0))
            .filter(|n| unit_counts.contains_key(&n.0))
            .collect();
        if read.nodes[index].is_forward {
            Self {
                id: read.id,
                index,
                unit,
                cluster,
                forward,
                backward,
            }
        } else {
            Self {
                id: read.id,
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
    center_unit: u64,
    cluster_num: HashMap<u64, u64>,
    fraction: Vec<f64>,
    models: Vec<RawModel>,
}

impl std::fmt::Display for EMModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{}", self.cluster_num[&self.center_unit])?;
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
        let (fm, bm) = contexts
            .iter()
            .map(|ctx| (ctx.forward.len(), ctx.backward.len()))
            .fold((0, 0), |(f, b), (x, y)| (f.max(x), b.max(y)));
        let mut models = vec![RawModel::new(fm, bm); cluster_num];
        for (weight, context) in weights.iter().zip(contexts) {
            models.iter_mut().zip(weight).for_each(|(model, &w)| {
                if SMALL < w {
                    model.add(w, context);
                }
            });
        }
        models.iter_mut().for_each(|model| model.normalize());
        let center_unit = contexts[0].unit;
        let mut cluster_num: HashMap<_, u64> =
            std::iter::once((contexts[0].unit, cluster_num as u64)).collect();
        for context in contexts {
            for &(unit, cluster) in context.forward.iter() {
                let x = cluster_num.entry(unit).or_default();
                *x = (*x).max(cluster + 1);
            }
            for &(unit, cluster) in context.backward.iter() {
                let x = cluster_num.entry(unit).or_default();
                *x = (*x).max(cluster + 1);
            }
        }
        Self {
            center_unit,
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
        {
            for weight in weights.iter() {
                for (i, w) in weight.iter().enumerate() {
                    self.fraction[i] += w;
                }
            }
            let sum: f64 = self.fraction.iter().sum();
            self.fraction.iter_mut().for_each(|x| *x /= sum);
        }
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
    #[allow(clippy::type_complexity)]
    fn get_weight(
        &self,
        context: &Context,
        weight: &mut [f64],
    ) -> (Vec<(Vec<usize>, Vec<usize>)>, f64) {
        let (lks, mapping_positions): (Vec<_>, Vec<_>) = self
            .models
            .iter()
            .zip(self.fraction.iter())
            .map(|(m, &f)| {
                // TODO: change here.
                let center_lk = (MISM_PROB * (self.cluster_num[&context.unit] as f64).recip()
                    + (1f64 - MISM_PROB) * m.center.get(&context.cluster).unwrap_or(&SMALL))
                .ln();
                let (forward_mapping, lkf) = m.forward_mapping(context);
                let (backward_mapping, lkb) = m.backward_mapping(context);
                let lk = f.max(SMALL).ln() + center_lk + lkf + lkb;
                assert!(
                    !lk.is_nan(),
                    "{}\n{}\n{}\t{}\t{}",
                    self,
                    context,
                    f,
                    lkf,
                    lkb
                );
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
}

const DEL_PROB: f64 = 0.02;
// const DEL_PROB: f64 = 0.1;
const MISM_PROB: f64 = 0.02;
const MATCH_PROB: f64 = 1f64 - DEL_PROB - MISM_PROB;
const SMALL: f64 = 0.00001;

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
    // fn num_param(&self) -> usize {
    //     let center = self.center.len().max(1) - 1;
    //     let count_param =
    //         |hm: &HashMap<(u64, u64), f64>| hm.values().filter(|&&x| 0.01 < x).count().max(1) - 1;
    //     let forward: usize = self.forward.iter().map(count_param).sum();
    //     let backward: usize = self.backward.iter().map(count_param).sum();
    //     center + forward + backward
    // }
    fn new(forward_len: usize, backward_len: usize) -> Self {
        Self {
            center: HashMap::new(),
            forward: vec![HashMap::new(); forward_len],
            backward: vec![HashMap::new(); backward_len],
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
    // Normalize this model to stat model.
    fn normalize(&mut self) {
        let sum: f64 = self.center.values().sum();
        if SMALL < sum {
            self.center.values_mut().for_each(|x| *x /= sum);
        }
        fn normalize(slots: &mut HashMap<(u64, u64), f64>) {
            let sum: f64 = slots.values().sum();
            slots.values_mut().for_each(|x| *x /= sum);
            slots.retain(|_, prob| 0.01 < *prob);
            // slots.retain(|_, prob| SMALL < *prob);
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
        let minimum = SMALL.ln() * (elements.len() * profiles.len()) as f64;
        let mut dp = vec![vec![minimum; profiles.len() + 1]; elements.len() + 1];
        // True for match, false for deletion.
        let mut is_matched = vec![vec![false; profiles.len() + 1]; elements.len() + 1];
        for j in 0..profiles.len() {
            dp[0][j] = j as f64 * DEL_PROB.ln();
        }
        for (i, elm) in elements.iter().enumerate().map(|(i, p)| (i + 1, p)) {
            for (j, slot) in profiles.iter().enumerate().map(|(j, p)| (j + 1, p)) {
                // If there is no unit, then we use uniformal distribution.
                let match_prob = MATCH_PROB * slot.get(elm).copied().unwrap_or(0.5);
                let mismatch_prob = MISM_PROB * 0.5;
                // Modified.
                let match_score = dp[i - 1][j - 1] + (match_prob + mismatch_prob).max(SMALL).ln();
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
            .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
            .unwrap();
        while 0 < i && 0 < j {
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
        while 0 < i {
            i -= 1;
            mapping_position[i] = j;
        }
        while 0 < j {
            j -= 1;
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

fn max_and_position(ws: &[f64]) -> (usize, f64) {
    ws.iter()
        .copied()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap()
}

// Sometimes un-informative reads would be clustered into the largest
// cluster. It is indeed a ML-predictor, but not good (Consider that
// all such reads are into the largest cluster. Then, the relative fraction of the
// largest cluster is no longer proportional to the EM-estimator!).
fn predict_from_weights<R: Rng>(
    contexts: &[Context],
    weights: &[Vec<f64>],
    _k: usize,
    _rng: &mut R,
) -> Vec<(u64, usize, u64)> {
    // let choises: Vec<_> = (0..k).collect();
    contexts
        .iter()
        .zip(weights.iter())
        .map(|(ctx, ws)| {
            let (cluster, _max_w) = max_and_position(ws);
            // let cluster = if max_w < 0.9 {
            //     cluster as u64
            // } else {
            //     use rand::seq::SliceRandom;
            //     *choises.choose_weighted(rng, |&k| ws[k]).unwrap_or(&0) as u64
            // };
            (ctx.id, ctx.index, cluster as u64)
        })
        .collect()
}

pub fn simple_clustering_inner<R: Rng>(
    contexts: &[Context],
    k: usize,
    rng: &mut R,
) -> (Vec<(u64, usize, u64)>, f64, f64) {
    // ID of this trial.
    let id: u64 = rng.gen::<u64>() % 1_000_000;
    let mut weights: Vec<_> = match rng.gen_bool(0.5) {
        true => contexts
            .iter()
            .map(|_| {
                let mut ws = vec![0f64; k];
                for _ in 0..100 {
                    ws[rng.gen_range(0..k)] += 1f64;
                }
                ws.iter_mut().for_each(|x| *x /= 100f64);
                ws
            })
            .collect(),
        false => initialize_weights(contexts, k, rng),
    };
    let mut model = SimpleModel::new(contexts, &weights, k);
    let mut lk: f64 = contexts.iter().map(|ctx| model.lk(ctx)).sum();
    trace!("LK\t{}\t{}", id, lk);
    for _ in 0..20 {
        // loop {
        model.update(&mut weights, contexts);
        model = SimpleModel::new(contexts, &weights, k);
        let next_lk = contexts.iter().map(|ctx| model.lk(ctx)).sum();
        trace!("LK\t{}\t{}", id, next_lk);
        if (next_lk - lk) < 0.001 {
            break;
        } else {
            lk = next_lk;
        }
    }
    trace!("MODEL\t{}\n{}", id, model);
    let (flen, blen) = contexts
        .iter()
        .map(|c| (c.forward.len(), c.backward.len()))
        .fold((0, 0), |(f, b), (x, y)| (f.max(x), b.max(y)));
    for ws in weights.iter() {
        trace!("WEIGHTS\t{}\t{:?}", id, ws);
    }
    let predictions = predict_from_weights(contexts, &weights, k, rng);
    for ((_, _, cluster), context) in predictions.iter().zip(contexts.iter()) {
        let unit = context.unit;
        let dump = context.dump_with(flen, blen);
        trace!(
            "DUMP\t{}\t{}\t{}\t{}\t{}",
            id,
            unit,
            context.id,
            cluster,
            dump,
        );
    }
    let num_param = model.num_parameters() as f64;
    trace!("Final\t{}\t{}\t{:.4}\t{}", id, k, lk - num_param, num_param);
    (predictions, lk, num_param)
}

#[derive(Debug, Clone)]
struct SimpleModel {
    fraction: Vec<f64>,
    consensus: Vec<Consensus>,
}

impl std::fmt::Display for SimpleModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dump: Vec<_> = self
            .fraction
            .iter()
            .zip(self.consensus.iter())
            .map(|(frac, m)| format!("{:.3}\t{}", frac, m))
            .collect();
        write!(f, "{}", dump.join("\n"))
    }
}

impl SimpleModel {
    // Create new model.
    fn new(contexts: &[Context], weights: &[Vec<f64>], cluster_num: usize) -> Self {
        let mut fraction = vec![0f64; cluster_num];
        for weight in weights.iter() {
            fraction.iter_mut().zip(weight).for_each(|(x, y)| *x += y);
        }
        let sum: f64 = fraction.iter().sum();
        fraction.iter_mut().for_each(|x| *x /= sum);
        let flen = contexts
            .iter()
            .map(|x| x.forward.len())
            .max()
            .unwrap_or_default();
        let blen = contexts
            .iter()
            .map(|x| x.forward.len())
            .max()
            .unwrap_or_default();
        let consensus: Vec<_> = (0..cluster_num)
            .map(|cl| {
                let ws: Vec<_> = weights.iter().map(|ws| ws[cl]).collect();
                Consensus::new(contexts, &ws, (flen, blen))
            })
            .collect();
        Self {
            fraction,
            consensus,
        }
    }
    fn num_parameters(&self) -> usize {
        let mut used_units: HashSet<_> = HashSet::new();
        for cons in self.consensus.iter() {
            used_units.insert(cons.center);
            used_units.extend(cons.forward.iter().copied());
            used_units.extend(cons.backward.iter().copied());
        }
        let model_param = used_units.len();
        // let model_param = self.consensus.iter().map(|x| x.num_parameters()).max();
        // let model_param = self.consensus.len() * model_param.unwrap_or_default();
        model_param + self.fraction.len() - 1
    }
    // Update weights.
    fn update(&mut self, weights: &mut [Vec<f64>], contexts: &[Context]) {
        for (weight, context) in weights.iter_mut().zip(contexts.iter()) {
            let _ = self.get_weight(context, weight);
        }
    }
    #[allow(clippy::type_complexity)]
    fn get_weight(&self, context: &Context, weight: &mut [f64]) -> (Vec<(Vec<u8>, Vec<u8>)>, f64) {
        let (scores, mapping_positions) = self.scores_and_mappings(context);
        let score = logsumexp(&scores);
        weight
            .iter_mut()
            .zip(scores.iter())
            .for_each(|(w, s)| *w = (s - score).exp());
        (mapping_positions, score)
    }
    fn scores_and_mappings(&self, context: &Context) -> (Vec<f64>, Vec<(Vec<u8>, Vec<u8>)>) {
        self.consensus
            .iter()
            .zip(self.fraction.iter())
            .map(|(m, &f)| {
                let (forward_map, backward_map, score) = m.map(context);
                let score = score as f64 + f.ln();
                (score, (forward_map, backward_map))
            })
            .unzip()
    }
    fn lk(&self, context: &Context) -> f64 {
        let (scores, _) = self.scores_and_mappings(context);
        logsumexp(&scores)
        // let score = logsumexp(&scores);
        // let scores: Vec<_> = scores.iter().map(|x| format!("{:.3}", x)).collect();
        // trace!("SCORE\t{}\t{}", scores.join(","), score);
        // score
    }
}
#[derive(Debug, Clone)]
struct Consensus {
    center: (u64, u64),
    forward: Vec<(u64, u64)>,
    backward: Vec<(u64, u64)>,
}

impl std::fmt::Display for Consensus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let backward = self
            .backward
            .iter()
            .map(|(u, c)| format!("{}-{}", u, c))
            .rev();
        let forward = self.forward.iter().map(|(u, c)| format!("{}-{}", u, c));
        let center = std::iter::once(format!("{}-{}", self.center.0, self.center.1));
        let dump: Vec<_> = backward.chain(center).chain(forward).collect();
        write!(f, "{}", dump.join("\t"))
    }
}

// Alignment parameters.
const MAT: i32 = 1;
const MIS_CL: i32 = -1;
const MIS_UN: i32 = -1000;
const INS: i32 = -10;
const DEL: i32 = -5;

impl Consensus {
    // fn num_parameters(&self) -> usize {
    //     // 1 for the center unit.
    //     1 + self.forward.len() + self.backward.len()
    // }
    // flen = max length of the forward, ble = max length of the minimum.
    // Cut consensus if the `coverage` is dropped below 1.
    // TODO: Maybe we should align reads onto the consensus to update...?
    fn new(xs: &[Context], ws: &[f64], (flen, blen): (usize, usize)) -> Self {
        // Column sum, take the max value.
        let mut center: HashMap<_, f64> = HashMap::new();
        let mut forward: Vec<HashMap<(u64, u64), f64>> = vec![HashMap::new(); flen];
        let mut backward: Vec<HashMap<(u64, u64), f64>> = vec![HashMap::new(); blen];
        for (w, ctx) in ws.iter().zip(xs).filter(|&(&w, _)| SMALL < w) {
            *center.entry(ctx.cluster).or_default() += w;
            for (fslot, &f) in forward.iter_mut().zip(ctx.forward.iter()) {
                *fslot.entry(f).or_default() += w;
            }
            for (bslot, &b) in backward.iter_mut().zip(ctx.backward.iter()) {
                *bslot.entry(b).or_default() += w;
            }
        }
        let center: u64 = center
            .iter()
            .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
            .map(|x| *(x.0))
            .unwrap_or_default();
        let center = (xs[0].unit, center);
        let forward: Vec<_> = forward
            .iter()
            .filter_map(|slots| slots.iter().max_by(|x, y| (x.1).partial_cmp(y.1).unwrap()))
            .map(|x| x.0)
            .copied()
            .collect();
        let backward: Vec<(u64, u64)> = backward
            .iter()
            .filter_map(|slots| slots.iter().max_by(|x, y| (x.1).partial_cmp(y.1).unwrap()))
            .map(|x| x.0)
            .copied()
            .collect();
        Self {
            center,
            forward,
            backward,
        }
    }
    // Align query to context, return the alignment positions of the query.
    // The match score is 1, mismatch score for different cluster is -1, different unit is -1000, and the
    // insertion score is -10, deletion score is -5.
    // TODO: Justify these alignment scores.
    // 0->mat, 1->Del to consensus 2-> Ins to consensus, 3-> Mismatch
    fn map(&self, context: &Context) -> (Vec<u8>, Vec<u8>, i32) {
        let (fm, fscore) = Self::align(&self.forward, &context.forward);
        let (bm, bscore) = Self::align(&self.backward, &context.backward);
        let center = match self.center.1 == context.cluster {
            true => MAT,
            false => MIS_CL,
        };
        (fm, bm, fscore + bscore + center)
    }
    // Align query to cons.
    // It is "tail-free" alignment. In other words,
    // the un-aligned reigon in the end of the `cons` would not be penalized.
    fn align(cons: &[(u64, u64)], query: &[(u64, u64)]) -> (Vec<u8>, i32) {
        let minimum = MIS_UN * (cons.len() + query.len()) as i32;
        let mut dp = vec![vec![minimum; cons.len() + 1]; query.len() + 1];
        for (i, row) in dp.iter_mut().enumerate() {
            row[0] = i as i32 * INS;
        }
        for j in 0..cons.len() {
            dp[0][j] = j as i32 * DEL;
        }
        for (i, q) in query.iter().enumerate().map(|(i, q)| (i + 1, q)) {
            // Because the last deletion would not be penalized.
            let del = if i == query.len() { 0 } else { DEL };
            for (j, c) in cons.iter().enumerate().map(|(j, c)| (j + 1, c)) {
                let mat_score = match (q, c) {
                    (x, y) if x == y => MAT,
                    ((qu, _), (cu, _)) if qu == cu => MIS_CL,
                    _ => MIS_UN,
                };
                dp[i][j] = (dp[i - 1][j - 1] + mat_score)
                    .max(dp[i - 1][j] + INS)
                    .max(dp[i][j - 1] + del);
            }
        }
        // Start from the corner.
        let mut i = query.len();
        let (mut j, score) = dp[i]
            .iter()
            .enumerate()
            .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
            .unwrap();
        let mut ops = vec![1; cons.len() - j];
        while 0 < i && 0 < j {
            let (q, c) = (query[i - 1], cons[j - 1]);
            let mat_score = match (q, c) {
                (x, y) if x == y => MAT,
                ((qu, _), (cu, _)) if qu == cu => MIS_CL,
                _ => MIS_UN,
            };
            let current = dp[i][j];
            let del = if i == query.len() { 0 } else { DEL };
            if current == dp[i - 1][j - 1] + mat_score {
                i -= 1;
                j -= 1;
                let op = if mat_score == MAT { 0 } else { 3 };
                ops.push(op);
            } else if current == dp[i][j - 1] + del {
                j -= 1;
                ops.push(1);
            } else {
                assert_eq!(current, dp[i - 1][j] + INS);
                i -= 1;
                ops.push(2);
            }
        }
        ops.extend(std::iter::repeat(1).take(j));
        ops.extend(std::iter::repeat(2).take(i));
        (ops, *score)
    }
}
