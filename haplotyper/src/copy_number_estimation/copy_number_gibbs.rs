use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use std::collections::HashMap;
// Each edge has its direction to either "plus" or "minus" direction of the node.
type Edge = (usize, bool, usize, bool, f64);
type GbsEdge = (usize, bool, usize, bool, u64);
use log::*;

// Out of the haploid coverage, how much fraction would be observed in the 0-coverage nodes.
// For example, if the haploid coverage is 20, we assume that there would be 20 * 0.25 = 5 occurence of edge
// even in the zero-copy elements.
const ERROR_FRAC: f64 = 0.25;
// If true, condidence is set to the MAX_CONFIDENCE from the beggining.
const CONST_CONFIDENCE: bool = false;
const MAX_CONFIDENCE: f64 = 0.95;
// sample copy-number estimation in  `BURN_IN` times from confidence=0 to condidence=MAX_CONFIDENCE,
// then keep sampling `BURN_IN` times at confidence=MAX_CONFIDENCE to reach the stationaly distribution.
const BURN_IN: usize = 2_000;
// After burn-in, sample `SAMPLE_LEN` to get max a posterior estimation.
const SAMPLE_LEN: usize = 2_000;

#[derive(Debug, Clone)]
pub struct GibbsSampler {
    nodes: Vec<u64>,
    edges: Vec<GbsEdge>,
    // i->(indices of edges with (i,true), indices of edges with (i,false));
    node_terminals: Vec<(Vec<usize>, Vec<usize>)>,
    haploid_coverage: f64,
}

impl std::fmt::Display for GibbsSampler {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "N:{}\tE:{}", self.nodes.len(), self.edges.len())
    }
}

/// The structure to configure the parameters in copy number estimation.
#[derive(Debug, Clone)]
pub struct Config {
    seed: u64,
}

impl Config {
    pub fn new(seed: u64) -> Self {
        Self { seed }
    }
}

impl std::default::Default for Config {
    fn default() -> Self {
        Self { seed: 24309 }
    }
}

impl GibbsSampler {
    pub fn new(nodes: &[f64], edges: &[Edge], cov: f64) -> Self {
        let mut node_terminals = vec![(vec![], vec![]); nodes.len()];
        for (idx, &(from, fd, to, td, _)) in edges.iter().enumerate() {
            match fd {
                true => node_terminals[from].0.push(idx),
                false => node_terminals[from].1.push(idx),
            }
            match td {
                true => node_terminals[to].0.push(idx),
                false => node_terminals[to].1.push(idx),
            }
        }
        let nodes: Vec<_> = nodes.iter().map(|x| x.round() as u64).collect();
        let edges: Vec<_> = edges
            .iter()
            .map(|&(f, fd, t, td, w)| (f, fd, t, td, w.round() as u64))
            .collect();
        Self {
            nodes,
            edges,
            node_terminals,
            haploid_coverage: cov,
        }
    }
    /// Sample copy number of nodes and edges.
    pub fn sample_copy_numer(&self, config: &Config) -> (Vec<usize>, Vec<usize>) {
        let mut node_cp: Vec<usize> = self
            .nodes
            .iter()
            .map(|&x| (x as f64 / self.haploid_coverage).round() as usize)
            .collect();
        let mut edge_cp: Vec<_> = self
            .edges
            .iter()
            .map(|x| (x.4 as f64 / self.haploid_coverage).round() as usize)
            .collect();
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(config.seed);
        for i in 0..BURN_IN {
            let confidence = match CONST_CONFIDENCE {
                true => MAX_CONFIDENCE,
                false => i as f64 * MAX_CONFIDENCE / BURN_IN as f64,
            };
            self.update_nodes(&mut node_cp, &edge_cp, confidence, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, confidence, &mut rng);
        }
        // spin loop to reach stationaly distribution.
        for _ in 0..BURN_IN {
            self.update_nodes(&mut node_cp, &edge_cp, MAX_CONFIDENCE, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, MAX_CONFIDENCE, &mut rng);
        }
        // Let's calculate the posterior distributions on each nodes/edges.
        let mut node_cp_dist: Vec<_> = node_cp.iter().map(|&c| vec![0; 2 * (c + 1)]).collect();
        let mut edge_cp_dist: Vec<_> = edge_cp.iter().map(|&c| vec![0; 2 * (c + 1)]).collect();
        for _ in 0..SAMPLE_LEN {
            self.update_nodes(&mut node_cp, &edge_cp, MAX_CONFIDENCE, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, MAX_CONFIDENCE, &mut rng);
            for (buf, &x) in node_cp_dist.iter_mut().zip(node_cp.iter()) {
                if let Some(slot) = buf.get_mut(x) {
                    *slot += 1;
                }
            }
            for (buf, &x) in edge_cp_dist.iter_mut().zip(edge_cp.iter()) {
                if let Some(slot) = buf.get_mut(x) {
                    *slot += 1;
                }
            }
        }
        let argmax = |buf: &Vec<u32>| buf.iter().enumerate().max_by_key(|x| x.1).unwrap().0;
        let node_cp: Vec<_> = node_cp_dist.iter().map(argmax).collect();
        let edge_cp: Vec<_> = edge_cp_dist.iter().map(argmax).collect();
        (node_cp, edge_cp)
    }
    fn update_nodes<R: Rng>(
        &self,
        node_cp: &mut [usize],
        edge_cp: &[usize],
        confidence: f64,
        rng: &mut R,
    ) {
        let mut order: Vec<_> = (0..self.nodes.len()).collect();
        order.shuffle(rng);
        for node in order {
            self.update_node(node, &mut node_cp[node], edge_cp, confidence, rng);
        }
    }
    fn update_node<R: Rng>(
        &self,
        node: usize,
        node_cp: &mut usize,
        edge_cp: &[usize],
        confidence: f64,
        rng: &mut R,
    ) {
        let mut cps: Vec<usize> = vec![];
        let (down, up) = &self.node_terminals[node];
        if !down.is_empty() {
            cps.push(down.iter().map(|&j| edge_cp[j]).sum());
        }
        if !up.is_empty() {
            cps.push(up.iter().map(|&j| edge_cp[j]).sum());
        }
        let w = self.nodes[node];
        *node_cp = choose_copy_num(w, &cps, self.haploid_coverage, confidence, rng);
    }
    fn update_edges<R: Rng>(
        &self,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        confidence: f64,
        rng: &mut R,
    ) {
        let mut order: Vec<_> = (0..self.edges.len()).collect();
        order.shuffle(rng);
        for edge in order {
            self.update_edge(edge, node_cp, edge_cp, confidence, rng);
        }
    }
    fn update_edge<R: Rng>(
        &self,
        edge: usize,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        confidence: f64,
        rng: &mut R,
    ) {
        let (from, fd, to, td, w) = self.edges[edge];
        let from = self.edge_copy_num(edge, from, fd, node_cp, edge_cp);
        let to = self.edge_copy_num(edge, to, td, node_cp, edge_cp);
        let cps = [from, to];
        edge_cp[edge] = choose_copy_num(w, &cps, self.haploid_coverage, confidence, rng);
    }
    fn edge_copy_num(
        &self,
        edge: usize,
        node: usize,
        is_plus: bool,
        node_cp: &[usize],
        edge_cp: &[usize],
    ) -> usize {
        let (plus, minus) = &self.node_terminals[node];
        let edge_cps_sum: usize = match is_plus {
            true => plus.iter().map(|&e| edge_cp[e]).sum(),
            false => minus.iter().map(|&e| edge_cp[e]).sum(),
        };
        (node_cp[node] + edge_cp[edge]).saturating_sub(edge_cps_sum)
    }
}

// Poisson(obs|copy_num*coverage)~Norm(obs|mean=var=copy_num*coverage)
// If copy number is zero, It is Norm(obs|mean=0,var=coverage*ERROR_FRAC), this is rough heuristics.
fn poisson(obs: u64, copy_num: usize, coverage: f64) -> f64 {
    let lambda = match copy_num {
        0 => coverage * ERROR_FRAC,
        _ => copy_num as f64 * coverage,
    };
    let denom: f64 = (1..obs + 1).map(|i| (i as f64).ln()).sum();
    (obs as f64 * lambda.ln() - lambda - denom).exp()
}

fn choose_copy_num<R: Rng>(
    w: u64,
    cps: &[usize],
    hap_cov: f64,
    confidence: f64,
    rng: &mut R,
) -> usize {
    let mut choises = vec![];
    for &cp in cps {
        if cp == 0 {
            let trust_prob = 0.5 + confidence / 2f64;
            choises.push((cp, trust_prob * poisson(w, cp, hap_cov)));
            choises.push((cp + 1, (1f64 - trust_prob) * poisson(w, cp + 1, hap_cov)));
        } else {
            let trust_prob = 3f64.recip() + 2f64 / 3f64 * confidence;
            choises.push((cp, trust_prob * poisson(w, cp, hap_cov)));
            let minus_one = (1f64 - trust_prob) / 2f64 * poisson(w, cp - 1, hap_cov);
            choises.push((cp - 1, minus_one));
            let plus_one = (1f64 - trust_prob) / 2f64 * poisson(w, cp + 1, hap_cov);
            choises.push((cp + 1, plus_one));
        }
    }
    if !choises.is_empty() {
        let sum: f64 = choises.iter().map(|x| x.1).sum();
        choises.choose_weighted(rng, |&(_, w)| w / sum).unwrap().0
    } else {
        // This node is isolated.
        let copy_num = (w as f64 / hap_cov).floor() as usize;
        (copy_num.saturating_sub(2)..copy_num + 2)
            .map(|cp| (cp, poisson(w, cp, hap_cov)))
            .max_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
            .unwrap()
            .0
    }
}

const CHAIN_NUM: u64 = 56;
//Input:haploid coverage.
pub fn estimate_copy_number_mcmc<R: Rng>(
    nodes: &[(f64, usize)],
    edges: &[Edge],
    cov: f64,
    rng: &mut R,
) -> (Vec<usize>, Vec<usize>) {
    let edges: Vec<_> = edges
        .iter()
        .map(|&(u, u_is_head, v, v_is_head, _)| (u, u_is_head, v, v_is_head))
        .collect();
    let coverages: Vec<_> = nodes
        .iter()
        .map(|&(x, len)| (x.round() as u64, len))
        .collect();
    let graph = crate::copy_number_estimation::copy_number_mrf::Graph::with(&edges, &coverages);
    debug!("COPYNUM\tGraph\t{}", graph);
    use rayon::prelude::*;
    let seeds: Vec<_> = (0..CHAIN_NUM).map(|_| rng.gen::<u64>()).collect();
    let (argmin, min) = seeds
        .into_par_iter()
        .map(|seed| {
            let config = crate::copy_number_estimation::copy_number_mrf::Config::new(cov);
            let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(seed);
            graph.map_estimate_copy_numbers(&mut rng, &config)
        })
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .unwrap();
    trace!("MCMCRES\t{min}");
    argmin
}

pub fn get_potential(
    nodes: &[(f64, usize)],
    node_cp: &[usize],
    edges: &[Edge],
    edge_cp: &[usize],
    cov: f64,
) -> f64 {
    let edges: Vec<_> = edges
        .iter()
        .map(|&(u, u_is_head, v, v_is_head, _)| (u, u_is_head, v, v_is_head))
        .collect();
    let coverages: Vec<_> = nodes
        .iter()
        .map(|&(x, len)| (x.round() as u64, len))
        .collect();
    let graph = crate::copy_number_estimation::copy_number_mrf::Graph::with(&edges, &coverages);
    let config = crate::copy_number_estimation::copy_number_mrf::MCMCConfig::new(1f64, 20f64, cov);
    graph.total_energy(node_cp, edge_cp, &config)
}

pub fn estimate_copy_number_gbs(
    nodes: &[f64],
    edges: &[Edge],
    cov: f64,
) -> (Vec<usize>, Vec<usize>) {
    let graph = GibbsSampler::new(nodes, edges, cov);
    debug!("COPYNUM\tGraph\t{}", graph);
    let config = Config::new(4382094);
    graph.sample_copy_numer(&config)
}

// Optimizer.
#[derive(Debug, Clone)]
struct Optimizer {
    nodes: Vec<(f64, usize)>,
    edges: Vec<Edge>,
    is_tip: Vec<bool>,
    gradient: Vec<f64>,
    momentum: Vec<f64>,
}

impl Optimizer {
    fn new(nodes: &[(f64, usize)], edges: &[Edge], is_tip: &[bool]) -> Self {
        Self {
            nodes: nodes.to_vec(),
            edges: edges.to_vec(),
            is_tip: is_tip.to_vec(),
            gradient: vec![0f64; edges.len()],
            momentum: vec![0f64; edges.len()],
        }
    }
    fn optimize(&mut self, copy_num: &mut [f64]) {
        //debug!("LOSS\tEdgeRes\tNodeRes\tNodePenalty\tIntPenalty\tTotal");
        use std::collections::VecDeque;
        let heat = (1..=100).map(|x| (1f64, 0.01 * x as f64, 0.01 * x as f64));
        let chill = (1..=100).map(|x| (1f64 - 0.01 * x as f64, 1f64, 1f64));
        let re_heat = (1..=100).map(|x| (0.01 * x as f64, 1f64, 1f64));
        let re_chill = (1..=100).map(|x| (1f64 - 0.01 * x as f64, 1f64, 1f64 - 0.01 * x as f64));
        let finalize = (1..=50).map(|x| (0.1, 1f64 - 0.01 * x as f64, 0.01 * x as f64));
        let schedule = heat
            .chain(chill)
            .chain(re_heat)
            .chain(re_chill)
            .chain(finalize);
        // debug!("PARAMS\tEPOCH\tRESIDUAL\tCONSISTENCY\tINTEGER");
        for (epoch, (alpha, beta, gamma)) in schedule.enumerate() {
            //debug!("PARAMS\t{}\t{}\t{}\t{}", epoch, alpha, beta, gamma);
            self.gradient.iter_mut().for_each(|x| *x = 0f64);
            self.momentum.iter_mut().for_each(|x| *x = 0f64);
            copy_num
                .iter_mut()
                .for_each(|x| *x += (epoch as f64 + 1f64).recip());
            let mut loss_log: VecDeque<_> =
                std::iter::repeat(std::f64::INFINITY).take(100).collect();
            for t in 1.. {
                let total_loss = self.update(copy_num, alpha, beta, gamma, t);
                let old = loss_log.pop_front().unwrap();
                loss_log.push_back(total_loss);
                if old - total_loss < 0.0001 {
                    break;
                };
            }
        }
    }
    // Momentum method
    fn update(&mut self, copy_num: &mut [f64], alpha: f64, beta: f64, gamma: f64, t: usize) -> f64 {
        let total_loss = self.update_gradient(copy_num, alpha, beta, gamma);
        let learn_rate = 0.02 / t as f64;
        let moment_coef = 0.9;
        // let (learn_rate, moment_coef) = (0.01, 0f64);
        for ((x, d), moment) in copy_num
            .iter_mut()
            .zip(self.gradient.iter())
            .zip(self.momentum.iter_mut())
        {
            let diff = moment_coef * *moment - learn_rate * d;
            *x += diff;
            *moment = diff;
        }
        total_loss
    }
    fn update_gradient(&mut self, copy_num: &[f64], alpha: f64, beta: f64, gamma: f64) -> f64 {
        // Node residual for copy numbers
        let mut node_residual = vec![0f64; self.nodes.len()];
        // Consisntency penalty.
        let mut node_penalty = vec![0f64; self.nodes.len()];
        for (&(from, fplus, to, tplus, _), cp) in self.edges.iter().zip(copy_num.iter()) {
            let coef = if fplus { -1f64 } else { 1f64 };
            node_penalty[from] += cp * coef;
            let coef = if tplus { -1f64 } else { 1f64 };
            node_penalty[to] += cp * coef;
            node_residual[from] += cp / 2f64;
            node_residual[to] += cp / 2f64;
        }
        for (idx, _) in self.is_tip.iter().enumerate().filter(|(_, &x)| x) {
            node_residual[idx] *= 2f64;
            node_penalty[idx] = 0f64;
        }
        node_residual
            .iter_mut()
            .zip(self.nodes.iter())
            .for_each(|(x, y)| *x -= y.0);
        let integer_penalty: Vec<_> = copy_num
            .iter()
            .map(|x| x - x.round())
            .map(|x| if x.abs() < 0.5 - 0.00001 { x } else { 0f64 })
            .collect();
        for ((g, edge), cp) in self
            .gradient
            .iter_mut()
            .zip(self.edges.iter())
            .zip(copy_num.iter())
        {
            *g = cp - edge.4;
        }
        let edge_res: f64 = self.gradient.iter().map(|x| x * x).sum();
        let node_res: f64 = node_residual
            .iter()
            .zip(self.nodes.iter())
            //.map(|(x, &(_, _))| x * x)
            .map(|(x, &(_, len))| len as f64 * x * x)
            .sum();
        let node_pen: f64 = node_penalty.iter().map(|x| x * x).sum();
        let int_pen: f64 = integer_penalty.iter().map(|x| x * x).sum();
        self.gradient.iter_mut().for_each(|g| *g *= 2f64);
        for ((grad, &(from, fplus, to, tplus, _)), int_pen) in self
            .gradient
            .iter_mut()
            .zip(self.edges.iter())
            .zip(integer_penalty.iter())
        {
            // Grad is currently the residual of the edge.
            let from_grad = node_residual[from] * self.nodes[from].0;
            let to_grad = node_residual[to] * self.nodes[to].0;
            // let from_grad = node_residual[from];
            // let to_grad = node_residual[to];
            let residual = from_grad + to_grad + *grad;
            let from_coef = if fplus { -1f64 } else { 1f64 };
            let to_coef = if tplus { -1f64 } else { 1f64 };
            let node_consist = 2f64 * (from_coef * node_penalty[from] + to_coef * node_penalty[to]);
            *grad = alpha * residual + beta * node_consist + gamma * int_pen;
        }
        alpha * (edge_res + node_res) + beta * node_pen + gamma * int_pen
    }
}

/// Estimate copy number of GFA file.
/// Each segment record should have `cv:i:` samtag and each edge segment record
/// should have `cv:i:` tag and `ln:i:` tag.
/// The `cov` parameter is the haplotype coverage,
/// the `len` parameter is the average length of the raw reads,
/// and `chunk_len` parameter is the length of the chunk.
/// If the assembly graph is gapless, `len` and `chunk_len` would be 0.
/// After estimation, the estimated copy number would be added to gfa as `cp:i:` tag.
pub fn estimate_copy_number_on_gfa(gfa: &mut gfa::GFA, cov: f64, lens: &[usize], chunk_len: usize) {
    let node_index: HashMap<_, _> = gfa
        .iter()
        .filter_map(|record| match &record.content {
            gfa::Content::Seg(node) => Some(&node.sid),
            _ => None,
        })
        .enumerate()
        .map(|(idx, id)| (id.clone(), idx))
        .collect();
    let calibrator = super::CoverageCalibrator::new(lens);
    let old_cov = cov;
    let cov = calibrator.calib_f64(cov, chunk_len);
    debug!("Coverage\t{}\t{}", old_cov, cov);
    let nodes: Vec<_> = gfa
        .iter()
        .filter_map(|record| {
            if let gfa::Content::Seg(seg) = &record.content {
                let coverage: usize = record
                    .tags
                    .iter()
                    .find(|tag| tag.inner.starts_with("cv"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|x| x.parse().ok())?;
                let weight = calibrator.calib(coverage as u64, chunk_len) / cov;
                let len = seg.slen as usize / chunk_len;
                Some((weight, len))
            } else {
                None
            }
        })
        .collect();
    assert_eq!(nodes.len(), node_index.len());
    let edges: Vec<_> = gfa
        .iter()
        .filter_map(|record| match &record.content {
            gfa::Content::Edge(edge) => {
                let from = node_index[&edge.sid1.id];
                let from_plus = edge.beg1.pos == 0;
                let to = node_index[&edge.sid2.id];
                let to_plus = edge.beg2.pos == 0;
                let tags = &record.tags;
                let coverage: usize = tags
                    .iter()
                    .find(|x| x.inner.starts_with("cv"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|cov| cov.parse().ok())
                    .unwrap_or_else(|| panic!("{:?}", record.tags));
                let gap_len: isize = tags
                    .iter()
                    .find(|x| x.inner.starts_with("ln"))
                    .and_then(|tag| tag.inner.split(':').nth(2))
                    .and_then(|cov| cov.parse().ok())
                    .unwrap_or_else(|| panic!("{:?}", record.tags));
                let gap_len = (gap_len + 2 * chunk_len as isize).max(0) as usize;
                let weight = calibrator.calib(coverage as u64, gap_len) / cov;
                Some((from, from_plus, to, to_plus, weight))
            }
            _ => None,
        })
        .collect();
    let (node_cp, edge_cp) = estimate_copy_number(&nodes, &edges);
    assert_eq!(nodes.len(), node_cp.len());
    let nodes = gfa
        .iter_mut()
        .filter(|record| matches!(&record.content, gfa::Content::Seg(_)));
    for (record, cp) in nodes.zip(node_cp) {
        record.tags.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
    }
    assert_eq!(edges.len(), edge_cp.len());
    let edges = gfa
        .iter_mut()
        .filter(|record| matches!(&record.content, gfa::Content::Edge(_)));
    for (record, cp) in edges.zip(edge_cp) {
        record.tags.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
    }
}

/// Estimate copy number.
/// Return copy number of nodes, and copy number of edges.
pub fn estimate_copy_number(nodes: &[(f64, usize)], edges: &[Edge]) -> (Vec<usize>, Vec<usize>) {
    let (is_tip, is_isolated): (Vec<_>, Vec<_>) = {
        // 2 * node_index + is_head.
        let mut degree = vec![0; nodes.len() * 2];
        for &(from, f_plus, to, t_plus, _) in edges.iter() {
            degree[2 * from + f_plus as usize] += 1;
            degree[2 * to + t_plus as usize] += 1;
        }
        // If either position is 0, it is a tip.
        degree
            .chunks_exact(2)
            .map(|w| (w[0] == 0 || w[1] == 0, w[0] == 0 && w[1] == 0))
            .unzip()
    };
    {
        let tip = is_tip.iter().filter(|&&x| x).count();
        let isolated = is_isolated.iter().filter(|&&x| x).count();
        let all = nodes.len();
        debug!("ISOLATED\tTIP\tNORMAL");
        debug!("{}\t{}\t{}", isolated, tip - isolated, all - tip);
    }
    let mut copy_num: Vec<_> = edges.iter().map(|f| f.4).collect();
    let mut optimizer = Optimizer::new(nodes, edges, &is_tip);
    optimizer.optimize(&mut copy_num);
    let mut node_cp = vec![0f64; nodes.len()];
    for (&(from, _, to, _, _), cp) in edges.iter().zip(copy_num.iter()) {
        node_cp[from] += cp / 2f64;
        node_cp[to] += cp / 2f64;
    }
    for (i, _) in is_tip.iter().enumerate().filter(|(_, &x)| x) {
        node_cp[i] *= 2f64;
    }
    for (i, _) in is_isolated.iter().enumerate().filter(|(_, &x)| x) {
        node_cp[i] = nodes[i].0;
    }
    let copy_num: Vec<_> = copy_num.iter().map(|&x| x.round() as usize).collect();
    let node_cp: Vec<_> = node_cp.iter().map(|&x| x.round() as usize).collect();
    (node_cp, copy_num)
}

#[cfg(test)]
mod cpe_test {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
    fn dataset1() -> (Vec<(f64, usize)>, Vec<Edge>) {
        let nodes = vec![(2f64, 1), (1f64, 1), (1f64, 1), (2f64, 1)];
        let edges = vec![
            (0, false, 1, true, 1f64),
            (0, false, 2, true, 1f64),
            (1, false, 3, true, 1f64),
            (2, false, 3, true, 1f64),
        ];
        (nodes, edges)
    }
    #[test]
    fn dataset1_test() {
        let (nodes, edges) = dataset1();
        let (nodes_cp, edges_cp) = estimate_copy_number(&nodes, &edges);
        assert_eq!(nodes_cp, vec![2, 1, 1, 2]);
        assert_eq!(edges_cp, vec![1, 1, 1, 1]);
    }
    type DataSet = (Vec<(f64, usize)>, Vec<usize>, Vec<Edge>, Vec<usize>);
    fn dataset2() -> DataSet {
        let nodes = vec![(1f64, 1), (1f64, 1), (1f64, 1), (1f64, 1)];
        let nodes_answer = vec![1; 4];
        let edges = vec![
            (0, false, 1, true, 1f64),
            (1, false, 2, true, 1f64),
            (2, false, 3, true, 1f64),
            (3, false, 0, true, 1f64),
        ];
        let edges_answer = vec![1; 4];
        (nodes, nodes_answer, edges, edges_answer)
    }
    #[test]
    fn dataset2_test() {
        let (n, na, e, ea) = dataset2();
        let (nc, ec) = estimate_copy_number(&n, &e);
        assert_eq!(nc, na);
        assert_eq!(ec, ea);
    }
    fn dataset3() -> DataSet {
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(4324);
        let nodes_answer = vec![2, 1, 2, 1, 1, 2, 4, 2];
        let nodes: Vec<_> = nodes_answer
            .iter()
            .map(|&cp| {
                let diff = rng.gen_range(-0.3..0.3);
                (cp as f64 + diff, rng.gen_range(1..3))
            })
            .collect();
        let edge_answer = vec![1, 1, 1, 1, 1, 1, 1, 2, 2, 2];
        let edges: Vec<_> = vec![
            (0, false, 1, true, 1),
            (0, false, 2, true, 1),
            (1, false, 2, true, 1),
            (2, false, 3, true, 1),
            (2, false, 4, true, 1),
            (3, false, 5, true, 1),
            (4, false, 5, true, 1),
            (5, false, 6, true, 2),
            (6, false, 6, true, 2),
            (6, false, 7, true, 2),
        ]
        .into_iter()
        .map(|(from, f, to, t, cov)| {
            let diff = rng.gen_range(-0.3..0.3);
            (from, f, to, t, cov as f64 + diff)
        })
        .collect();
        (nodes, nodes_answer, edges, edge_answer)
    }
    #[test]
    fn dataset3_test() {
        let (n, na, e, ea) = dataset3();
        let (nc, ec) = estimate_copy_number(&n, &e);
        assert_eq!(nc, na);
        assert_eq!(ec, ea);
    }
    #[test]
    fn poisson_test() {
        let test_p = poisson(10, 1, 10f64);
        assert!((test_p - 0.12511).abs() < 0.0001, "{}", test_p);
        let test_p = poisson(1, 2, 0.3f64);
        assert!((test_p - 0.329287).abs() < 0.0001, "{}", test_p);
        let test_p = poisson(5, 0, 14f64);
        assert!((test_p - 0.1321686).abs() < 0.0001, "{}", test_p);
    }
    #[test]
    fn simple_case() {
        let nodes: Vec<_> = vec![10f64; 4];
        let edges: Vec<_> = (0..3)
            .map(|from| (from, true, from + 1, false, 10f64))
            .collect();
        let cov = 10f64;
        let graph = GibbsSampler::new(&nodes, &edges, cov);
        let config = Config::default();
        let (node_cp, edge_cp) = graph.sample_copy_numer(&config);
        assert_eq!(node_cp, vec![1; 4]);
        assert_eq!(edge_cp, vec![1; 3]);
    }
    #[test]
    fn long_case_with_rand() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20f64;
        let div = 5f64;
        let loss = 4f64;
        let nodes: Vec<_> = (0..100)
            .map(|_| rng.gen_range(mean_cov - div..mean_cov + div))
            .collect();
        let edges: Vec<_> = (0..99)
            .map(|from| {
                let w = rng.gen_range(mean_cov - div - loss..mean_cov + div - loss);
                (from, true, from + 1, false, w)
            })
            .collect();
        let graph = GibbsSampler::new(&nodes, &edges, mean_cov);
        let config = Config::default();
        let (node_cp, edge_cp) = graph.sample_copy_numer(&config);
        for (i, cp) in node_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("{}\t{}", cp, graph.nodes[i]);
        }
        for (i, cp) in edge_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("{}\t{}", cp, edges[i].4);
        }
        assert_eq!(node_cp, vec![1; 100]);
        assert_eq!(edge_cp, vec![1; 99]);
    }
    #[test]
    fn branching_case() {
        // Two branchings and reducings.
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20f64;
        let div = 5f64;
        let node_cp: Vec<usize> =
            [vec![2; 2], vec![1; 4], vec![2; 3], vec![1; 9], vec![2; 3]].concat();
        let nodes: Vec<f64> = node_cp
            .iter()
            .map(|&copy| {
                let copy = copy as f64;
                rng.gen_range(mean_cov * copy - div..mean_cov * copy + div)
            })
            .collect();
        let edge_cp: Vec<_> = vec![
            (0, true, 1, false, 2),
            (1, true, 2, false, 1),
            (2, true, 3, false, 1),
            (3, true, 6, false, 1),
            (1, true, 4, false, 1),
            (4, true, 5, false, 1), // 5
            (5, true, 6, false, 1),
            (6, true, 7, false, 2),
            (7, true, 8, false, 2),
            (8, true, 9, false, 1),
            (9, true, 10, false, 1), //10
            (10, true, 11, false, 1),
            (11, true, 12, false, 1),
            (12, true, 13, false, 1),
            (8, true, 14, false, 1),
            (14, true, 15, false, 1), //15
            (15, true, 16, false, 1),
            (16, true, 17, false, 1),
            (17, true, 18, false, 2),
            (18, true, 19, false, 2),
            (19, true, 20, false, 2), //20
            (13, true, 18, false, 1), //21
        ];
        let edges: Vec<_> = edge_cp
            .iter()
            .map(|&(f, fd, t, td, copy)| {
                let copy = copy as f64;
                let w = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div);
                (f, fd, t, td, w)
            })
            .collect();
        let graph = GibbsSampler::new(&nodes, &edges, mean_cov);
        let config = Config::default();
        let (node_cp_e, edge_cp_e) = graph.sample_copy_numer(&config);
        let edge_cp: Vec<_> = edge_cp.iter().map(|x| x.4).collect();
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
    #[test]
    fn looping_case() {
        let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20f64;
        let div = 5f64;
        let nodes: Vec<_> = node_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy as f64 - div..mean_cov * copy as f64 + div))
            .collect();
        let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 4, 4, 4, 2, 2];
        let edges = vec![
            (0, true, 1, false),
            (1, true, 2, false),
            (2, false, 3, true),
            (4, true, 3, false),
            (2, true, 4, false),
            (5, true, 2, false), //5
            (6, true, 5, false),
            (2, true, 6, false),
            (2, true, 7, false),
            (7, true, 8, false),
        ];
        let edges: Vec<_> = edges
            .iter()
            .zip(edge_cp.iter())
            .map(|(&(f, fd, t, td), &cp)| {
                let copy = cp as f64;
                let w = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div);
                (f, fd, t, td, w)
            })
            .collect();
        println!("{:?}", nodes);
        println!("{:?}", edges);
        let graph = GibbsSampler::new(&nodes, &edges, mean_cov);
        let config = Config::default();
        let (node_cp_e, edge_cp_e) = graph.sample_copy_numer(&config);
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
    #[test]
    fn complex_case() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 1, 1, 1, 2, 2];
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20f64;
        let div = 5f64;
        let edges = vec![
            (0, true, 1, false),
            (1, true, 2, false),
            (2, true, 3, false),
            (4, true, 4, false),
            (3, true, 5, false),
            (1, true, 3, false), //5
            (3, true, 5, false),
            (5, true, 1, false),
            (5, true, 6, false),
            (6, true, 7, false),
        ];
        let nodes: Vec<_> = node_cp
            .iter()
            .map(|&copy| {
                let copy = copy as f64;
                rng.gen_range(mean_cov * copy - div..mean_cov * copy + div)
            })
            .collect();
        let edges: Vec<_> = edges
            .iter()
            .zip(edge_cp.iter())
            .map(|(&(f, fd, t, td), &cp)| {
                let copy = cp as f64;
                let w = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div);
                (f, fd, t, td, w)
            })
            .collect();
        let graph = GibbsSampler::new(&nodes, &edges, mean_cov);
        let config = Config::default();
        let (node_cp_e, edge_cp_e) = graph.sample_copy_numer(&config);
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
}
