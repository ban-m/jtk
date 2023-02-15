//! *DISCARDED* 2021/12/06
//! Estimate copy number of chunks based on coverages.
//! # Example
//! ```rust
//! use haplotyper::copy_number_estimation::*;
//! let mut ds = definitions::DataSet::new();
//! let config = Config::default();
//! // Update copy numbers.
//! ds.update_copy_numbers(&config);
//! // Get copy numbers.
//! let (node_copies, edge_copies) = ds.estimate_copy_numbers(&config);
//! ```
// Out of the haploid coverage, how much fraction would be observed in the 0-coverage nodes.
// For example, if the haploid coverage is 20, we assume that there would be 20 * 0.25 = 5 occurence of edge
// even in the zero-copy elements.
const ERROR_FRAC: f64 = 0.25;
// If true, condidence is set to the MAX_CONFIDENCE from the beggining.
const CONST_CONFIDENCE: bool = false;
const MAX_CONFIDENCE: f64 = 0.95;
// sample copy-number estimation in  `BURN_IN` times from confidence=0 to condidence=MAX_CONFIDENCE,
// then keep sampling `BURN_IN` times at confidence=MAX_CONFIDENCE to reach the stationaly distribution.
const BURN_IN: usize = 1_000;
// After burn-in, sample `SAMPLE_LEN` to get max a posterior estimation.
const SAMPLE_LEN: usize = 1_000;
/// The structure to configure the parameters in copy number estimation.
/// Currently, one can configure the seed used in a random variable generater
/// and the haploid coverage. The latter, haploid coverage, can be estimated automatically
/// by using `Config::estimate_coverage(seed)`.
/// Also, `Config::default()` would do the same thing.
#[derive(Debug, Clone)]
pub struct Config {
    // If None, automatically estimate the haploid coverage.
    haploid_coverage: Option<f64>,
    seed: u64,
}

impl Config {
    pub fn estimate_coverage(seed: u64) -> Self {
        Self {
            haploid_coverage: None,
            seed,
        }
    }
    pub fn new(haploid_coverage: f64, seed: u64) -> Self {
        Self {
            haploid_coverage: Some(haploid_coverage),
            seed,
        }
    }
}

impl std::default::Default for Config {
    fn default() -> Self {
        Self {
            haploid_coverage: None,
            seed: 24309,
        }
    }
}

type Node = (u64, u64);
type Edge = (Node, Node);
type CopyNumResult = (Vec<(Node, usize)>, Vec<(Edge, usize)>);
pub trait CopyNumberEstimation {
    fn update_copy_numbers(&mut self, config: &Config);
    fn estimate_copy_numbers(&self, config: &Config) -> CopyNumResult;
}

use definitions::DataSet;
use definitions::EncodedRead;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use std::collections::HashMap;

impl CopyNumberEstimation for DataSet {
    fn update_copy_numbers(&mut self, config: &Config) {
        let (node_cp, _) = self.estimate_copy_numbers(config);
        // Reset copy numbers.
        self.selected_chunks.iter_mut().for_each(|c| c.copy_num = 0);
        // Update copy numbers.
        let mut chunks: HashMap<u64, &mut definitions::Chunk> =
            self.selected_chunks.iter_mut().map(|c| (c.id, c)).collect();
        for ((chunk, _), cp) in node_cp {
            if let Some(chunk) = chunks.get_mut(&chunk) {
                chunk.copy_num += cp;
            }
        }
    }
    fn estimate_copy_numbers(&self, config: &Config) -> (Vec<(Node, usize)>, Vec<(Edge, usize)>) {
        let chunk_len = 2_000;
        let graph = Graph::new(&self.encoded_reads, chunk_len);
        graph.estimate_copy_numbers(config)
    }
}

#[derive(Debug, Clone)]
pub struct Graph {
    // Node->Node Index.
    node_to_idx: HashMap<Node, usize>,
    // Edge to Edge Index.
    edge_to_idx: HashMap<Edge, usize>,
    // weight of the nodes.
    nodes: Vec<u64>,
    // weight of the edges.
    edges: Vec<u64>,
    // edge index i -> (u,v). No direction needed.
    edge_terminals: Vec<(usize, usize)>,
    // node index i -> set of edge index j such that the j-th edge is connected to the "downstream" direction of nodes[i].
    downstream_edges: Vec<Vec<usize>>,
    // node index i -> set of edge index j such that the j-th edge is connected to the "upstream" direction of nodes[i].
    upstream_edges: Vec<Vec<usize>>,
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "N:{}\tE:{}", self.nodes.len(), self.edges.len(),)?;
        let mut dump_ord: Vec<_> = self.node_to_idx.iter().collect();
        dump_ord.sort_by_key(|x| x.0 .0);
        for ((u, c), &idx) in dump_ord {
            writeln!(f, "{}:{}\t{}", u, c, self.nodes[idx])?;
            for (up1, up2) in self.upstream_edges[idx]
                .iter()
                .map(|&e| self.edge_terminals[e])
            {
                assert!(idx == up1 || idx == up2);
            }
            for (down1, down2) in self.downstream_edges[idx]
                .iter()
                .map(|&e| self.edge_terminals[e])
            {
                assert!(idx == down1 || idx == down2);
            }
        }
        for (&((from, _), (to, _)), &idx) in self.edge_to_idx.iter() {
            writeln!(f, "{}\t{}\t{}", from, to, self.edges[idx])?;
        }
        Ok(())
    }
}

impl Graph {
    fn serialize_node<T: std::borrow::Borrow<EncodedRead>>(reads: &[T]) -> HashMap<Node, usize> {
        let mut nodes: HashMap<Node, usize> = HashMap::new();
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            let len = nodes.len();
            nodes.entry((node.chunk, node.cluster)).or_insert(len);
        }
        nodes
    }
    fn serialize_edge<T: std::borrow::Borrow<EncodedRead>>(reads: &[T]) -> HashMap<Edge, usize> {
        let mut edges: HashMap<(Node, Node), usize> = HashMap::new();
        for w in reads.iter().flat_map(|r| r.borrow().nodes.windows(2)) {
            let len = edges.len();
            edges.entry(Self::normalize(&w[0], &w[1])).or_insert(len);
        }
        edges
    }
    fn normalize(u: &definitions::Node, v: &definitions::Node) -> Edge {
        let u = (u.chunk, u.cluster);
        let v = (v.chunk, v.cluster);
        (u.min(v), u.max(v))
    }
    fn get_connections<T: std::borrow::Borrow<EncodedRead>>(
        node_to_idx: &HashMap<Node, usize>,
        edge_to_idx: &HashMap<Edge, usize>,
        reads: &[T],
    ) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let mut downstream = vec![vec![]; node_to_idx.len()];
        let mut upstream = vec![vec![]; node_to_idx.len()];
        for w in reads.iter().flat_map(|r| r.borrow().nodes.windows(2)) {
            let from = node_to_idx[&(w[0].chunk, w[0].cluster)];
            let to = node_to_idx[&(w[1].chunk, w[1].cluster)];
            let edge_idx = edge_to_idx[&Self::normalize(&w[0], &w[1])];
            let slot = match w[0].is_forward {
                true => downstream.get_mut(from).unwrap(),
                false => upstream.get_mut(from).unwrap(),
            };
            // I do not think this would too slow, even though it takes O(D) times for each insertion..
            if !slot.contains(&edge_idx) {
                slot.push(edge_idx);
            }
            let slot = match w[1].is_forward {
                true => upstream.get_mut(to).unwrap(),
                false => downstream.get_mut(to).unwrap(),
            };
            if !slot.contains(&edge_idx) {
                slot.push(edge_idx);
            }
        }
        (downstream, upstream)
    }
    fn get_terminals(
        node_to_idx: &HashMap<Node, usize>,
        edge_to_idx: &HashMap<Edge, usize>,
    ) -> Vec<(usize, usize)> {
        let mut edges: Vec<_> = edge_to_idx
            .iter()
            .map(|((from, to), idx)| (idx, (node_to_idx[from], node_to_idx[to])))
            .collect();
        edges.sort_by_key(|x| x.0);
        edges.iter().map(|x| x.1).collect()
    }
    fn new<T: std::borrow::Borrow<EncodedRead>>(reads: &[T], chunk_len: usize) -> Self {
        let node_to_idx = Self::serialize_node(reads);
        let edge_to_idx = Self::serialize_edge(reads);
        let (downstream_edges, upstream_edges) =
            Self::get_connections(&node_to_idx, &edge_to_idx, reads);
        let edge_terminals = Self::get_terminals(&node_to_idx, &edge_to_idx);
        let mut nodes: Vec<u64> = vec![0; node_to_idx.len()];
        let mut edges_len: Vec<_> = vec![0; edge_to_idx.len()];
        let mut edges: Vec<u64> = vec![0; edge_to_idx.len()];
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            nodes[node_to_idx[&(node.chunk, node.cluster)]] += 1;
        }
        for read in reads.iter().map(|r| r.borrow()) {
            for (i, edge) in read.edges.iter().enumerate() {
                assert_eq!(edge.from, read.nodes[i].chunk);
                assert_eq!(edge.to, read.nodes[i + 1].chunk);
            }
        }
        for read in reads.iter().map(|r| r.borrow()) {
            for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
                let idx = edge_to_idx[&Self::normalize(&w[0], &w[1])];
                edges[idx] += 1;
                edges_len[idx] += edge.offset;
            }
        }
        // Tune the count of the edges.
        let lens: Vec<_> = reads.iter().map(|r| r.borrow().original_length).collect();
        let calibrator = CoverageCalibrator::new(&lens);
        nodes
            .iter_mut()
            .for_each(|x| *x = calibrator.calib(*x, chunk_len));
        edges.iter_mut().zip(edges_len.iter()).for_each(|(x, len)| {
            let gap_len = (*len / *x as i64 + 2 * chunk_len as i64).max(0) as usize;
            *x = calibrator.calib(*x, gap_len);
        });
        Self {
            node_to_idx,
            edge_to_idx,
            nodes,
            edges,
            edge_terminals,
            downstream_edges,
            upstream_edges,
        }
    }
    // Rounding p into p.trunc() + 1/0 depending on the p.fract().
    fn round<R: Rng>(rng: &mut R, f: f64) -> usize {
        f.trunc() as usize + rng.gen_bool(f.fract()) as usize
    }
    fn estimate_coverage(&self) -> f64 {
        if self.nodes.is_empty() {
            0f64
        } else {
            let mut weights: Vec<_> = self.nodes.clone();
            let position = weights.len() / 2;
            (*weights.select_nth_unstable(position).1) as f64 / 2f64
        }
    }
    fn estimate_copy_numbers(&self, config: &Config) -> CopyNumResult {
        let (node_cp, edge_cp) = self.estimate_copy_numbers_inner(config);
        let mut node_cp: Vec<_> = self
            .node_to_idx
            .iter()
            .map(|(&node, &i)| (node, node_cp[i]))
            .collect();
        node_cp.sort_unstable_by_key(|x| x.0);
        let mut edge_cp: Vec<_> = self
            .edge_to_idx
            .iter()
            .map(|(&edge, &i)| (edge, edge_cp[i]))
            .collect();
        edge_cp.sort_unstable_by_key(|x| x.0);
        (node_cp, edge_cp)
    }
    // Return vector of (estimated) copy numbers of nodes and those of edges.
    fn estimate_copy_numbers_inner(&self, config: &Config) -> (Vec<usize>, Vec<usize>) {
        let hap_cov = config
            .haploid_coverage
            .unwrap_or_else(|| self.estimate_coverage());
        debug!("COPYNUMBER\tHapCov\t{:.3}", hap_cov);
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(config.seed);
        let mut node_cp: Vec<_> = self
            .nodes
            .iter()
            .map(|&x| Self::round(&mut rng, x as f64 / hap_cov))
            .collect();
        let mut edge_cp: Vec<_> = self
            .edges
            .iter()
            .map(|&x| Self::round(&mut rng, x as f64 / hap_cov))
            .collect();
        for i in 0..BURN_IN {
            let confidence = match CONST_CONFIDENCE {
                true => MAX_CONFIDENCE,
                false => i as f64 * MAX_CONFIDENCE / BURN_IN as f64,
            };
            self.update_nodes(&mut node_cp, &edge_cp, hap_cov, confidence, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, hap_cov, confidence, &mut rng);
        }
        // spin loop to reach stationaly distribution.
        for _ in 0..BURN_IN {
            self.update_nodes(&mut node_cp, &edge_cp, hap_cov, MAX_CONFIDENCE, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, hap_cov, MAX_CONFIDENCE, &mut rng);
        }
        // Let's calculate the posterior distributions on each nodes/edges.
        let mut node_cp_dist: Vec<_> = node_cp.iter().map(|&c| vec![0; 2 * (c + 1)]).collect();
        let mut edge_cp_dist: Vec<_> = edge_cp.iter().map(|&c| vec![0; 2 * (c + 1)]).collect();
        for _ in 0..SAMPLE_LEN {
            self.update_nodes(&mut node_cp, &edge_cp, hap_cov, MAX_CONFIDENCE, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, hap_cov, MAX_CONFIDENCE, &mut rng);
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
        hap_cov: f64,
        confidence: f64,
        rng: &mut R,
    ) {
        let edges = self.downstream_edges.iter().zip(self.upstream_edges.iter());
        for (((down_edge, up_edge), &w), cp) in edges.zip(self.nodes.iter()).zip(node_cp.iter_mut())
        {
            let down_cp: usize = down_edge.iter().map(|&j| edge_cp[j]).sum();
            let up_cp: usize = up_edge.iter().map(|&j| edge_cp[j]).sum();
            let cps = [down_cp, up_cp];
            let new_cp = choose_copy_num(w, &cps, hap_cov, confidence, rng);
            *cp = new_cp;
        }
    }
    fn update_edges<R: Rng>(
        &self,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        hap_cov: f64,
        confidence: f64,
        rng: &mut R,
    ) {
        let mut ordering: Vec<_> = (0..edge_cp.len()).collect();
        ordering.shuffle(rng);
        for edge in ordering {
            self.update_edge(edge, node_cp, edge_cp, hap_cov, confidence, rng);
        }
    }
    fn update_edge<R: Rng>(
        &self,
        edge: usize,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        hap_cov: f64,
        confidence: f64,
        rng: &mut R,
    ) {
        let (from, to) = self.edge_terminals[edge];
        let from_cp = self.edge_copy_num(edge, from, node_cp, edge_cp);
        let to_cp = self.edge_copy_num(edge, to, node_cp, edge_cp);
        let cps = [from_cp, to_cp];
        let new_asn = choose_copy_num(self.edges[edge], &cps, hap_cov, confidence, rng);
        edge_cp[edge] = new_asn;
    }
    // Return the current estimation of the copy number of the `edge`.
    // It is calculated by node_cp[node] minus total copy number of edges from `node`.
    fn edge_copy_num(
        &self,
        edge: usize,
        node: usize,
        node_cp: &[usize],
        edge_cp: &[usize],
    ) -> usize {
        let node_cp = node_cp[node];
        let adj: usize = if self.downstream_edges[node].contains(&edge) {
            self.downstream_edges[node]
                .iter()
                .map(|&edge| edge_cp[edge])
                .sum()
        } else {
            assert!(self.upstream_edges[node].contains(&edge), "{}", edge);
            self.upstream_edges[node]
                .iter()
                .map(|&edge| edge_cp[edge])
                .sum()
        };
        (node_cp + edge_cp[edge]).saturating_sub(adj)
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
            //let trust_prob = 0.5 + CONFIDENCE / 2f64;
            let trust_prob = 0.5 + confidence / 2f64;
            choises.push((cp, trust_prob * poisson(w, cp, hap_cov)));
            choises.push((cp + 1, (1f64 - trust_prob) * poisson(w, cp + 1, hap_cov)));
        } else {
            //let trust_prob = 3f64.recip() + 2f64 / 3f64 * CONFIDENCE;
            let trust_prob = 3f64.recip() + 2f64 / 3f64 * confidence;
            choises.push((cp, trust_prob * poisson(w, cp, hap_cov)));
            let minus_one = (1f64 - trust_prob) / 2f64 * poisson(w, cp - 1, hap_cov);
            choises.push((cp - 1, minus_one));
            let plus_one = (1f64 - trust_prob) / 2f64 * poisson(w, cp + 1, hap_cov);
            choises.push((cp + 1, plus_one));
        }
    }
    let sum: f64 = choises.iter().map(|x| x.1).sum();
    choises.choose_weighted(rng, |&(_, w)| w / sum).unwrap().0
}

#[derive(Debug, Clone)]
pub struct CoverageCalibrator {
    // Sorted.
    lengths: Vec<usize>,
    // i->sum of the inverse value from i to the end.
    // cum_inverse_sum: Vec<f64>,
    // i->sum of the length of the reads longer than the i-th read.
    // Note that u64 is as big as 2^64 > total bp possibly.
    cum_sum: Vec<usize>,
    // Mean length.
    mean: f64,
}

impl CoverageCalibrator {
    pub fn new(lens: &[usize]) -> Self {
        let mut lengths = lens.to_vec();
        lengths.sort_unstable();
        let (mut cum_sum, sum): (Vec<_>, _) =
            lengths
                .iter()
                .rev()
                .fold((vec![], 0), |(mut sums, cum), &x| {
                    sums.push(cum + x as usize);
                    (sums, cum + x as usize)
                });
        cum_sum.reverse();
        let mean = sum as f64 / lengths.len() as f64;
        Self {
            lengths,
            cum_sum,
            mean,
        }
    }
    /// Calibrate the observed coverage in the gap-len bp region
    /// into the "actual" coverage, regarding that region as a single-bp region.
    /// Note that if the gap_len is longer than the longest reads in the dataset,
    /// it returns zero.
    pub fn calib(&self, observed: u64, gap_len: usize) -> u64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0,
            Err(idx) => idx,
        };
        let factor = self.cum_sum[idx] - gap_len * (self.lengths.len() - idx);
        let calibed = observed as f64 * self.mean / (factor as f64 / self.lengths.len() as f64);
        calibed.round().max(0f64) as u64
    }
    /// Return the probability that a read span `gap_len` gap at a specific position.
    pub fn prob_spanning(&self, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        (self.cum_sum[idx] - gap_len * (self.lengths.len() - idx)) as f64
            / self.lengths.len() as f64
            / self.mean
    }
}

#[cfg(test)]
mod test {
    use super::*;
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
    fn it_works() {}
    #[test]
    fn simple_case() {
        let nodes: Vec<_> = vec![10; 4];
        let edges: Vec<_> = vec![10; 3];
        let edge_terminals: Vec<_> = vec![(0, 1), (1, 2), (2, 3)];
        let node_to_idx: HashMap<_, _> = HashMap::new();
        let edge_to_idx: HashMap<_, _> = HashMap::new();
        let downstream_edges = vec![vec![0], vec![1], vec![2], vec![]];
        let upstream_edges = vec![vec![], vec![0], vec![1], vec![2]];
        let graph = Graph {
            nodes,
            edges,
            node_to_idx,
            edge_to_idx,
            edge_terminals,
            downstream_edges,
            upstream_edges,
        };
        let config = Config::new(10f64, 4329804);
        let (node_cp, edge_cp) = graph.estimate_copy_numbers_inner(&config);
        assert_eq!(node_cp, vec![1; 4]);
        assert_eq!(edge_cp, vec![1; 3]);
    }
    #[test]
    fn long_case_with_rand() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let loss = 4;
        let nodes: Vec<_> = (0..100)
            .map(|_| rng.gen_range(mean_cov - div..mean_cov + div))
            .collect();
        let edges: Vec<_> = (0..99)
            .map(|_| rng.gen_range(mean_cov - div - loss..mean_cov + div - loss))
            .collect();
        let edge_terminals: Vec<_> = (0..99).map(|i| (i, i + 1)).collect();
        let node_to_idx = HashMap::new();
        let edge_to_idx = HashMap::new();
        let downstream_edges: Vec<_> = (0..100)
            .map(|i| if i < 99 { vec![i] } else { vec![] })
            .collect();
        let upstream_edges: Vec<_> = (0..100)
            .map(|i| if 0 < i { vec![i - 1] } else { vec![] })
            .collect();
        let graph = Graph {
            nodes,
            edges,
            node_to_idx,
            edge_to_idx,
            upstream_edges,
            downstream_edges,
            edge_terminals,
        };
        let config = Config::new(mean_cov as f64, 392480);
        let (node_cp, edge_cp) = graph.estimate_copy_numbers_inner(&config);
        for (i, cp) in node_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("{}\t{}", cp, graph.nodes[i]);
        }
        for (i, cp) in edge_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("{}\t{}", cp, graph.edges[i]);
        }
        assert_eq!(node_cp, vec![1; 100]);
        assert_eq!(edge_cp, vec![1; 99]);
    }
    #[test]
    fn branching_case() {
        // Two branchings and reducings.
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let node_cp: Vec<usize> =
            vec![vec![2; 2], vec![1; 4], vec![2; 3], vec![1; 9], vec![2; 3]].concat();
        let nodes: Vec<_> = node_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edge_cp: Vec<usize> = vec![
            vec![2],
            vec![1; 6],
            vec![2; 2],
            vec![1; 10],
            vec![2; 2],
            vec![1],
        ]
        .concat();
        let edges: Vec<_> = edge_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edge_terminals: Vec<_> = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 6),
            (1, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10), // 10
            (10, 11),
            (11, 12),
            (12, 13),
            (8, 14),
            (14, 15), // 15
            (15, 16),
            (16, 17),
            (17, 18),
            (18, 19),
            (19, 20), //20
            (13, 18),
        ];
        assert_eq!(edge_terminals.len(), edges.len());
        let downstream_edges = vec![
            vec![0],
            vec![1, 4],
            vec![2],
            vec![3],
            vec![5],
            vec![6], // 5
            vec![7],
            vec![8],
            vec![9, 14],
            vec![10],
            vec![11], // 10
            vec![12],
            vec![13],
            vec![21],
            vec![15],
            vec![16], // 15
            vec![17],
            vec![18],
            vec![19],
            vec![20],
            vec![], // 20
        ];
        let upstream_edges = vec![
            vec![],
            vec![0],
            vec![1],
            vec![2],
            vec![4],
            vec![5], // 5
            vec![3, 6],
            vec![7],
            vec![8],
            vec![9],
            vec![10], // 10
            vec![11],
            vec![12],
            vec![13],
            vec![14],
            vec![15], // 15
            vec![16],
            vec![17],
            vec![18, 21],
            vec![19],
            vec![20], // 20
        ];
        let graph = Graph {
            nodes,
            edges,
            node_to_idx: HashMap::new(),
            edge_to_idx: HashMap::new(),
            upstream_edges,
            downstream_edges,
            edge_terminals,
        };
        let config = Config::new(mean_cov as f64, 392480);
        let (node_cp_e, edge_cp_e) = graph.estimate_copy_numbers_inner(&config);
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
    #[test]
    fn looping_case() {
        let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
        let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 4, 4, 4, 2, 2];
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let nodes: Vec<_> = node_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edges: Vec<_> = edge_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edge_terminals = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 2),
            (2, 5),
            (5, 6),
            (6, 2),
            (2, 7),
            (7, 8),
        ];
        let downstream_edges = vec![
            vec![0],
            vec![1],
            vec![4, 7, 8],
            vec![2],
            vec![3],
            vec![5],
            vec![6],
            vec![9],
            vec![],
        ];
        let upstream_edges = vec![
            vec![],
            vec![0],
            vec![1, 2, 5],
            vec![3],
            vec![4],
            vec![6],
            vec![7],
            vec![8],
            vec![9],
        ];
        println!("{:?}", nodes);
        println!("{:?}", edges);
        let graph = Graph {
            nodes,
            edges,
            node_to_idx: HashMap::new(),
            edge_to_idx: HashMap::new(),
            upstream_edges,
            downstream_edges,
            edge_terminals,
        };
        let config = Config::new(mean_cov as f64, 392480);
        let (node_cp_e, edge_cp_e) = graph.estimate_copy_numbers_inner(&config);
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
    #[test]
    fn complex_case() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 1, 1, 1, 2, 2];
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let nodes: Vec<_> = node_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edges: Vec<_> = edge_cp
            .iter()
            .map(|&copy| rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64)
            .collect();
        let edge_terminals = vec![
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (1, 3),
            (3, 5),
            (1, 5),
            (5, 6),
            (6, 7),
        ];
        let downstream_edges = vec![
            vec![0],
            vec![1, 5],
            vec![2],
            vec![3, 6],
            vec![4],
            vec![8, 7], //5
            vec![9],
            vec![],
        ];
        let upstream_edges = vec![
            vec![],
            vec![0, 7],
            vec![1],
            vec![2, 5],
            vec![3],
            vec![4, 6], //5
            vec![8],
            vec![9],
        ];
        println!("{:?}", nodes);
        println!("{:?}", edges);
        let graph = Graph {
            nodes,
            edges,
            node_to_idx: HashMap::new(),
            edge_to_idx: HashMap::new(),
            upstream_edges,
            downstream_edges,
            edge_terminals,
        };
        let config = Config::new(mean_cov as f64, 392480);
        let (node_cp_e, edge_cp_e) = graph.estimate_copy_numbers_inner(&config);
        assert_eq!(node_cp_e, node_cp);
        assert_eq!(edge_cp_e, edge_cp);
    }
    fn gen_read<R: Rng>(id: u64, rng: &mut R, hap: &[u64]) -> EncodedRead {
        let original_length = 5_000 + rng.gen_range(0..10_000);
        let start_pos = rng.gen_range(0..hap.len());
        let seq = vec![b'A'; 2_000];
        let nodes: Vec<_> = hap
            .iter()
            .cycle()
            .skip(start_pos)
            .take(original_length / 2_000)
            .enumerate()
            .map(|(idx, &chunk)| {
                let position = idx as usize * 2_000;
                let cigar = vec![definitions::Op::Match(2_000)];
                definitions::Node::new(chunk, true, seq.clone(), cigar, position, 2)
            })
            .collect();
        let edges = nodes
            .windows(2)
            .map(|ns| {
                let (from, to) = match *ns {
                    [ref from, ref to] => (from, to),
                    _ => unreachable!(),
                };
                let end = from.position_from_start + from.query_length();
                let start = to.position_from_start;
                let label = vec![];
                definitions::Edge {
                    from: from.chunk,
                    to: to.chunk,
                    offset: start as i64 - end as i64,
                    label: label.into(),
                }
            })
            .collect();
        let rem = original_length - original_length / 2_000 * 2_000;
        EncodedRead {
            id,
            original_length,
            leading_gap: vec![].into(),
            trailing_gap: vec![b'A'; rem].into(),
            edges,
            nodes,
        }
    }
    #[test]
    fn read_to_graph_test() {
        let nodes: Vec<_> = (0..2)
            .map(|chunk| {
                let position = chunk as usize * 2000;
                let cigar = vec![definitions::Op::Match(2_000)];
                definitions::Node::new(chunk, true, vec![b'A'; 2_000], cigar, position, 2)
            })
            .collect();
        let edges = nodes
            .windows(2)
            .map(|ns| {
                let (from, to) = match *ns {
                    [ref from, ref to] => (from, to),
                    _ => unreachable!(),
                };
                let end = from.position_from_start + from.query_length();
                let start = to.position_from_start;
                let label = vec![];
                definitions::Edge {
                    from: from.chunk,
                    to: to.chunk,
                    offset: start as i64 - end as i64,
                    label: label.into(),
                }
            })
            .collect();
        let read = EncodedRead {
            id: 0,
            original_length: 0,
            leading_gap: vec![].into(),
            trailing_gap: vec![].into(),
            edges,
            nodes,
        };
        let graph = Graph::new(&[read], 2_000);
        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(graph.edges.len(), 1);
        let nodes: Vec<_> = (0..2).map(|x| graph.node_to_idx[&(x, 0)]).collect();
        assert_eq!(graph.edge_terminals[0], (nodes[0], nodes[1]));
        if nodes[0] == 0 {
            assert_eq!(graph.downstream_edges[0], vec![0]);
            assert!(graph.upstream_edges[0].is_empty());
        } else {
            assert!(graph.downstream_edges[0].is_empty());
            assert!(graph.upstream_edges[0].is_empty());
        }
    }
    #[test]
    fn from_reads_1() {
        // Generating reads from looping_case
        let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
        let hap: Vec<_> = vec![0, 1, 2, 4, 3, 2, 6, 5, 2, 6, 5, 2, 7, 8];
        let read_num = 2 * 2_000 * 30 * hap.len() / 10_000;
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let reads: Vec<_> = (0..read_num)
            .map(|i| gen_read(i as u64, &mut rng, &hap))
            .collect();
        let total_chunks: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let mean_cov = total_chunks / hap.len() / 2;
        let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let calibrator = CoverageCalibrator::new(&lens);
        let ajd = calibrator.calib(mean_cov as u64, 2_000);
        println!("{}\t{}\t{}\t{}", total_chunks, hap.len(), mean_cov, ajd);
        let graph = Graph::new(&reads, 2_000);
        println!("{}", graph);
        let config = Config::estimate_coverage(392480);
        let (node_cp_e, _edge_cp_e) = graph.estimate_copy_numbers_inner(&config);
        let mut node_cp_e: Vec<_> = graph
            .node_to_idx
            .iter()
            .map(|(&(u, _), &idx)| (u, node_cp_e[idx]))
            .collect();
        node_cp_e.sort_by_key(|x| x.0);
        let node_cp_e: Vec<_> = node_cp_e.iter().map(|x| x.1).collect();
        assert_eq!(node_cp_e, node_cp);
    }
    #[test]
    fn from_reads_2() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let hap1: Vec<_> = vec![0, 1, 3, 5, 6, 7];
        let hap2: Vec<_> = vec![0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7];
        let prob_1 = (hap1.len() as f64) / ((hap1.len() + hap2.len()) as f64);
        let read_num = 2_000 * 30 * (hap1.len() + hap2.len()) / 10_000;
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(82304);
        let mut hap1c = 0;
        let reads: Vec<_> = (0..read_num)
            .map(|i| {
                if rng.gen_bool(prob_1) {
                    hap1c += 1;
                    gen_read(i as u64, &mut rng, &hap1)
                } else {
                    gen_read(i as u64, &mut rng, &hap2)
                }
            })
            .collect();
        println!("{},{},{}", read_num, hap1c, read_num - hap1c);
        let count: usize = reads
            .iter()
            .map(|r| {
                r.nodes
                    .windows(2)
                    .filter(|w| w[0].chunk == 1 && w[1].chunk == 3)
                    .count()
            })
            .sum();
        println!("(1,3)\t{}", count);
        let count: usize = reads
            .iter()
            .map(|r| {
                r.nodes
                    .windows(2)
                    .filter(|w| w[0].chunk == 1 && w[1].chunk == 2)
                    .count()
            })
            .sum();
        println!("(1,2)\t{}", count);
        let total_chunks: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let mean_cov = total_chunks / (hap1.len() + hap2.len());
        let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let calibrator = CoverageCalibrator::new(&lens);
        let ajd = calibrator.calib(mean_cov as u64, 2_000);
        println!(
            "{}\t{}\t{}\t{}",
            total_chunks,
            hap1.len() + hap2.len(),
            mean_cov,
            ajd
        );
        // Haplotype coverage should be estimated automatically.
        let graph = Graph::new(&reads, 2_000);
        println!("{}", graph);
        let config = Config::estimate_coverage(392480);
        let (node_cp_e, _edge_cp_e) = graph.estimate_copy_numbers_inner(&config);
        let mut node_cp_e: Vec<_> = graph
            .node_to_idx
            .iter()
            .map(|(&(u, _), &idx)| (u, node_cp_e[idx]))
            .collect();
        node_cp_e.sort_by_key(|x| x.0);
        let node_cp_e: Vec<_> = node_cp_e.iter().map(|x| x.1).collect();
        assert_eq!(node_cp_e, node_cp);
    }
}
