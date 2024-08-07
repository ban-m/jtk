use definitions::EncodedRead;
use log::*;
use rand::Rng;
use std::collections::HashMap;
use std::collections::HashSet;

type Node = (u64, u64);
type Edge = ((Node, bool), (Node, bool));
type CopyNumResult = (Vec<(Node, usize)>, Vec<(Edge, usize)>);
pub trait CopyNumberEstimation {
    fn update_copy_numbers(&mut self, config: &Config);
    fn estimate_copy_numbers(&self, config: &Config) -> CopyNumResult;
}

#[derive(Debug, Clone, Default)]
pub struct Config {
    haploid_coverage: Option<f64>,
}

impl Config {
    pub fn new(haploid_coverage: f64) -> Self {
        Self {
            haploid_coverage: Some(haploid_coverage),
        }
    }
}
const ERROR_FRAC: f64 = 0.05;
const BURN_IN: usize = 20_000;
const TARGET: f64 = 20f64;
const CHOICES: [usize; 3] = [0, 1, 2];

#[derive(Debug, Clone)]
pub struct Graph {
    // (u,u_pos, v, v_pos).
    edges: Vec<(usize, bool, usize, bool)>,
    // coverage of nodes, length(number of chunks) in nodes.
    coverages: Vec<(u64, usize)>,
    // for each node, for each position, return the set of edge indices.
    edge_lists: Vec<[Vec<usize>; 2]>,
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "N:{}\tE:{}", self.coverages.len(), self.edges.len(),)?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct MCMCConfig {
    temprature: f64,
    consist_factor: f64,
    coverage: f64,
}
impl MCMCConfig {
    pub fn new(temprature: f64, consist_factor: f64, coverage: f64) -> Self {
        Self {
            temprature,
            consist_factor,
            coverage,
        }
    }
    // Smaller is better.
    pub fn node_potential(&self, w: u64, copy: usize) -> f64 {
        let cov = self.coverage;
        let lambda = (copy as f64 * cov).max(cov * ERROR_FRAC);
        -(w as f64) * lambda.ln() + lambda
    }
}

impl Graph {
    pub fn edges(&self) -> &[(usize, bool, usize, bool)] {
        self.edges.as_slice()
    }
    pub fn nodes(&self) -> &[(u64, usize)] {
        self.coverages.as_slice()
    }
    pub fn with(edges: &[(usize, bool, usize, bool)], coverages: &[(u64, usize)]) -> Self {
        let mut edge_lists = vec![[vec![], vec![]]; coverages.len()];
        for (idx, &(u, u_is_head, v, v_is_head)) in edges.iter().enumerate() {
            edge_lists[u][u_is_head as usize].push(idx);
            edge_lists[v][v_is_head as usize].push(idx);
        }
        Self {
            coverages: coverages.to_vec(),
            edges: edges.to_vec(),
            edge_lists,
        }
    }
    fn serialize_node<T: std::borrow::Borrow<EncodedRead>>(reads: &[T]) -> HashMap<Node, usize> {
        let nodes: HashSet<_> = reads
            .iter()
            .flat_map(|r| r.borrow().nodes.iter().map(|n| (n.chunk, n.cluster)))
            .collect();
        let mut nodes: Vec<_> = nodes.into_iter().collect();
        nodes.sort_unstable();
        nodes.into_iter().enumerate().map(|(i, k)| (k, i)).collect()
    }
    fn serialize_edge<T: std::borrow::Borrow<EncodedRead>>(reads: &[T]) -> HashMap<Edge, usize> {
        let edges: HashSet<_> = reads
            .iter()
            .flat_map(|r| {
                r.borrow()
                    .nodes
                    .windows(2)
                    .map(|w| Self::normalize(&w[0], &w[1]))
            })
            .collect();
        let mut edges: Vec<_> = edges.into_iter().collect();
        edges.sort_unstable();
        edges.into_iter().enumerate().map(|(i, k)| (k, i)).collect()
    }
    fn normalize(u: &definitions::Node, v: &definitions::Node) -> Edge {
        let u_is_head = !u.is_forward;
        let v_is_head = v.is_forward;
        let u = ((u.chunk, u.cluster), u_is_head);
        let v = ((v.chunk, v.cluster), v_is_head);
        (u.min(v), u.max(v))
    }
    pub fn new<T: std::borrow::Borrow<EncodedRead>>(
        reads: &[T],
    ) -> (Self, HashMap<Node, usize>, HashMap<Edge, usize>) {
        let node_to_idx = Self::serialize_node(reads);
        let edge_to_idx = Self::serialize_edge(reads);
        let mut coverages: Vec<u64> = vec![0; node_to_idx.len()];
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            coverages[node_to_idx[&(node.chunk, node.cluster)]] += 1;
        }
        for read in reads.iter().map(|r| r.borrow()) {
            for (i, edge) in read.edges.iter().enumerate() {
                assert_eq!(edge.from, read.nodes[i].chunk);
                assert_eq!(edge.to, read.nodes[i + 1].chunk);
            }
        }
        let mut edges: Vec<_> = edge_to_idx.keys().collect();
        edges.sort();
        let edges: Vec<_> = edges
            .into_iter()
            .map(|&((u, u_is_head), (v, v_is_head))| {
                trace!("EDGE\t{}\t{}\t{}\t{}", u.0, u_is_head, v.0, v_is_head);
                (node_to_idx[&u], u_is_head, node_to_idx[&v], v_is_head)
            })
            .collect();
        let mut edge_lists = vec![[vec![], vec![]]; node_to_idx.len()];
        for (idx, &(u, u_is_head, v, v_is_head)) in edges.iter().enumerate() {
            edge_lists[u][u_is_head as usize].push(idx);
            edge_lists[v][v_is_head as usize].push(idx);
        }
        let coverages: Vec<_> = coverages.iter().map(|&x| (x, 1)).collect();
        let graph = Self {
            coverages,
            edges,
            edge_lists,
        };
        (graph, node_to_idx, edge_to_idx)
    }
    fn estimate_coverage(&self) -> f64 {
        if self.coverages.is_empty() {
            0f64
        } else {
            let mut weights: Vec<_> = self.coverages.clone();
            let position = weights.len() / 2;
            (weights.select_nth_unstable(position).1).0 as f64 / 2f64
        }
    }
    // Rounding p into p.trunc() + 1/0 depending on the p.fract().
    fn round<R: Rng>(rng: &mut R, f: f64) -> usize {
        f.trunc() as usize + rng.gen_bool(f.fract()) as usize
    }
    fn initial_guess<R: Rng>(&self, config: &MCMCConfig, rng: &mut R) -> (Vec<usize>, Vec<usize>) {
        let node_cp: Vec<_> = self
            .coverages
            .iter()
            .map(|&(x, _)| Self::round(rng, x as f64 / config.coverage))
            .collect();
        let mut edge_cov_dist = vec![0f64; self.edges.len()];
        for (&(cov, _), edges) in self.coverages.iter().zip(self.edge_lists.iter()) {
            for &to in edges[0].iter() {
                edge_cov_dist[to] += cov as f64 / edges[0].len() as f64;
            }
            for &to in edges[1].iter() {
                edge_cov_dist[to] += cov as f64 / edges[1].len() as f64;
            }
        }
        let edge_cp: Vec<_> = edge_cov_dist
            .iter()
            .map(|&cov| Self::round(rng, cov / 2f64 / config.coverage))
            .collect();
        (node_cp, edge_cp)
    }
    // Return vector of (MAP-estimated) copy numbers of nodes and those of edges.
    pub fn map_estimate_copy_numbers<R: Rng>(
        &self,
        rng: &mut R,
        config: &Config,
    ) -> ((Vec<usize>, Vec<usize>), f64) {
        let hap_cov = config
            .haploid_coverage
            .unwrap_or_else(|| self.estimate_coverage());
        let mut mcmc_config = MCMCConfig::new(1f64, 1f64, hap_cov);
        let (mut node_cp, mut edge_cp) = self.initial_guess(&mcmc_config, rng);
        // To get MAP estimates, get the lowest potential combination.
        // We gradually increase the consistency factor so that after BURN_IN period,
        // the consistency factor would be TARGET.
        let total_step = 2 * (node_cp.len() + edge_cp.len()) * BURN_IN;
        let grad = (TARGET.ln() / total_step as f64).exp();
        mcmc_config.temprature = 100f64;
        let chill = (100f64.ln() / total_step as f64).exp();
        for _t in 1..total_step {
            let _ = self.update(&mut node_cp, &mut edge_cp, &mcmc_config, rng);
            mcmc_config.consist_factor *= grad;
            mcmc_config.temprature /= chill;
        }
        let mut minpot = self.total_energy(&node_cp, &edge_cp, &mcmc_config);
        let mut argmin = (node_cp.clone(), edge_cp.clone());
        for _t in 0..1000 {
            let (accept, _) = self.update(&mut node_cp, &mut edge_cp, &mcmc_config, rng);
            if accept {
                let pot = self.total_energy(&node_cp, &edge_cp, &mcmc_config);
                // trace!("DUMP\t{}\t{}\t{}", config.seed, _t, pot);
                if pot < minpot {
                    minpot = pot;
                    argmin = (node_cp.clone(), edge_cp.clone());
                }
            }
        }
        (argmin, minpot)
    }
    // return (if the proposal is accepted, the potential difference between now-prev.
    fn update<R: Rng>(
        &self,
        node_cp: &mut [usize],
        edge_cp: &mut [usize],
        config: &MCMCConfig,
        rng: &mut R,
    ) -> (bool, f64) {
        use rand::seq::SliceRandom;
        (0..)
            .find_map(|_| match CHOICES.choose(rng).unwrap() {
                0 => self.update_node(node_cp, edge_cp, config, rng),
                1 => self.update_edge(node_cp, edge_cp, config, rng),
                2 => self.update_neighbor(node_cp, edge_cp, config, rng),
                _ => unreachable!(),
            })
            .unwrap()
    }
    fn update_node<R: Rng>(
        &self,
        node_cp: &mut [usize],
        edge_cp: &[usize],
        config: &MCMCConfig,
        rng: &mut R,
    ) -> Option<(bool, f64)> {
        if node_cp.is_empty() {
            return None;
        }
        let flip = rng.gen_range(0..node_cp.len());
        let to_decrease = rng.gen_bool(0.5);
        if node_cp[flip] == 0 && to_decrease {
            Some((true, 0f64))
        } else {
            let potential_diff = self.energy_diff_node(node_cp, edge_cp, config, flip, to_decrease);
            let accept_prob_ln = -potential_diff / config.temprature;
            let accept = 0f64 <= accept_prob_ln || rng.gen_bool(accept_prob_ln.exp());
            if accept {
                if to_decrease {
                    node_cp[flip] -= 1;
                } else {
                    node_cp[flip] += 1;
                }
            };
            Some((accept, potential_diff))
        }
    }
    fn update_edge<R: Rng>(
        &self,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        config: &MCMCConfig,
        rng: &mut R,
    ) -> Option<(bool, f64)> {
        if edge_cp.is_empty() {
            return None;
        }
        let flip = rng.gen_range(0..edge_cp.len());
        let to_decrease = rng.gen_bool(0.5);
        if edge_cp[flip] == 0 && to_decrease {
            Some((true, 0f64))
        } else {
            let potential_diff = self.energy_diff_edge(node_cp, edge_cp, config, flip, to_decrease);
            let prob_ln = -potential_diff / config.temprature;
            // let prob = ().exp().min(1f64);
            // let accept = rng.gen_bool(prob);
            let accept = 0f64 <= prob_ln || rng.gen_bool(prob_ln.exp());
            if accept {
                if to_decrease {
                    edge_cp[flip] -= 1;
                } else {
                    edge_cp[flip] += 1;
                }
            }
            Some((accept, potential_diff))
        }
    }
    fn update_neighbor<R: Rng>(
        &self,
        node_cp: &mut [usize],
        edge_cp: &mut [usize],
        config: &MCMCConfig,
        rng: &mut R,
    ) -> Option<(bool, f64)> {
        if node_cp.is_empty() {
            return None;
        }
        use rand::seq::IteratorRandom;
        let flip = rng.gen_range(0..node_cp.len());
        let to_decrease = rng.gen_bool(0.5);
        if to_decrease && node_cp[flip] == 0 {
            return Some((true, 0f64));
        }
        let targets = {
            let mut targets = self.edge_lists[flip]
                .iter()
                .map(|eds| eds.iter().choose_stable(rng));
            // false -> 0, true -> 1, this is very important.
            let tail_target = targets.next().unwrap().copied();
            let tail_target = tail_target.filter(|&n| !to_decrease || 0 < edge_cp[n]);
            let head_target = targets.next().unwrap().copied();
            let head_target = head_target.filter(|&n| !to_decrease || 0 < edge_cp[n]);
            // If this is a loop edge, then merge them.
            if tail_target == head_target {
                (flip, tail_target, None)
            } else {
                (flip, tail_target, head_target)
            }
        };
        let potential_diff =
            self.energy_diff_neighbor(node_cp, edge_cp, config, targets, to_decrease);
        let prob_ln = -potential_diff / config.temprature;
        let accept = 0f64 < prob_ln || rng.gen_bool(prob_ln.exp());
        if accept {
            match to_decrease {
                true => node_cp[flip] -= 1,
                false => node_cp[flip] += 1,
            }
            match (targets.1, to_decrease) {
                (Some(tail), true) => {
                    assert!(edge_cp[tail] > 0);
                    edge_cp[tail] -= 1
                }
                (Some(tail), false) => edge_cp[tail] += 1,
                _ => {}
            }
            match (targets.2, to_decrease) {
                (Some(head), true) => {
                    assert!(edge_cp[head] > 0);
                    edge_cp[head] -= 1
                }
                (Some(head), false) => edge_cp[head] += 1,
                _ => {}
            }
        }
        Some((accept, potential_diff))
    }
    fn energy_diff_node(
        &self,
        node_cp: &[usize],
        edge_cp: &[usize],
        config: &MCMCConfig,
        flip: usize,
        to_decrease: bool,
    ) -> f64 {
        // Energy difference by changing copy number of node.
        let (node_cov, len) = self.coverages[flip];
        let old_cp = node_cp[flip];
        let new_cp = if to_decrease { old_cp - 1 } else { old_cp + 1 };
        let node_diff =
            config.node_potential(node_cov, new_cp) - config.node_potential(node_cov, old_cp);
        let node_diff = node_diff * len as f64;
        let edge_consistency: f64 = self.edge_lists[flip]
            .iter()
            .filter(|eds| !eds.is_empty()) // If there's no edge, then there's no penalty.
            .map(|eds| {
                let total_cp: usize = eds.iter().map(|&e| edge_cp[e]).sum();
                abs_diff(new_cp, total_cp).pow(2) as f64 - abs_diff(old_cp, total_cp).pow(2) as f64
            })
            .sum();
        node_diff + edge_consistency * config.consist_factor
    }
    fn energy_diff_edge(
        &self,
        node_cp: &[usize],
        edge_cp: &[usize],
        config: &MCMCConfig,
        flip: usize,
        to_decrease: bool,
    ) -> f64 {
        // Energy difference by consistency factor.
        let (u, u_is_head, v, v_is_head) = self.edges[flip];
        let consis_diff_on_u: f64 = {
            let u_edges = &self.edge_lists[u][u_is_head as usize];
            let old_cp: usize = u_edges.iter().map(|&e| edge_cp[e]).sum();
            let new_cp = if to_decrease { old_cp - 1 } else { old_cp + 1 };
            if !u_edges.is_empty() {
                abs_diff(node_cp[u], new_cp).pow(2) as f64
                    - abs_diff(node_cp[u], old_cp).pow(2) as f64
            } else {
                0f64
            }
        };
        let consis_diff_on_v: f64 = {
            let v_edges = &self.edge_lists[v][v_is_head as usize];
            let old_cp: usize = v_edges.iter().map(|&e| edge_cp[e]).sum();
            let new_cp = if to_decrease { old_cp - 1 } else { old_cp + 1 };
            if !v_edges.is_empty() {
                abs_diff(node_cp[v], new_cp).pow(2) as f64
                    - abs_diff(node_cp[v], old_cp).pow(2) as f64
            } else {
                0f64
            }
        };
        (consis_diff_on_v + consis_diff_on_u) * config.consist_factor
    }
    fn energy_diff_neighbor(
        &self,
        node_cp: &[usize],
        edge_cp: &[usize],
        config: &MCMCConfig,
        (node, tail_edge, head_edge): (usize, Option<usize>, Option<usize>),
        to_decrease: bool,
    ) -> f64 {
        // potential diff by node.
        let (node_cov, len) = self.coverages[node];
        let old_cp = node_cp[node];
        let new_cp = if to_decrease { old_cp - 1 } else { old_cp + 1 };
        // Potential diff by edges.
        // LK diff by edge potential
        let edges = &self.edge_lists[node];
        let terminals = {
            let mut terminals: Vec<((usize, bool), usize)> = Vec::with_capacity(5);
            match tail_edge {
                None if !edges[0].is_empty() => {
                    if !terminals.iter().any(|x| x.0 == (node, false)) {
                        terminals.push(((node, false), 0));
                    }
                }
                Some(edge_idx) => {
                    let (u, u_is_head, v, v_is_head) = self.edges[edge_idx];
                    match terminals.iter_mut().find(|x| x.0 == (u, u_is_head)) {
                        Some(res) => res.1 += 1,
                        None => terminals.push(((u, u_is_head), 1)),
                    }
                    match terminals.iter_mut().find(|x| x.0 == (v, v_is_head)) {
                        Some(res) => res.1 += 1,
                        None => terminals.push(((v, v_is_head), 1)),
                    }
                }
                _ => {}
            }
            match head_edge {
                None if !edges[1].is_empty() => {
                    if !terminals.iter().any(|x| x.0 == (node, true)) {
                        terminals.push(((node, true), 0));
                    }
                }
                Some(edge_idx) => {
                    let (u, u_is_head, v, v_is_head) = self.edges[edge_idx];
                    match terminals.iter_mut().find(|x| x.0 == (u, u_is_head)) {
                        Some(res) => res.1 += 1,
                        None => terminals.push(((u, u_is_head), 1)),
                    }
                    match terminals.iter_mut().find(|x| x.0 == (v, v_is_head)) {
                        Some(res) => res.1 += 1,
                        None => terminals.push(((v, v_is_head), 1)),
                    }
                }
                _ => {}
            }
            terminals
        };
        let edge_consistency_diff: f64 = terminals
            .iter()
            .map(|&((n, is_head), diff)| {
                let node_old = node_cp[n];
                let node_new = match (n == node, to_decrease) {
                    (true, true) => node_old - 1,
                    (true, false) => node_old + 1,
                    (false, _) => node_old,
                };
                let edge_old: usize = self.edge_lists[n][is_head as usize]
                    .iter()
                    .map(|&i| edge_cp[i])
                    .sum();
                let edge_new = match to_decrease {
                    true => edge_old - diff,
                    false => edge_old + diff,
                };
                let new_diff = abs_diff(node_new, edge_new).pow(2);
                let old_diff = abs_diff(node_old, edge_old).pow(2);
                new_diff as f64 - old_diff as f64
            })
            .sum();
        let node_diff =
            config.node_potential(node_cov, new_cp) - config.node_potential(node_cov, old_cp);
        let node_diff = len as f64 * node_diff;
        node_diff + edge_consistency_diff * config.consist_factor
    }

    // Loged version of the total energy.
    // SMALLER IS BETTER
    pub fn total_energy(&self, node_cp: &[usize], edge_cp: &[usize], config: &MCMCConfig) -> f64 {
        let node_potential: f64 = self
            .coverages
            .iter()
            .zip(node_cp.iter())
            .map(|(&(w, len), &c)| len as f64 * config.node_potential(w, c))
            .sum();
        let edge_consistency: usize = self
            .edge_lists
            .iter()
            .zip(node_cp.iter())
            .map(|(eds_both_sides, &cp)| {
                let tail_eds = &eds_both_sides[0];
                let tail_cp: usize = tail_eds.iter().map(|&idx| edge_cp[idx]).sum();
                let tail_potential: usize = match tail_eds.is_empty() {
                    true => 0,
                    false => abs_diff(cp, tail_cp).pow(2),
                };
                let head_eds = &eds_both_sides[1];
                let head_cp: usize = head_eds.iter().map(|&idx| edge_cp[idx]).sum();
                let head_potential: usize = match head_eds.is_empty() {
                    true => 0,
                    false => abs_diff(cp, head_cp).pow(2),
                };
                head_potential + tail_potential
            })
            .sum();
        // let edpen = edge_consistency as f64 * config.consist_factor;
        // trace!("{node_potential}\t{edpen}");
        node_potential + edge_consistency as f64 * config.consist_factor
    }
}

fn abs_diff(x: usize, y: usize) -> usize {
    x.max(y) - x.min(y)
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128StarStar;
    #[test]
    fn it_works() {}
    #[test]
    fn simple_case() {
        let coverages: Vec<_> = vec![(10, 1); 4];
        let edges = vec![
            (0, false, 1, true),
            (1, false, 2, true),
            (2, false, 3, true),
        ];
        let graph = Graph::with(&edges, &coverages);
        let config = Config::new(10f64);
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp, edge_cp), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
        assert_eq!(node_cp, vec![1; 4]);
        assert_eq!(edge_cp, vec![1; 3]);
    }
    #[test]
    fn long_case_with_rand() {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let coverages: Vec<_> = (0..100)
            .map(|_| (rng.gen_range(mean_cov - div..mean_cov + div), 1))
            .collect();
        let edges: Vec<_> = (0..99).map(|from| (from, false, from + 1, true)).collect();
        let graph = Graph::with(&edges, &coverages);
        for eds in graph.edge_lists.iter().skip(1).take(98) {
            assert_eq!(eds[0].len(), 1);
            assert_eq!(eds[1].len(), 1);
        }
        let config = Config::new(mean_cov as f64);

        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp, edge_cp), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
        for (i, cp) in node_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("ND\t{}\t{}", cp, graph.coverages[i].0);
        }
        for (i, cp) in edge_cp.iter().enumerate().filter(|&(_, &cp)| cp != 1) {
            println!("ED\t{}\t{:?}", cp, graph.edges[i]);
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
            [vec![2; 2], vec![1; 4], vec![2; 3], vec![1; 9], vec![2; 3]].concat();
        let coverages: Vec<_> = node_cp
            .iter()
            .map(|&copy| {
                let cov = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64;
                (cov, 1)
            })
            .collect();
        let _edge_cp: Vec<usize> = [vec![2],
            vec![1; 6],
            vec![2; 2],
            vec![1; 10],
            vec![2; 2],
            vec![1]]
        .concat();
        let edges = vec![
            (0, false, 1, true),
            (1, false, 2, true),
            (2, false, 3, true),
            (3, false, 6, true),
            (1, false, 4, true),
            (4, false, 5, true), //5
            (5, false, 6, true),
            (6, false, 7, true),
            (7, false, 8, true),
            (8, false, 9, true),
            (9, false, 10, true), //10
            (10, false, 11, true),
            (11, false, 12, true),
            (12, false, 13, true),
            (8, false, 14, true),
            (14, false, 15, true), //15
            (15, false, 16, true),
            (16, false, 17, true),
            (17, false, 18, true),
            (18, false, 19, true),
            (19, false, 20, true), //20
            (13, false, 18, true), //21
        ];
        let graph = Graph::with(&edges, &coverages);
        let config = Config::new(mean_cov as f64);

        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp_e, _edge_cp_e), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
        assert_eq!(node_cp_e, node_cp);
        // assert_eq!(edge_cp_e, edge_cp);
    }
    // #[test]
    // fn looping_case() {
    //     let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
    //     let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 4, 4, 4, 2, 2];
    //     let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
    //     let mean_cov = 20;
    //     let div = 5;
    //     let coverages: Vec<_> = node_cp
    //         .iter()
    //         .map(|&copy| {
    //             let cov = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64;
    //             (cov, 1)
    //         })
    //         .collect();
    //     let edges = vec![
    //         (0, false, 1, true),
    //         (1, false, 2, true),
    //         (3, false, 2, true),
    //         (4, false, 3, true),
    //         (2, false, 4, true),
    //         (5, false, 2, true),
    //         (6, false, 5, true),
    //         (2, false, 6, true),
    //         (2, false, 7, true),
    //         (7, false, 8, true),
    //     ];
    //     let graph = Graph::with(&edges, &coverages);
    //     let config = Config::new(mean_cov as f64);

    //     let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
    //     let ((node_cp_e, edge_cp_e), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
    //     println!("{}\t{:?}", mean_cov, coverages);
    //     println!("{:?}\n{:?}", node_cp_e, edge_cp_e);
    //     assert_eq!(node_cp_e, node_cp);
    //     assert_eq!(edge_cp_e, edge_cp);
    // }
    #[test]
    fn complex_case() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let edge_cp: Vec<_> = vec![2, 2, 2, 2, 2, 1, 1, 1, 2, 2];
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let mean_cov = 20;
        let div = 5;
        let coverages: Vec<_> = node_cp
            .iter()
            .map(|&copy| {
                let cov = rng.gen_range(mean_cov * copy - div..mean_cov * copy + div) as u64;
                (cov, 1)
            })
            .collect();
        let edges: Vec<_> = vec![
            (0, false, 1, true),
            (1, false, 2, true),
            (2, false, 3, true),
            (3, false, 4, true),
            (4, false, 5, true),
            (1, false, 3, true),
            (3, false, 5, true),
            (5, false, 1, true),
            (5, false, 6, true),
            (6, false, 7, true),
        ];
        let graph = Graph::with(&edges, &coverages);
        let config = Config::new(mean_cov as f64);

        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp_e, edge_cp_e), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
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
                let position = idx * 2_000;
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
    fn from_reads_1() {
        // Generating reads from looping_case
        let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
        let hap: Vec<_> = vec![0, 1, 2, 4, 3, 2, 6, 5, 2, 6, 5, 2, 7, 8];
        let read_num = 2 * 2_000 * 30 * hap.len() / 10_000;
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let reads: Vec<_> = (0..read_num)
            .map(|i| gen_read(i as u64, &mut rng, &hap))
            .collect();
        let total_occs: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let mean_cov = total_occs as f64 / hap.len() as f64 / 2f64;
        let (graph, node_to_idx, _) = Graph::new(&reads);
        println!("{}", graph);
        let config = Config::new(mean_cov);

        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp_e, _edge_cp_e), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
        let mut node_cp_e: Vec<_> = node_to_idx
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
        let total_occs: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let mean_cov = total_occs as f64 / (hap1.len() + hap2.len()) as f64;
        let (graph, node_to_idx, _) = Graph::new(&reads);
        println!("{}", graph);
        let config = Config::new(mean_cov);

        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(3482309);
        let ((node_cp_e, _edge_cp_e), _) = graph.map_estimate_copy_numbers(&mut rng, &config);
        let mut node_cp_e: Vec<_> = node_to_idx
            .iter()
            .map(|(&(u, _), &idx)| (u, node_cp_e[idx]))
            .collect();
        node_cp_e.sort_by_key(|x| x.0);
        let node_cp_e: Vec<_> = node_cp_e.iter().map(|x| x.1).collect();
        assert_eq!(node_cp_e, node_cp);
    }
}
