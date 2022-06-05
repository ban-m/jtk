// To me: do not forget update lightedge if you update fatedges
// Maybe it's better to explicitly state the original index of the nodes...?
use rand::{prelude::SliceRandom, Rng};
use std::collections::HashMap;
const LARGE_VALUE: f64 = 1000000f64;
#[derive(Debug, Clone)]
pub struct MSTConfig {
    temperature: f64,
}

impl std::default::Default for MSTConfig {
    fn default() -> Self {
        Self { temperature: 1f64 }
    }
}

#[derive(Debug, Clone)]
pub struct Graph {
    hap_coverage: f64,
    edges: Vec<FatEdge>,
    // nodes: HashMap<Node, usize>,
    // rev_index[nodes[node]] = node.
    // rev_index: Vec<Node>,
    graph: Vec<Vec<LightEdge>>,
    self_loops: Vec<FatEdge>,
    one_degree_nodes: Vec<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Processing,
    Processed,
    Undiscovered,
}

#[derive(Debug, Clone)]
struct LightEdge {
    to: usize,
    is_in_mst: bool,
    penalty_diff_by_decrease: f64,
    penalty_diff_by_increase: f64,
}

impl LightEdge {
    fn new(to: usize) -> Self {
        Self {
            to,
            is_in_mst: false,
            penalty_diff_by_increase: 0f64,
            penalty_diff_by_decrease: 0f64,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FatEdge {
    pub from: usize,
    pub to: usize,
    len: usize,
    target: usize,
    pub copy_number: usize,
    is_in_mst: bool,
    penalty_diff: f64,
}

impl FatEdge {
    pub fn new(from: usize, to: usize, target: usize, len: usize) -> Self {
        let (from, to) = (from.min(to), from.max(to));
        Self {
            from,
            to,
            target,
            copy_number: 0,
            is_in_mst: false,
            penalty_diff: 0f64,
            len,
        }
    }
    pub fn key(&self) -> (usize, usize) {
        (self.from, self.to)
    }
    pub fn penalty_diff(&self, to_increase: bool, cov: f64) -> f64 {
        let current_penalty = penalty(self.target, self.copy_number, cov);
        let proposed_penalty = match to_increase {
            true => penalty(self.target, self.copy_number + 1, cov),
            false if self.copy_number == 0 => LARGE_VALUE,
            false => penalty(self.target, self.copy_number - 1, cov),
        };
        (proposed_penalty - current_penalty) * self.len as f64
    }
}

// sq_error / copy_num / copy_num,
// This is the same as the negative log-likelihood of the Norm(mean=copy_num * hap_cov, var = copy_num * hap_cov)
fn penalty(x: usize, copy_num: usize, hap_cov: f64) -> f64 {
    const ZERO_COPY: f64 = 0.15;
    let mean = hap_cov * copy_num as f64;
    let denom = match copy_num {
        0 => ZERO_COPY * hap_cov,
        _ => copy_num as f64 * hap_cov,
    };
    (x as f64 - mean).powi(2) / denom
}

// type Node = ((u64, u64), super::Position);
// type Edge = (Node, Node);
impl Graph {
    pub fn edges(&self) -> &[FatEdge] {
        self.edges.as_slice()
    }
    pub fn self_loops(&self) -> &[FatEdge] {
        self.self_loops.as_slice()
    }
    pub fn new(
        hap_coverage: f64,
        // nodes: HashMap<Node, usize>,
        edges: Vec<FatEdge>,
        self_loops: Vec<FatEdge>,
    ) -> Self {
        let max_idx = edges.iter().map(|x| x.from.max(x.to)).max().unwrap() + 1;
        let mut graph = vec![Vec::with_capacity(2); max_idx];
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        graph.iter().for_each(|edges| assert!(!edges.is_empty()));
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then(|| i))
            .collect();
        // let mut rev_index = vec![((0, 0), super::Position::Head); nodes.len()];
        // for (&key, &idx) in nodes.iter() {
        //     rev_index[idx] = key;
        // }
        // for edge in edges.iter() {
        //     let from = rev_index[edge.from];
        //     let to = rev_index[edge.to];
        //     debug!("{from:?}<->{to:?}\t{}", edge.target);
        // }
        Self {
            hap_coverage,
            // nodes,
            // rev_index,
            edges,
            graph,
            one_degree_nodes,
            self_loops,
        }
    }
    fn update_lightedges(&mut self) {
        for edge in self.edges.iter() {
            let ledge = self.graph[edge.from]
                .iter_mut()
                .find(|e| e.to == edge.to)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff_by_decrease = edge.penalty_diff(false, self.hap_coverage);
            ledge.penalty_diff_by_increase = edge.penalty_diff(true, self.hap_coverage);
            let ledge = self.graph[edge.to]
                .iter_mut()
                .find(|e| e.to == edge.from)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff_by_decrease = edge.penalty_diff(false, self.hap_coverage);
            ledge.penalty_diff_by_increase = edge.penalty_diff(true, self.hap_coverage);
        }
    }
    const LOOPTIMES: usize = 1000;
    // TODO:Which is better, minimum penalty vs mcmc posterior prob?
    // I think this is not proper MCMC (as proposed distribution does not satisfy
    // detailed balanced condition), so I would like to retin the minimum penalty assignment.
    pub fn update_copy_numbers<R: Rng>(&mut self, rng: &mut R, config: &MSTConfig) {
        let mut dfs_stack = vec![];
        let mut dfs_arrived_status = vec![Status::Undiscovered; self.graph.len()];
        let dfs_stack = &mut dfs_stack;
        let dfs_arrived_status = dfs_arrived_status.as_mut_slice();
        let (mut current_min, mut argmin) = (self.penalty(), self.edges.clone());
        let mut count = 0;
        loop {
            self.update_mst(true);
            self.update_lightedges();
            let (optimal_cycle, penalty_diff) =
                self.find_optimal_cycle(dfs_stack, dfs_arrived_status);
            let prob = (-penalty_diff / config.temperature).exp().min(1f64);
            // trace!("PROPOSED\t{prob:.4}\t{optimal_cycle:?}\t{penalty_diff}");
            if rng.gen_bool(prob) {
                self.update_copy_number_by_cycle(optimal_cycle);
                self.update_self_loops();
                let penalty = self.penalty();
                if penalty < current_min {
                    (current_min, argmin) = (penalty, self.edges.clone());
                }
                // debug!("COPYNUM\t{count}\t{current_min:.3}");
                count += 1;
            } else {
                break;
            }
        }
        for _ in 0..Self::LOOPTIMES {
            let to_increase = rng.gen_bool(0.5);
            self.update_mst(to_increase);
            self.update_lightedges();
            if let Some((_, cycle)) = self.sample_cycle(rng, dfs_stack, dfs_arrived_status) {
                // trace!("UPDATE\t{diff:.3}");
                self.update_copy_number_by_cycle(cycle);
            }
            let penalty = self.penalty();
            if penalty < current_min {
                (current_min, argmin) = (penalty, self.edges.clone());
            }
            // debug!("COPYNUM\t{count}\t{current_min:.3}");
            count += 1;
            self.update_self_loops();
        }
        trace!("MST\t{count}");
        self.edges = argmin;
    }
    fn update_self_loops(&mut self) {
        let self_loop_num = self.self_loops.len();
        for i in 0..self_loop_num {
            for to_increase in [true, false] {
                self.tune_self_loop(i, to_increase);
            }
        }
    }
    fn tune_self_loop(&mut self, i: usize, to_increase: bool) {
        let self_loop = self.self_loops.get_mut(i).unwrap();
        let main_edge = self
            .edges
            .iter_mut()
            .find(|edge| edge.key() == self_loop.key())
            .unwrap();
        let main_edge_diff = main_edge.penalty_diff(to_increase, self.hap_coverage);
        let self_loop_diff = self_loop.penalty_diff(to_increase, self.hap_coverage);
        if main_edge_diff + self_loop_diff < 0f64 {
            if to_increase {
                self_loop.copy_number += 1;
                main_edge.copy_number += 1;
            } else {
                self_loop.copy_number -= 1;
                main_edge.copy_number -= 1;
            }
        }
    }
    fn penalty(&self) -> f64 {
        self.edges
            .iter()
            .map(|e| penalty(e.target, e.copy_number, self.hap_coverage) * e.len as f64)
            .sum()
        // self.edges
        //     .iter()
        //     .map(|e| (e.target as f64 - e.copy_number as f64 * self.hap_coverage).powi(2))
        //     .sum()
    }
    // Search minimum-spanning-tree.
    fn update_mst(&mut self, to_increase: bool) {
        self.graph.iter().for_each(|eds| assert!(!eds.is_empty()));
        // Remove tempolary values
        let hap_coverage = self.hap_coverage;
        self.edges.iter_mut().for_each(|edge| {
            edge.is_in_mst = false;
            edge.penalty_diff = edge.penalty_diff(to_increase, hap_coverage);
        });
        // Find minimum spanning tree.
        self.edges
            .sort_by(|e, f| e.penalty_diff.partial_cmp(&f.penalty_diff).unwrap());
        use crate::find_union;
        let mut fu = find_union::FindUnion::new(self.graph.len());
        for edge in self.edges.iter_mut() {
            let from_parent = fu.find(edge.from).unwrap();
            let to_parent = fu.find(edge.to).unwrap();
            if from_parent != to_parent {
                edge.is_in_mst = true;
                fu.unite(edge.from, edge.to);
            }
        }
        let cl_num: usize = (0..self.graph.len())
            .filter(|&i| fu.find(i).unwrap() == i)
            .count();
        let tree_edges = self.edges.iter().filter(|e| e.is_in_mst).count();
        assert_eq!(tree_edges + cl_num, self.graph.len(), "{}", cl_num);
    }
    fn find_optimal_cycle(
        &self,
        stack: &mut Vec<usize>,
        status: &mut [Status],
    ) -> (Vec<usize>, f64) {
        let (mut current_min, mut argmin) = (LARGE_VALUE, vec![]);
        for edge in self.edges.iter().filter(|e| !e.is_in_mst) {
            let cycle = self
                .find_cycle_between(edge.from, edge.to, stack, status)
                .unwrap_or_else(|| {
                    // let from = self.nodes.iter().find(|&(_, &i)| i == edge.from).unwrap();
                    // let to = self.nodes.iter().find(|&(_, &i)| i == edge.to).unwrap();
                    panic!("{:?}\t{:?}", edge.from, edge.to);
                });
            let penalty = self.penalty_of_cycle(&cycle);
            if penalty < current_min {
                current_min = penalty;
                argmin = cycle;
            }
        }
        for (i, &from) in self.one_degree_nodes.iter().enumerate() {
            for &to in self.one_degree_nodes.iter().skip(i + 1) {
                if let Some(cycle) = self.find_cycle_between(from, to, stack, status) {
                    let penalty = self.penalty_of_cycle(&cycle);
                    if penalty < current_min {
                        current_min = penalty;
                        argmin = cycle;
                    }
                }
            }
        }
        (argmin, current_min)
    }
    fn sample_cycle<R: Rng>(
        &self,
        rng: &mut R,
        stack: &mut Vec<usize>,
        status: &mut [Status],
    ) -> Option<(f64, Vec<usize>)> {
        let onedeg = self.one_degree_nodes.len();
        let mut cycles = Vec::with_capacity(self.edges.len() + onedeg * onedeg);
        let in_cycles = self
            .edges
            .iter()
            .filter(|e| !e.is_in_mst)
            .filter_map(|edge| {
                self.find_cycle_between(edge.from, edge.to, stack, status)
                    .map(|cycle| (self.penalty_of_cycle(&cycle), cycle))
            });
        cycles.extend(in_cycles);
        for (i, &from) in self.one_degree_nodes.iter().enumerate() {
            let cs = self.one_degree_nodes.iter().skip(i + 1).filter_map(|&to| {
                self.find_cycle_between(from, to, stack, status)
                    .map(|cycle| (self.penalty_of_cycle(&cycle), cycle))
            });
            cycles.extend(cs);
        }
        const UPPER: f64 = 50f64;
        cycles.iter_mut().for_each(|x| {
            x.0 = (-x.0).min(UPPER).exp();
        });
        cycles
            .choose_weighted(rng, |x| x.0)
            .ok()
            .and_then(|(x, cycle)| rng.gen_bool(x.min(1f64)).then(|| (*x, cycle.clone())))
    }

    fn find_cycle_between(
        &self,
        from: usize,
        to: usize,
        stack: &mut Vec<usize>,
        status: &mut [Status],
    ) -> Option<Vec<usize>> {
        // DFS...
        stack.clear();
        status.iter_mut().for_each(|x| *x = Status::Undiscovered);
        stack.push(from);
        'dfs: while !stack.is_empty() {
            let last = *stack.last().unwrap();
            if last == to {
                // Find target node.
                let mut cycle = Vec::with_capacity(stack.len() + 1);
                cycle.extend(stack.iter());
                cycle.push(from);
                return Some(cycle);
            }
            status[last] = Status::Processing;
            let edge_in_tree = self.graph[last].iter().filter(|e| e.is_in_mst);
            let edge_undiscovered = edge_in_tree.filter(|e| status[e.to] == Status::Undiscovered);
            for take_edge in edge_undiscovered {
                stack.push(take_edge.to);
                continue 'dfs;
            }
            let fin_node = stack.pop().unwrap();
            status[fin_node] = Status::Processed;
        }
        None
    }
    fn penalty_of_cycle(&self, cycle: &[usize]) -> f64 {
        let from_up = self.penalty_of_cycle_from(cycle, true);
        let from_down = self.penalty_of_cycle_from(cycle, false);
        from_up.min(from_down)
    }
    fn update_copy_number_by_cycle(&mut self, cycle: Vec<usize>) {
        let from_up = self.penalty_of_cycle_from(&cycle, true);
        let from_down = self.penalty_of_cycle_from(&cycle, false);
        if from_up < from_down {
            self.update_copy_number_by_cycle_from(&cycle, true);
        } else {
            self.update_copy_number_by_cycle_from(&cycle, false);
        }
    }
    fn penalty_of_cycle_from(&self, cycle: &[usize], mut change_direction: bool) -> f64 {
        let mut score = 0f64;
        let mut is_prev_e_edge = false;
        for (from, to) in cycle.windows(2).map(|w| (w[0], w[1])) {
            let edge = match self.graph[from].iter().find(|e| e.to == to) {
                Some(edge) => edge,
                None => continue,
            };
            let is_e_edge = from / 2 != to / 2;
            if is_prev_e_edge && is_e_edge {
                change_direction = !change_direction;
            }
            score += match change_direction {
                true => edge.penalty_diff_by_increase,
                false => edge.penalty_diff_by_decrease,
            };
            is_prev_e_edge = is_e_edge;
        }
        score
    }
    fn update_copy_number_by_cycle_from(&mut self, cycle: &[usize], mut change_direction: bool) {
        let mut mod_edges: HashMap<_, bool> = HashMap::new();
        let mut is_prev_e_edge = false;
        for (from, to) in cycle.windows(2).map(|w| (w[0], w[1])) {
            let is_e_edge = from / 2 != to / 2;
            if is_prev_e_edge && is_e_edge {
                change_direction = !change_direction;
            }
            mod_edges.insert((from.min(to), from.max(to)), change_direction);
            is_prev_e_edge = is_e_edge;
        }
        self.edges
            .iter_mut()
            .for_each(|e| match mod_edges.get(&(e.from, e.to)) {
                Some(true) => e.copy_number += 1,
                Some(false) => e.copy_number -= 1,
                None => {}
            });
    }
    // pub fn self_loop_copy_numbers(&self) -> Vec<((u64, u64), usize)> {
    //     self.self_loops
    //         .iter()
    //         .map(|edge| (self.rev_index[edge.from].0, edge.copy_number))
    //         .collect()
    // }
    // pub fn node_copy_numbers(&self) -> HashMap<(u64, u64), usize> {
    //     self.edges
    //         .iter()
    //         .filter_map(|edge| {
    //             let from = self.rev_index[edge.from].0;
    //             let to = self.rev_index[edge.to].0;
    //             (from == to).then(|| (from, edge.copy_number))
    //         })
    //         .collect()
    // }
    // pub fn edge_copy_numbers(&self) -> HashMap<Edge, usize> {
    //     let mut edge_copy = HashMap::new();
    //     for edge in self.edges.iter() {
    //         let from = self.rev_index[edge.from];
    //         let to = self.rev_index[edge.to];
    //         let cp = edge.copy_number;
    //         if from.0 != to.0 {
    //             edge_copy.insert((from, to), cp);
    //             edge_copy.insert((to, from), cp);
    //         }
    //     }
    //     edge_copy
    // }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    fn mock_data_1() -> Graph {
        // let nodes = HashMap::new();
        let mut edges: Vec<_> = vec![];
        edges.push(FatEdge::new(0, 1, 2, 1));
        edges.push(FatEdge::new(1, 2, 1, 1));
        edges.push(FatEdge::new(1, 4, 1, 1));
        edges.push(FatEdge::new(2, 3, 1, 1));
        edges.push(FatEdge::new(4, 5, 1, 1));
        edges.push(FatEdge::new(3, 6, 1, 1));
        edges.push(FatEdge::new(5, 6, 1, 1));
        edges.push(FatEdge::new(6, 7, 2, 1));
        edges.push(FatEdge::new(7, 8, 1, 1));
        edges.push(FatEdge::new(7, 10, 1, 1));
        edges.push(FatEdge::new(8, 9, 1, 1));
        edges.push(FatEdge::new(10, 11, 1, 1));
        edges.push(FatEdge::new(9, 12, 1, 1));
        edges.push(FatEdge::new(12, 11, 1, 1));
        edges.push(FatEdge::new(13, 12, 2, 1));
        let mut graph = vec![Vec::with_capacity(2); 14];
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then(|| i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 13]);
        Graph {
            hap_coverage: 1f64,
            edges,
            // rev_index: vec![],
            // nodes,
            graph,
            self_loops: vec![],
            one_degree_nodes,
        }
    }
    fn mock_data_2() -> Graph {
        // let nodes = HashMap::new();
        let mut edges: Vec<_> = vec![];
        edges.push(FatEdge::new(0, 1, 2, 1));
        edges.push(FatEdge::new(1, 2, 2, 1));
        edges.push(FatEdge::new(2, 3, 3, 1));
        edges.push(FatEdge::new(3, 4, 1, 1));
        edges.push(FatEdge::new(3, 6, 2, 1));
        edges.push(FatEdge::new(4, 5, 1, 1));
        edges.push(FatEdge::new(5, 6, 1, 1));
        edges.push(FatEdge::new(6, 7, 3, 1));
        edges.push(FatEdge::new(7, 8, 2, 1));
        edges.push(FatEdge::new(7, 10, 1, 1));
        edges.push(FatEdge::new(8, 9, 2, 1));
        edges.push(FatEdge::new(9, 10, 2, 1));
        edges.push(FatEdge::new(10, 11, 3, 1));
        edges.push(FatEdge::new(11, 2, 1, 1));
        edges.push(FatEdge::new(11, 12, 2, 1));
        edges.push(FatEdge::new(12, 13, 2, 1));
        let mut graph = vec![Vec::with_capacity(2); 14];
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then(|| i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 13]);
        Graph {
            hap_coverage: 1f64,
            // rev_index: vec![],
            edges,
            // nodes,
            graph,
            self_loops: vec![],
            one_degree_nodes,
        }
    }
    fn mock_data_3() -> Graph {
        // let nodes = HashMap::new();
        let mut edges: Vec<_> = vec![];
        edges.push(FatEdge::new(0, 1, 20, 1));
        edges.push(FatEdge::new(1, 2, 5, 1));
        edges.push(FatEdge::new(2, 3, 5, 1));
        edges.push(FatEdge::new(6, 3, 5, 1));
        edges.push(FatEdge::new(1, 4, 5, 1));
        edges.push(FatEdge::new(4, 5, 5, 1));
        edges.push(FatEdge::new(5, 6, 5, 1));
        edges.push(FatEdge::new(6, 7, 20, 1));
        let mut graph = vec![Vec::with_capacity(2); 8];
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then(|| i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 7]);
        Graph {
            hap_coverage: 10f64,
            edges,
            // nodes,
            graph,
            self_loops: vec![],
            one_degree_nodes,
        }
    }

    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128PlusPlus;
    #[test]
    fn mock_data_1_test() {
        let mut g = mock_data_1();
        let mut rng = Xoroshiro128PlusPlus::seed_from_u64(20329);
        let config = MSTConfig::default();
        g.update_copy_numbers(&mut rng, &config);
        let two_copy = g.edges.iter().filter(|e| e.copy_number == 2).count();
        assert_eq!(two_copy, 3);
        for e in g.edges.iter() {
            println!("{e:?}");
            assert!([1, 2].contains(&e.copy_number));
        }
    }
    #[test]
    fn mock_data_2_test() {
        let mut g = mock_data_2();
        let mut rng = Xoroshiro128PlusPlus::seed_from_u64(20329);
        let config = MSTConfig::default();
        g.update_copy_numbers(&mut rng, &config);
        for edge in g.edges.iter() {
            assert_eq!(edge.copy_number, edge.target);
        }
    }
    #[test]
    fn mock_data_3_test() {
        let mut g = mock_data_3();
        let mut rng = Xoroshiro128PlusPlus::seed_from_u64(20329);
        let config = MSTConfig::default();
        g.update_copy_numbers(&mut rng, &config);
        for edge in g.edges.iter() {
            if edge.from == 0 || edge.to == 7 {
                assert_eq!(edge.copy_number, 2);
            } else {
                assert_eq!(edge.copy_number, 1);
            }
        }
    }
}
