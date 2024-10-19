// To me: do not forget update lightedge if you update fatedges
// Maybe it's better to explicitly state the original index of the nodes...?
use rand::{prelude::SliceRandom, Rng};
use rayon::prelude::*;
use std::collections::HashMap;
const LARGE_VALUE: f64 = 1000000f64;
use log::*;

const INIT_NEG_COPY_NUM_PEN: f64 = 100f64;
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
pub struct FitParam {
    hap_coverage: f64,
    negative_copy_num_pen: f64,
}

impl FitParam {
    fn new(hap_coverage: f64, negative_copy_num_pen: f64) -> Self {
        Self {
            hap_coverage,
            negative_copy_num_pen,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Graph {
    hap_coverage: f64,
    edges: Vec<FatEdge>,
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
    pub copy_number: isize,
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
    // Penalty difference. Lower is better.
    pub fn penalty_diff(&self, to_increase: bool, param: &FitParam) -> f64 {
        let cov = param.hap_coverage;
        let neg_pen = param.negative_copy_num_pen;
        let current_penalty = match self.copy_number {
            x if x < 0 => neg_pen * -self.copy_number as f64,
            _ => penalty(self.target, self.copy_number, cov),
        };
        let next_copy_num = match to_increase {
            true => self.copy_number + 1,
            false => self.copy_number - 1,
        };
        let proposed_penalty = match next_copy_num.cmp(&0) {
            std::cmp::Ordering::Less => neg_pen * -next_copy_num as f64,
            _ => penalty(self.target, next_copy_num, cov),
        };
        (proposed_penalty - current_penalty) * self.len as f64
    }
}

// sq_error / copy_num / copy_num,
// This is the same as the negative log-likelihood of the Norm(mean=copy_num * hap_cov, var = copy_num * hap_cov)
fn penalty(x: usize, copy_num: isize, hap_cov: f64) -> f64 {
    assert!(0 <= copy_num);
    const ZERO_COPY: f64 = 0.15;
    let mean = hap_cov * copy_num as f64;
    let denom = match copy_num {
        0 => ZERO_COPY * hap_cov,
        _ => copy_num as f64 * hap_cov,
    };
    (x as f64 - mean).powi(2) / denom
}

type CopyNumber = (Vec<FatEdge>, Vec<FatEdge>);
impl Graph {
    pub fn sanity_check(&self) -> bool {
        let mut diff_cp = vec![0isize; self.graph.len()];
        for edge in self.edges.iter() {
            if edge.from / 2 == edge.to / 2 {
                diff_cp[edge.from] -= edge.copy_number;
                diff_cp[edge.to] -= edge.copy_number;
            } else {
                diff_cp[edge.from] += edge.copy_number;
                diff_cp[edge.to] += edge.copy_number;
            }
        }
        for edge in self.self_loops.iter() {
            diff_cp[edge.from] += edge.copy_number;
            diff_cp[edge.to] += edge.copy_number;
        }
        for (i, &e) in diff_cp.iter().enumerate() {
            if e != 0 && !self.one_degree_nodes.contains(&i) {
                eprintln!("ERR\t{i}\t{e}");
                for edge in self.edges.iter().filter(|e| e.from == i || e.to == i) {
                    eprintln!("\tEDGE\t{edge:?}");
                }
                for edge in self.self_loops.iter().filter(|e| e.from == i || e.to == i) {
                    eprintln!("\tEDGE\t{edge:?}");
                }
            }
        }
        diff_cp
            .iter()
            .enumerate()
            .all(|(i, &e)| e == 0 || self.one_degree_nodes.contains(&i))
    }
    pub fn clone_copy_number(&self) -> CopyNumber {
        (self.edges.clone(), self.self_loops.clone())
    }
    pub fn set_copy_number(&mut self, (edge_cp, self_loop): CopyNumber) {
        assert_eq!(self.edges.len(), edge_cp.len());
        assert_eq!(self.self_loops.len(), self_loop.len());
        self.edges = edge_cp;
        self.self_loops = self_loop;
    }
    pub fn edges(&self) -> &[FatEdge] {
        self.edges.as_slice()
    }
    pub fn self_loops(&self) -> &[FatEdge] {
        self.self_loops.as_slice()
    }
    pub fn new(hap_coverage: f64, edges: Vec<FatEdge>, self_loops: Vec<FatEdge>) -> Self {
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
            .filter_map(|(i, eds)| (eds.len() == 1).then_some(i))
            .collect();
        Self {
            hap_coverage,
            edges,
            graph,
            one_degree_nodes,
            self_loops,
        }
    }
    fn update_lightedges(&mut self, param: &FitParam) {
        for edge in self.edges.iter() {
            let ledge = self.graph[edge.from]
                .iter_mut()
                .find(|e| e.to == edge.to)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff_by_decrease = edge.penalty_diff(false, param);
            ledge.penalty_diff_by_increase = edge.penalty_diff(true, param);
            let ledge = self.graph[edge.to]
                .iter_mut()
                .find(|e| e.to == edge.from)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff_by_decrease = edge.penalty_diff(false, param);
            ledge.penalty_diff_by_increase = edge.penalty_diff(true, param);
        }
    }
    const LOOPTIMES: usize = 500;
    pub fn update_copy_numbers<R: Rng>(&mut self, rng: &mut R, config: &MSTConfig) {
        let (min, argmin) = (0..10)
            .map(|_| self.update_copy_numbers_inner(rng, config))
            .min_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
            .unwrap();
        self.set_copy_number(argmin);
        debug!("MIN\t{min:.1}");
        let neg_copy_num = self.edges.iter().filter(|x| x.copy_number < 0).count();
        debug!("NEGCOPY\t{neg_copy_num}");
    }
    fn update_copy_numbers_inner<R: Rng>(
        &mut self,
        rng: &mut R,
        config: &MSTConfig,
    ) -> (f64, CopyNumber) {
        for edge in self.edges.iter() {
            let count = self.graph[edge.from]
                .iter()
                .filter(|e| e.to == edge.to)
                .count();
            assert_eq!(count, 1);
            let count = self.graph[edge.to]
                .iter()
                .filter(|e| e.to == edge.from)
                .count();
            assert_eq!(count, 1);
        }
        self.edges.iter_mut().for_each(|e| e.copy_number = 0);
        self.self_loops.iter_mut().for_each(|e| e.copy_number = 0);
        let mut parameters = FitParam::new(self.hap_coverage, INIT_NEG_COPY_NUM_PEN);
        let test_param = FitParam::new(self.hap_coverage, LARGE_VALUE);
        let (mut current_min, mut argmin) = (self.penalty(&test_param), self.clone_copy_number());
        loop {
            self.update_mst(true, &parameters);
            parameters.negative_copy_num_pen =
                (parameters.negative_copy_num_pen * 1.05f64).min(LARGE_VALUE);
            self.update_lightedges(&parameters);
            let (optimal_cycle, penalty_diff) = self.find_optimal_cycle();
            // -0.1 to make sure that even if the penalty is the same, the have non-zero prob to break.
            let penalty_diff = penalty_diff + 0.01;
            let prob = (-penalty_diff / config.temperature).exp().min(1f64);
            let len = optimal_cycle.len();
            let pen = parameters.negative_copy_num_pen;
            let (cmin, neg) = (current_min, self.count_neg());
            trace!("PROPOSED\t{prob:.1}\t{penalty_diff:.1}\t{len}\t{pen:.1}\t{cmin:.1}\t{neg}");
            if rng.gen_bool(prob) {
                self.update_copy_number_by_cycle(optimal_cycle);
                self.update_self_loops(&parameters);
                let penalty = self.penalty(&test_param);
                if penalty < current_min {
                    (current_min, argmin) = (penalty, self.clone_copy_number());
                }
            } else {
                break;
            }
        }
        trace!("RANDOM Mode");
        for _t in 0..Self::LOOPTIMES {
            parameters.negative_copy_num_pen =
                (parameters.negative_copy_num_pen * 1.05).min(LARGE_VALUE);
            let to_increase = rng.gen_bool(0.5);
            // self.update_mst(to_increase, &parameters);
            self.update_mst_random(to_increase, rng, &parameters);
            self.update_lightedges(&parameters);
            if let Some((_diff, cycle)) = self.sample_cycle(rng) {
                let len = cycle.len();
                let pen = parameters.negative_copy_num_pen;
                let (cmin, neg) = (current_min, self.count_neg());
                trace!("PROPOSED\t1.0\t{_diff:.1}\t{len}\t{pen:.1}\t{cmin:.1}\t{neg}");
                self.update_copy_number_by_cycle(cycle);
            }
            self.update_self_loops(&parameters);
            let penalty = self.penalty(&test_param);
            if penalty < current_min {
                (current_min, argmin) = (penalty, self.clone_copy_number());
            }
        }
        (current_min, argmin)
    }
    fn count_neg(&self) -> usize {
        self.edges.iter().filter(|e| e.copy_number < 0).count()
    }
    fn update_self_loops(&mut self, param: &FitParam) {
        let self_loop_num = self.self_loops.len();
        for i in 0..self_loop_num {
            self.tune_self_loop(i, true, param);
            self.tune_self_loop(i, false, param);
        }
    }
    fn tune_self_loop(&mut self, i: usize, to_increase: bool, param: &FitParam) {
        let self_loop = self.self_loops.get_mut(i).unwrap();
        let num_main_edge = self
            .edges
            .iter()
            .filter(|e| e.key() == self_loop.key())
            .count();
        assert_eq!(num_main_edge, 1);
        let main_edge = self
            .edges
            .iter_mut()
            .find(|edge| edge.key() == self_loop.key())
            .unwrap();
        let self_loop_diff = self_loop.penalty_diff(to_increase, param);
        let main_edge_diff = main_edge.penalty_diff(to_increase, param);
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
    fn penalty(&self, param: &FitParam) -> f64 {
        self.edges
            .iter()
            .map(|e| match e.copy_number {
                x if x < 0 => -e.copy_number as f64 * param.negative_copy_num_pen,
                _ => penalty(e.target, e.copy_number, param.hap_coverage) * e.len as f64,
            })
            .sum()
    }
    // Search minimum-spanning-tree.
    fn update_mst(&mut self, to_increase: bool, param: &FitParam) {
        self.graph.iter().for_each(|eds| assert!(!eds.is_empty()));
        // Remove tempolary values
        self.edges.iter_mut().for_each(|edge| {
            edge.is_in_mst = false;
            edge.penalty_diff = edge.penalty_diff(to_increase, param);
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
    // Randomly construct spanning tree.
    fn update_mst_random<R: Rng>(&mut self, to_increase: bool, rng: &mut R, param: &FitParam) {
        self.graph.iter().for_each(|eds| assert!(!eds.is_empty()));
        // Remove tempolary values
        self.edges.iter_mut().for_each(|edge| {
            edge.is_in_mst = false;
            edge.penalty_diff = edge.penalty_diff(to_increase, param);
        });
        // Find random spanning tree.
        use crate::find_union;
        let mut fu = find_union::FindUnion::new(self.graph.len());
        let weight = |idx: usize| 1f64 - self.edges[idx].penalty_diff.min(0f64);
        let edge_num = self.edges.len();
        let sampled = rand::seq::index::sample_weighted(rng, edge_num, weight, edge_num).unwrap();
        for edge_idx in sampled {
            let edge = self.edges.get_mut(edge_idx).unwrap();
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
    fn find_optimal_cycle(&self) -> (Vec<usize>, f64) {
        let (mut current_min, mut argmin) = self
            .edges
            .par_iter()
            .filter(|e| !e.is_in_mst)
            .map(|edge| {
                let mut stack = vec![];
                let mut status = vec![Status::Undiscovered; self.graph.len()];
                self.find_cycle_between(edge.from, edge.to, &mut stack, &mut status)
                    .map(|cycle| (self.penalty_of_cycle(&cycle), cycle))
                    .unwrap_or_else(|| panic!("{:?}\t{:?}", edge.from, edge.to))
            })
            .min_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
            .unwrap_or((LARGE_VALUE, vec![]));
        let mut stack = vec![];
        let mut status = vec![Status::Undiscovered; self.graph.len()];
        // If there are too many one-degree nodes, the inference gets slower and slower.
        // (from,to) should be resides in the same connected components....
        // Use find union tree to check if two node is in the same ....
        // No, we should first enumerate reachable one-degree nodes at first.
        // But we usually see no problem here...
        for (i, &from) in self.one_degree_nodes.iter().enumerate() {
            for &to in self.one_degree_nodes.iter().skip(i + 1) {
                if let Some(cycle) = self.find_cycle_between(from, to, &mut stack, &mut status) {
                    let penalty = self.penalty_of_cycle(&cycle);
                    trace!("CYCLE\t{penalty:.2}");
                    for idx in cycle.iter() {
                        let (node, tip) = (idx / 2, idx % 2);
                        trace!("CYCLE\t{node}\t{tip}");
                    }
                    if penalty < current_min {
                        current_min = penalty;
                        argmin = cycle;
                    }
                }
            }
        }
        (argmin, current_min)
    }
    fn sample_cycle<R: Rng>(&self, rng: &mut R) -> Option<(f64, Vec<usize>)> {
        let onedeg = self.one_degree_nodes.len();
        let mut cycles = Vec::with_capacity(self.edges.len() + onedeg * onedeg);
        let in_cycles = self
            .edges
            .par_iter()
            .filter(|e| !e.is_in_mst)
            .filter_map(|edge| {
                let mut stack = vec![];
                let mut status = vec![Status::Undiscovered; self.graph.len()];
                self.find_cycle_between(edge.from, edge.to, &mut stack, &mut status)
                    .map(|cycle| (self.penalty_of_cycle(&cycle), cycle))
            });
        cycles.par_extend(in_cycles);
        let mut stack = vec![];
        let mut status = vec![Status::Undiscovered; self.graph.len()];
        for (i, &from) in self.one_degree_nodes.iter().enumerate() {
            let cs = self.one_degree_nodes.iter().skip(i + 1).filter_map(|&to| {
                self.find_cycle_between(from, to, &mut stack, &mut status)
                    .map(|cycle| (self.penalty_of_cycle(&cycle), cycle))
            });
            cycles.extend(cs);
        }
        if cycles.iter().any(|x| x.0 < -0.01) {
            let min = cycles.iter().min_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
            return min.cloned();
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
            let mut edge_undiscovered =
                edge_in_tree.filter(|e| status[e.to] == Status::Undiscovered);
            if let Some(take_edge) = edge_undiscovered.next() {
                stack.push(take_edge.to);
                continue 'dfs;
            }
            // for take_edge in edge_undiscovered {
            //     stack.push(take_edge.to);
            //     continue 'dfs;
            // }
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
        if !self.sanity_check() {
            panic!("AF\t{:?}", cycle);
        }
    }
    fn penalty_of_cycle_from(&self, cycle: &[usize], start_direction: bool) -> f64 {
        let mut change_direction = start_direction;
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
        // Check consistency.
        assert!(3 <= cycle.len());
        assert_eq!(cycle.first(), cycle.last());
        let (from, to) = (cycle[0], cycle[cycle.len() - 2]);
        let first_node = cycle[1];
        let is_between_onedegree =
            self.one_degree_nodes.contains(&from) && self.one_degree_nodes.contains(&to);
        let is_between_node = from / 2 == to / 2;
        let starts_with_node = from / 2 == first_node / 2;
        let is_consistent = match is_between_onedegree || is_between_node || starts_with_node {
            true => start_direction == change_direction,
            false => start_direction != change_direction,
        };
        match is_consistent {
            true => score,
            false => score + LARGE_VALUE,
        }
    }
    fn update_copy_number_by_cycle_from(&mut self, cycle: &[usize], start_direction: bool) {
        let mut change_direction = start_direction;
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
}

#[cfg(test)]
pub mod tests {
    use super::*;
    fn mock_data_1() -> Graph {
        // let nodes = HashMap::new();
        let edges = vec![
            FatEdge::new(0, 1, 2, 1),
            FatEdge::new(1, 2, 1, 1),
            FatEdge::new(1, 4, 1, 1),
            FatEdge::new(2, 3, 1, 1),
            FatEdge::new(4, 5, 1, 1),
            FatEdge::new(3, 6, 1, 1),
            FatEdge::new(5, 6, 1, 1),
            FatEdge::new(6, 7, 2, 1),
            FatEdge::new(7, 8, 1, 1),
            FatEdge::new(7, 10, 1, 1),
            FatEdge::new(8, 9, 1, 1),
            FatEdge::new(10, 11, 1, 1),
            FatEdge::new(9, 12, 1, 1),
            FatEdge::new(12, 11, 1, 1),
            FatEdge::new(13, 12, 2, 1),
        ];
        let mut graph = (0..14).map(|_| Vec::with_capacity(2)).collect::<Vec<_>>();
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then_some(i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 13]);
        Graph {
            hap_coverage: 1f64,
            edges,
            graph,
            self_loops: vec![],
            one_degree_nodes,
        }
    }
    fn mock_data_2() -> Graph {
        let edges = vec![
            FatEdge::new(0, 1, 2, 1),
            FatEdge::new(1, 2, 2, 1),
            FatEdge::new(2, 3, 3, 1),
            FatEdge::new(3, 4, 1, 1),
            FatEdge::new(3, 6, 2, 1),
            FatEdge::new(4, 5, 1, 1),
            FatEdge::new(5, 6, 1, 1),
            FatEdge::new(6, 7, 3, 1),
            FatEdge::new(7, 8, 2, 1),
            FatEdge::new(7, 10, 1, 1),
            FatEdge::new(8, 9, 2, 1),
            FatEdge::new(9, 10, 2, 1),
            FatEdge::new(10, 11, 3, 1),
            FatEdge::new(11, 2, 1, 1),
            FatEdge::new(11, 12, 2, 1),
            FatEdge::new(12, 13, 2, 1),
        ];
        let mut graph = (0..14).map(|_| Vec::with_capacity(2)).collect::<Vec<_>>();
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then_some(i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 13]);
        Graph {
            hap_coverage: 1f64,
            edges,
            graph,
            self_loops: vec![],
            one_degree_nodes,
        }
    }
    fn mock_data_3() -> Graph {
        // let nodes = HashMap::new();
        let edges: Vec<_> = vec![
            FatEdge::new(0, 1, 20, 1),
            FatEdge::new(1, 2, 5, 1),
            FatEdge::new(2, 3, 5, 1),
            FatEdge::new(6, 3, 5, 1),
            FatEdge::new(1, 4, 5, 1),
            FatEdge::new(4, 5, 5, 1),
            FatEdge::new(5, 6, 5, 1),
            FatEdge::new(6, 7, 20, 1),
        ];
        let mut graph = (0..8).map(|_| Vec::with_capacity(2)).collect::<Vec<_>>();
        for edge in edges.iter() {
            graph[edge.from].push(LightEdge::new(edge.to));
            graph[edge.to].push(LightEdge::new(edge.from));
        }
        let one_degree_nodes: Vec<_> = graph
            .iter()
            .enumerate()
            .filter_map(|(i, eds)| (eds.len() == 1).then_some(i))
            .collect();
        assert_eq!(one_degree_nodes, vec![0, 7]);
        Graph {
            hap_coverage: 10f64,
            edges,
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
        eprintln!("{:?}", g.clone_copy_number());
        for edge in g.edges.iter() {
            eprintln!("{edge:?}");
        }
        for edge in g.edges.iter() {
            assert_eq!(edge.copy_number as usize, edge.target);
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
