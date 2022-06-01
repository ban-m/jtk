// To me: do not forget update lightedge if you update fatedges
use rand::Rng;
use std::collections::{HashMap, HashSet};
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
    edges: Vec<FatEdge>,
    nodes: HashMap<Node, usize>,
    graph: Vec<Vec<LightEdge>>,
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
    penalty_diff: f64,
}

impl LightEdge {
    fn new(to: usize) -> Self {
        Self {
            to,
            is_in_mst: false,
            penalty_diff: 0f64,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FatEdge {
    from: usize,
    to: usize,
    target: f64,
    copy_number: usize,
    is_in_mst: bool,
    penalty_diff: f64,
}

impl FatEdge {
    pub fn new(from: usize, to: usize, target: f64) -> Self {
        let (from, to) = (from.min(to), from.max(to));
        Self {
            from,
            to,
            target,
            copy_number: 0,
            is_in_mst: false,
            penalty_diff: 0f64,
        }
    }
}

type Node = ((u64, u64), super::Position);
type Edge = (Node, Node);
impl Graph {
    pub fn new(nodes: HashMap<Node, usize>, edges: Vec<FatEdge>) -> Self {
        let mut graph = vec![Vec::with_capacity(2); nodes.len()];
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
        Self {
            nodes,
            edges,
            graph,
            one_degree_nodes,
        }
    }
    const LOOPTIMES: usize = 100;
    const LARGE_VALUE: f64 = 1000000f64;
    // TODO:Which is better, minimum penalty vs mcmc posterior prob?
    // I think this is not proper MCMC (as proposed distribution does not satisfy
    // detailed balanced condition), so I would like to retin the minimum penalty assignment.
    pub fn update_copy_numbers<R: Rng>(&mut self, rng: &mut R, config: &MSTConfig) {
        // for (i, edges) in self.graph.iter().enumerate() {
        //     let tos: Vec<_> = edges.iter().map(|e| e.to).collect();
        //     println!("{i}\t{tos:?}");
        // }
        let mut dfs_stack = vec![];
        let mut dfs_arrived_status = vec![Status::Undiscovered; self.graph.len()];
        let dfs_stack = &mut dfs_stack;
        let dfs_arrived_status = dfs_arrived_status.as_mut_slice();
        let (mut current_min, mut argmin) = (self.penalty(), self.edges.clone());
        for _ in 0..Self::LOOPTIMES {
            let to_increase = rng.gen_bool(0.5);
            self.update_mst(to_increase);
            let (optimal_cycle, penalty_diff) =
                self.find_optimal_cycle(dfs_stack, dfs_arrived_status);
            let prob = (-penalty_diff / config.temperature).exp().min(1f64);
            trace!("PROPOSED\t{prob:.4}\t{optimal_cycle:?}");
            if rng.gen_bool(prob) {
                self.update_copy_number_by_cycle(optimal_cycle, to_increase);
            }
            let penalty = self.penalty();
            if penalty < current_min {
                (current_min, argmin) = (penalty, self.edges.clone());
            }
        }
        self.edges = argmin;
    }
    fn penalty(&self) -> f64 {
        self.edges
            .iter()
            .map(|e| (e.target - e.copy_number as f64).powi(2))
            .sum()
    }
    // Search minimum-spanning-tree.
    fn update_mst(&mut self, to_increase: bool) {
        // Remove tempolary values
        self.edges.iter_mut().for_each(|edge| {
            edge.is_in_mst = false;
            let current_penalty = (edge.target - edge.copy_number as f64).powi(2);
            let proposed_penalty = match to_increase {
                true => (edge.target - edge.copy_number as f64 - 1f64).powi(2),
                false if edge.copy_number == 0 => Self::LARGE_VALUE,
                false => (edge.target + 1f64 - edge.copy_number as f64).powi(2),
            };
            edge.penalty_diff = proposed_penalty - current_penalty;
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
        let tree_edges = self.edges.iter().filter(|e| e.is_in_mst).count();
        assert_eq!(tree_edges + 1, self.graph.len());
        // Update lightedge
        for edge in self.edges.iter() {
            let ledge = self.graph[edge.from]
                .iter_mut()
                .find(|e| e.to == edge.to)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff = edge.penalty_diff;
            let ledge = self.graph[edge.to]
                .iter_mut()
                .find(|e| e.to == edge.from)
                .unwrap();
            ledge.is_in_mst = edge.is_in_mst;
            ledge.penalty_diff = edge.penalty_diff;
        }
    }
    fn find_optimal_cycle(
        &self,
        stack: &mut Vec<usize>,
        status: &mut [Status],
    ) -> (Vec<usize>, f64) {
        let (mut current_min, mut argmin) = (Self::LARGE_VALUE, vec![]);
        for edge in self.edges.iter().filter(|e| !e.is_in_mst) {
            let cycle = self
                .find_cycle_between(edge.from, edge.to, stack, status)
                .unwrap();
            let penalty = self.penalty_of_cycle(&cycle);
            if penalty < current_min {
                current_min = penalty;
                argmin = cycle;
            }
        }
        for (i, &from) in self.one_degree_nodes.iter().enumerate() {
            for &to in self.one_degree_nodes.iter().skip(i + 1) {
                let cycle = self.find_cycle_between(from, to, stack, status).unwrap();
                let penalty = self.penalty_of_cycle(&cycle);
                if penalty < current_min {
                    current_min = penalty;
                    argmin = cycle;
                }
            }
        }
        (argmin, current_min)
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
            status[last] = Status::Processing;
            let edge_in_tree = self.graph[last].iter().filter(|e| e.is_in_mst);
            let edge_undiscovered = edge_in_tree.filter(|e| status[e.to] == Status::Undiscovered);
            for take_edge in edge_undiscovered {
                if take_edge.to == to {
                    // Find target node.
                    let mut cycle = Vec::with_capacity(stack.len() + 2);
                    cycle.extend(stack.iter());
                    cycle.push(to);
                    cycle.push(from);
                    return Some(cycle);
                } else {
                    stack.push(take_edge.to);
                    continue 'dfs;
                }
            }
            let fin_node = stack.pop().unwrap();
            status[fin_node] = Status::Processed;
        }
        None
    }
    fn penalty_of_cycle(&self, cycle: &[usize]) -> f64 {
        cycle
            .windows(2)
            .map(|w| {
                let (from, to) = (w[0], w[1]);
                // If there's no (from,to) edge, it SHOULD be one-degree bridge.
                self.graph[from]
                    .iter()
                    .find(|e| e.to == to)
                    .map(|e| e.penalty_diff)
                    .unwrap_or(0f64)
            })
            .sum()
    }
    fn update_copy_number_by_cycle(&mut self, cycle: Vec<usize>, to_increase: bool) {
        let cycle: HashSet<_> = cycle
            .windows(2)
            .map(|w| (w[0].min(w[1]), w[0].max(w[1])))
            .collect();
        self.edges
            .iter_mut()
            .filter(|e| cycle.contains(&(e.from, e.to)))
            .for_each(|e| match to_increase {
                true => e.copy_number += 1,
                false => e.copy_number -= 1,
            });
    }
    pub fn edge_copy_numbers(&self) -> Vec<(Edge, usize)> {
        // We need reverse index ...
        let mut reverse_index = vec![((0, 0), super::Position::Head); self.nodes.len()];
        for (&node, &idx) in self.nodes.iter() {
            reverse_index[idx] = node;
        }
        self.edges
            .iter()
            .map(|edge| {
                let key = (reverse_index[edge.from], reverse_index[edge.to]);
                (key, edge.copy_number)
            })
            .collect()
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    fn mock_data_1() -> Graph {
        let nodes = HashMap::new();
        let mut edges: Vec<_> = vec![];
        edges.push(FatEdge::new(0, 1, 2f64));
        edges.push(FatEdge::new(1, 2, 1f64));
        edges.push(FatEdge::new(1, 4, 1f64));
        edges.push(FatEdge::new(2, 3, 1f64));
        edges.push(FatEdge::new(4, 5, 1f64));
        edges.push(FatEdge::new(3, 6, 1f64));
        edges.push(FatEdge::new(5, 6, 1f64));
        edges.push(FatEdge::new(6, 7, 2f64));
        edges.push(FatEdge::new(7, 8, 1f64));
        edges.push(FatEdge::new(7, 10, 1f64));
        edges.push(FatEdge::new(8, 9, 1f64));
        edges.push(FatEdge::new(10, 11, 1f64));
        edges.push(FatEdge::new(9, 12, 1f64));
        edges.push(FatEdge::new(12, 11, 1f64));
        edges.push(FatEdge::new(13, 12, 2f64));
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
            edges,
            nodes,
            graph,
            one_degree_nodes,
        }
    }
    fn mock_data_2() -> Graph {
        let nodes = HashMap::new();
        let mut edges: Vec<_> = vec![];
        edges.push(FatEdge::new(0, 1, 2f64));
        edges.push(FatEdge::new(1, 2, 2f64));
        edges.push(FatEdge::new(2, 3, 3f64));
        edges.push(FatEdge::new(3, 4, 1f64));
        edges.push(FatEdge::new(3, 6, 2f64));
        edges.push(FatEdge::new(4, 5, 1f64));
        edges.push(FatEdge::new(5, 6, 1f64));
        edges.push(FatEdge::new(6, 7, 3f64));
        edges.push(FatEdge::new(7, 8, 2f64));
        edges.push(FatEdge::new(7, 10, 1f64));
        edges.push(FatEdge::new(8, 9, 2f64));
        edges.push(FatEdge::new(9, 10, 2f64));
        edges.push(FatEdge::new(10, 11, 3f64));
        edges.push(FatEdge::new(11, 2, 1f64));
        edges.push(FatEdge::new(11, 12, 2f64));
        edges.push(FatEdge::new(12, 13, 2f64));
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
            edges,
            nodes,
            graph,
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
            assert!((edge.copy_number as f64 - edge.target).abs() < 0.01);
        }
    }
}
