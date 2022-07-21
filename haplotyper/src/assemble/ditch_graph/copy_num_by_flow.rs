use rand::Rng;

use crate::find_union::FindUnion;

// Four bit packing
// 4i, 4i+1, 4i+2, 4i+3 as follows
// Head ----------> Tail
// into
//     4i   ----->  4i+1
//     4i+2 <-----  4i+3
// In general,
// 4i, 4i+3: Nodes connected to nodes of other region.
// 4i+1, 4i+2: Nodes connected internally.
// Note that the reverse/forward edges can be decided by only see the indices.

const LARGE_VALUE: f64 = 10_000_000_000f64;
type RawNode = (f64, usize);
type RawEdge = (usize, bool, usize, bool, f64);

#[derive(Debug, Clone)]
pub struct Graph {
    nodes: Vec<RawNode>,
    node_copy_numbers: Vec<usize>,
    edges: Vec<RawEdge>,
    edge_copy_numbers: Vec<usize>,
    residual_graph: ResGraph,
    hap_cov: f64,
}

#[derive(Debug, Clone)]
struct ResGraph {
    nodes: Vec<ResEdges>,
}

impl ResGraph {
    fn new(dim: usize) -> Self {
        Self {
            nodes: vec![ResEdges::default(); dim],
        }
    }
    fn len(&self) -> usize {
        self.nodes.len()
    }
}

impl std::ops::Index<ResIndex> for ResGraph {
    type Output = ResEdges;
    fn index(&self, ResIndex(idx): ResIndex) -> &Self::Output {
        self.nodes.get(idx).unwrap()
    }
}

impl std::ops::IndexMut<ResIndex> for ResGraph {
    fn index_mut(&mut self, ResIndex(idx): ResIndex) -> &mut Self::Output {
        self.nodes.get_mut(idx).unwrap()
    }
}

#[derive(Debug, Clone, Default)]
struct ResEdges(Vec<ResEdge>);

impl std::ops::Index<EdgeIndex> for ResEdges {
    type Output = ResEdge;
    fn index(&self, EdgeIndex(idx): EdgeIndex) -> &Self::Output {
        self.0.get(idx).unwrap()
    }
}

impl std::ops::IndexMut<EdgeIndex> for ResEdges {
    fn index_mut(&mut self, EdgeIndex(idx): EdgeIndex) -> &mut Self::Output {
        self.0.get_mut(idx).unwrap()
    }
}

impl ResEdges {
    fn push(&mut self, edge: ResEdge) {
        self.0.push(edge);
    }
    fn len(&self) -> usize {
        self.0.len()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ResIndex(usize);

#[derive(Debug, Clone, Copy)]
struct EdgeIndex(usize);

#[derive(Debug, Clone, Copy)]
struct ResEdge {
    from: ResIndex,
    to: ResIndex,
    // Currently this is not needed, because the capacity/penalty is
    // stored in the `target`
    #[allow(dead_code)]
    rev_idx: EdgeIndex,
    target: RawPointer,
}

impl ResEdge {
    fn new(from: ResIndex, to: ResIndex, rev_idx: EdgeIndex, target: RawPointer) -> Self {
        Self {
            from,
            to,
            rev_idx,
            target,
        }
    }
    // Return true if this is a back edge.
    fn is_back(&self) -> bool {
        let from = matches!(self.from.0 % 4, 0 | 3);
        let to = matches!(self.to.0 % 4, 0 | 3);
        assert_ne!(from, to);
        match (self.target, from, to) {
            (RawPointer::Node(_), true, false) => false,
            (RawPointer::Node(_), false, true) => true,
            (RawPointer::Edge(_), false, true) => false,
            (RawPointer::Edge(_), true, false) => true,
            (_, f, t) => panic!("{:?}\t{}\t{}", self, f, t),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum RawPointer {
    Node(usize),
    Edge(usize),
}

impl std::fmt::Display for RawPointer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RawPointer::Node(n) => write!(f, "N({n})"),
            RawPointer::Edge(n) => write!(f, "E({n})"),
        }
    }
}

impl Graph {
    pub fn new(nodes: &[(f64, usize)], edges: &[RawEdge], hap_cov: f64) -> Self {
        let mut graph = Self {
            nodes: nodes.to_vec(),
            edges: edges.to_vec(),
            edge_copy_numbers: vec![0; edges.len()],
            node_copy_numbers: vec![0; nodes.len()],
            residual_graph: ResGraph::new(nodes.len() * 4),
            hap_cov,
        };
        for (i, _) in nodes.iter().enumerate() {
            graph.add_inner_edge(4 * i, 4 * i + 1, i);
            graph.add_inner_edge(4 * i + 3, 4 * i + 2, i);
        }
        for (i, (f, fdir, t, tdir, _)) in edges.iter().enumerate() {
            let (fsource, fdest) = match fdir {
                true => (4 * f + 2, 4 * f),
                false => (4 * f + 1, 4 * f + 3),
            };
            let (tsource, tdest) = match tdir {
                true => (4 * t + 2, 4 * t),
                false => (4 * t + 1, 4 * t + 3),
            };
            graph.add_intra_edge(fsource, tdest, i);
            if fsource != tsource || fdest != tdest {
                assert!(fsource != tsource && fdest != tdest);
                graph.add_intra_edge(tsource, fdest, i);
            }
        }
        graph
    }
    fn add_intra_edge(&mut self, from: usize, to: usize, edge_idx: usize) {
        let edge_idx = RawPointer::Edge(edge_idx);
        let from = ResIndex(from);
        let to = ResIndex(to);
        let to_idx = EdgeIndex(self.residual_graph[to].len());
        let from_idx = EdgeIndex(self.residual_graph[from].len());
        self.residual_graph[from].push(ResEdge::new(from, to, to_idx, edge_idx));
        self.residual_graph[to].push(ResEdge::new(to, from, from_idx, edge_idx));
    }
    fn add_inner_edge(&mut self, from: usize, to: usize, node_idx: usize) {
        let node_idx = RawPointer::Node(node_idx);
        let from = ResIndex(from);
        let to = ResIndex(to);
        let to_idx = EdgeIndex(self.residual_graph[to].len());
        let from_idx = EdgeIndex(self.residual_graph[from].len());
        self.residual_graph[from].push(ResEdge::new(from, to, to_idx, node_idx));
        self.residual_graph[to].push(ResEdge::new(to, from, from_idx, node_idx));
    }
    pub fn optimize<R: Rng>(&mut self, _: &mut R) {
        // Clear copy numbers.
        trace!("FLOW\tSOURCES\t{}", self.source_sinks().len());
        self.node_copy_numbers.iter_mut().for_each(|n| *n = 0);
        self.edge_copy_numbers.iter_mut().for_each(|n| *n = 0);
        trace!("PENALTY\t{:.1}", self.penalty());
        while self.update() {
            if log_enabled!(log::Level::Trace) {
                trace!("PENALTY\t{:.1}", self.penalty());
            }
        }
        trace!("PENALTY\t{:.1}", self.penalty());
    }
    fn penalty(&self) -> f64 {
        let nodes: f64 = self
            .nodes
            .iter()
            .zip(self.node_copy_numbers.iter())
            .map(|(&(target, len), &cp)| len as f64 * (target - cp as f64 * self.hap_cov).powi(2))
            .sum();
        let edges: f64 = self
            .edges
            .iter()
            .zip(self.edge_copy_numbers.iter())
            .map(|(edge, &cp)| (edge.4 - cp as f64 * self.hap_cov).powi(2))
            .sum();
        nodes + edges
    }
    fn update(&mut self) -> bool {
        // Here we have two way:
        // 1. Do Warshal Floyd to get the minimum penalty of each pair of nodes,
        // selecting the minimum pair of the nodes, and traceback.
        // It requires O(V^3) exactly.
        // 2. Do Bellman-Ford on each pair of the sources and sinks.
        // Usually, # of sources and sinks are small. So it is roughly O(VE). I choose this.
        let mut fu = self.connected_components();
        let source_sinks = self.source_sinks();
        // Maybe we need to pick the minimum one...
        for &source in source_sinks.iter() {
            for &sink in source_sinks.iter() {
                if !fu.same(source.0, sink.0).unwrap() {
                    continue;
                }
                // Bellman-Ford.
                if self.update_between(source, sink) {
                    return true;
                }
            }
        }
        if source_sinks.is_empty() {
            use std::collections::BTreeMap;
            let mut clusters: BTreeMap<_, Vec<_>> = BTreeMap::new();
            for i in 0..self.residual_graph.len() {
                clusters
                    .entry(fu.find(i).unwrap())
                    .or_default()
                    .push(ResIndex(i));
            }
            for nodes in clusters.values().filter(|ns| ns.len() > 2) {
                if self.update_between(nodes[0], nodes[1]) {
                    return true;
                }
            }
        }
        false
    }
    fn connected_components(&self) -> crate::find_union::FindUnion {
        let mut fu = FindUnion::new(self.residual_graph.len());
        for edges in self.residual_graph.nodes.iter() {
            for edge in edges.0.iter() {
                fu.unite(edge.from.0, edge.to.0);
            }
        }
        fu
    }
    fn update_between(&mut self, source: ResIndex, sink: ResIndex) -> bool {
        // TODO: First mark the nodes reachable from source and to sink, then
        // update only on these bi-directionally reachable nodes.
        // It would speed up some...?
        // TODO: Memoize the `self.score(edge)` function. It would make this function 2x faster?
        // Bellman Ford.
        let len = self.residual_graph.len();
        // Init.
        let mut dists = vec![LARGE_VALUE; len];
        let mut predv = vec![ResIndex(len + 1); len];
        let mut prede = vec![None; len];
        dists[source.0] = 0f64;
        let edge_scores: Vec<Vec<_>> = self
            .residual_graph
            .nodes
            .iter()
            .map(|edges| edges.0.iter().map(|e| self.score(e)).collect())
            .collect();
        // Loop.
        for _ in 0..len - 1 {
            for (edges, scores) in self.residual_graph.nodes.iter().zip(edge_scores.iter()) {
                assert_eq!(scores.len(), edges.len());
                for (edge, score) in std::iter::zip(edges.0.iter(), scores.iter()) {
                    let (from, to) = (edge.from.0, edge.to.0);
                    if dists[from] < LARGE_VALUE && dists[from] + score < dists[to] {
                        dists[to] = dists[from] + score;
                        predv[to] = edge.from;
                        prede[to] = Some(*edge);
                    }
                }
            }
        }
        // Check Negative cycle.
        for (edges, scores) in std::iter::zip(&self.residual_graph.nodes, &edge_scores) {
            assert_eq!(edges.len(), scores.len());
            for (edge, score) in std::iter::zip(&edges.0, scores) {
                if dists[edge.from.0] < LARGE_VALUE && dists[edge.from.0] + score < dists[edge.to.0]
                {
                    // negative loop. Return if it indeed decrease the penalty.
                    let cycle = Self::traverse_cycle(&predv, &prede, edge.from.0);
                    let score = self.eval(&cycle);
                    if score < 0f64 {
                        trace!("UPDATE\tCYCLE\t{score:.2}\t{}", cycle.len());
                        self.update_by(&cycle);
                        return true;
                    }
                }
            }
        }
        if dists[sink.0] < 0f64 {
            // Recover shortest path from the source to the sink.
            let mut current = sink;
            let mut path = vec![prede[current.0].unwrap()];
            current = predv[current.0];
            while current != source {
                path.push(prede[current.0].unwrap());
                current = predv[current.0];
            }
            path.reverse();
            let score = self.eval(&path);
            if score < 0f64 {
                trace!("UPDATE\tPATH\t{score:.2}\t{}", path.len());
                self.update_by(&path);
                return true;
            }
        }
        false
    }
    fn eval(&self, path: &[ResEdge]) -> f64 {
        let mut node_diff = vec![0; self.nodes.len()];
        let mut edge_diff = vec![0; self.edges.len()];
        for edge in path.iter() {
            let cp_diff = match edge.is_back() {
                true => -1,
                false => 1,
            };
            match edge.target {
                RawPointer::Node(idx) => node_diff[idx] += cp_diff,
                RawPointer::Edge(idx) => edge_diff[idx] += cp_diff,
            }
        }
        let node_score: f64 = node_diff
            .iter()
            .zip(self.node_copy_numbers.iter())
            .zip(self.nodes.iter())
            .filter(|&((&diff, _), _)| diff != 0)
            .map(|((diff, &current), &(target, len))| {
                let old_penalty = (target - current as f64 * self.hap_cov).powi(2);
                let cp = current as i64 + diff;
                let new_pen = match 0 <= cp {
                    true => (target - cp as f64 * self.hap_cov).powi(2),
                    false => LARGE_VALUE,
                };
                (new_pen - old_penalty) * len as f64
            })
            .sum();
        let edge_score: f64 = edge_diff
            .iter()
            .zip(self.edge_copy_numbers.iter())
            .zip(self.edges.iter())
            .filter(|&((&diff, _), _)| diff != 0)
            .map(|((diff, &current), target)| {
                let target = target.4;
                let old_pen = (target - current as f64 * self.hap_cov).powi(2);
                let cp = current as i64 + diff;
                let new_pen = match 0 <= cp {
                    true => (target - cp as f64 * self.hap_cov).powi(2),
                    false => LARGE_VALUE,
                };
                new_pen - old_pen
            })
            .sum();
        node_score + edge_score
    }
    fn update_by(&mut self, path: &[ResEdge]) {
        for edge in path {
            match (edge.target, edge.is_back()) {
                (RawPointer::Node(idx), true) => self.node_copy_numbers[idx] -= 1,
                (RawPointer::Node(idx), false) => self.node_copy_numbers[idx] += 1,
                (RawPointer::Edge(idx), true) => self.edge_copy_numbers[idx] -= 1,
                (RawPointer::Edge(idx), false) => self.edge_copy_numbers[idx] += 1,
            }
        }
        // If copy number becomes negative, underflow occurs.
        assert!(self.node_copy_numbers.iter().all(|&x| x < 100000));
        assert!(self.edge_copy_numbers.iter().all(|&x| x < 100000));
    }
    // Extract cycle from `pred`. Assume there is at least one.
    // The returned cycle does not always start from `start` node.
    fn traverse_cycle(predv: &[ResIndex], prede: &[Option<ResEdge>], start: usize) -> Vec<ResEdge> {
        // Trace back `V` times to get into the cycle.
        let mut current = ResIndex(start);
        for _ in 0..predv.len() + 3 {
            current = predv[current.0];
        }
        // Trace back until re-entry `current.`
        let root = current;
        let mut edges = vec![prede[current.0].unwrap()];
        current = predv[current.0];
        while current != root {
            edges.push(prede[current.0].unwrap());
            current = predv[current.0];
        }
        edges.reverse();
        edges
    }
    fn source_sinks(&self) -> Vec<ResIndex> {
        // If a node is connected only inner edge, it is a source/sink node.
        self.residual_graph
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, edges)| {
                assert!(!edges.0.is_empty());
                (edges.0.len() == 1).then(|| ResIndex(i))
            })
            .collect()
    }
    // Score of the edge.
    fn score(&self, edge: &ResEdge) -> f64 {
        let to_decrease = edge.is_back();
        match edge.target {
            RawPointer::Node(idx) => self.score_node(idx, to_decrease),
            RawPointer::Edge(idx) => self.score_edge(idx, to_decrease),
        }
    }
    fn score_node(&self, idx: usize, to_decrease: bool) -> f64 {
        let (target, len) = self.nodes[idx];
        let cp = self.node_copy_numbers[idx];
        let old_error = (target - cp as f64 * self.hap_cov).powi(2);
        let new_error = match to_decrease {
            false => (target - (cp + 1) as f64 * self.hap_cov).powi(2),
            true if cp == 0 => LARGE_VALUE,
            true => (target - (cp - 1) as f64 * self.hap_cov).powi(2),
        };
        (new_error - old_error) * len as f64
    }
    fn score_edge(&self, idx: usize, to_decrease: bool) -> f64 {
        let target = self.edges[idx].4;
        let cp = self.edge_copy_numbers[idx];
        let old_error = (target - cp as f64 * self.hap_cov).powi(2);
        let new_error = match to_decrease {
            false => (target - (cp + 1) as f64 * self.hap_cov).powi(2),
            true if cp == 0 => LARGE_VALUE,
            true => (target - (cp - 1) as f64 * self.hap_cov).powi(2),
        };
        new_error - old_error
    }
    pub fn copy_numbers(&self) -> (Vec<usize>, Vec<usize>) {
        (
            self.node_copy_numbers.clone(),
            self.edge_copy_numbers.clone(),
        )
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};
    use rand_xoshiro::Xoroshiro128PlusPlus;
    #[test]
    fn is_back_test() {
        let edge = ResEdge::new(ResIndex(0), ResIndex(1), EdgeIndex(0), RawPointer::Node(0));
        assert!(!edge.is_back());
        let edge = ResEdge::new(ResIndex(1), ResIndex(0), EdgeIndex(0), RawPointer::Node(0));
        assert!(edge.is_back());
        let edge = ResEdge::new(ResIndex(2), ResIndex(3), EdgeIndex(0), RawPointer::Node(0));
        assert!(edge.is_back());
        let edge = ResEdge::new(ResIndex(3), ResIndex(2), EdgeIndex(0), RawPointer::Node(0));
        assert!(!edge.is_back());
        let edge = ResEdge::new(ResIndex(1), ResIndex(4), EdgeIndex(0), RawPointer::Edge(0));
        assert!(!edge.is_back());
        let edge = ResEdge::new(ResIndex(4), ResIndex(2), EdgeIndex(0), RawPointer::Edge(0));
        assert!(edge.is_back());
        let edge = ResEdge::new(ResIndex(2), ResIndex(7), EdgeIndex(0), RawPointer::Edge(0));
        assert!(!edge.is_back());
        let edge = ResEdge::new(ResIndex(5), ResIndex(7), EdgeIndex(0), RawPointer::Edge(0));
        assert!(!edge.is_back());
    }
    #[test]
    fn mock_data_1() {
        let nodes_cp = vec![2, 1, 1, 2, 1, 1, 2];
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(349823094);
        let cov = 10f64;
        let nodes: Vec<_> = nodes_cp
            .iter()
            .map(|&cp| {
                let cov = Normal::new(cov * cp as f64, 1f64).unwrap().sample(&mut rng);
                (cov, 2)
            })
            .collect();
        let edges: Vec<_> = vec![
            (0, false, 1, true, cov),
            (0, false, 2, true, cov),
            (1, false, 3, true, cov),
            (2, false, 3, true, cov),
            (3, false, 4, true, cov),
            (3, false, 5, true, cov),
            (4, false, 6, true, cov),
            (5, false, 6, true, cov),
        ];
        let edge_cp = vec![1; edges.len()];
        let mut graph = Graph::new(&nodes, &edges, cov);
        graph.optimize(&mut rng);
        let (pred, edge_pred) = graph.copy_numbers();
        assert_eq!(pred, nodes_cp);
        assert_eq!(edge_pred, edge_cp);
    }
    #[test]
    fn mock_data_2() {
        let nodes_cp = vec![3, 1, 3, 2, 1, 3];
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(349823094);
        let cov = 10f64;
        let nodes: Vec<_> = nodes_cp
            .iter()
            .map(|&cp| {
                let cov = Normal::new(cov * cp as f64, 1f64).unwrap().sample(&mut rng);
                (cov, 10)
            })
            .collect();
        let edges: Vec<_> = vec![
            (0, false, 1, true, cov),
            (1, false, 2, true, cov),
            (0, false, 2, true, 2f64 * cov),
            (2, false, 3, true, 2f64 * cov),
            (3, false, 5, true, 2f64 * cov),
            (2, false, 4, true, cov),
            (4, false, 5, true, cov),
            (5, false, 0, true, 3f64 * cov),
        ];
        let edge_cp = vec![1, 1, 2, 2, 2, 1, 1, 3];
        let mut graph = Graph::new(&nodes, &edges, cov);
        graph.optimize(&mut rng);
        let (pred, edge_pred) = graph.copy_numbers();
        assert_eq!(pred, nodes_cp);
        assert_eq!(edge_pred, edge_cp);
    }
    #[test]
    fn mock_data_3() {
        let nodes_cp = vec![2, 4, 2, 1, 1, 2, 2];
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(349823094);
        let cov = 10f64;
        let nodes: Vec<_> = nodes_cp
            .iter()
            .map(|&cp| {
                let cov = Normal::new(cov * cp as f64, 1f64).unwrap().sample(&mut rng);
                (cov, 10)
            })
            .collect();
        let edges: Vec<_> = vec![
            (0, false, 1, true, 2f64 * cov),
            (1, false, 2, true, 2f64 * cov),
            (2, false, 3, true, cov),
            (3, false, 5, true, cov),
            (5, false, 1, false, 2f64 * cov),
            (2, false, 4, true, cov),
            (4, false, 5, true, cov),
            (6, false, 1, true, 2f64 * cov),
        ];
        let edge_cp = vec![2, 2, 1, 1, 2, 1, 1, 2];
        let mut graph = Graph::new(&nodes, &edges, cov);
        graph.optimize(&mut rng);
        let (pred, edge_pred) = graph.copy_numbers();
        assert_eq!(pred, nodes_cp);
        assert_eq!(edge_pred, edge_cp);
    }
    #[test]
    fn mock_data_4() {
        let nodes_cp = vec![2, 1, 1, 2];
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(349823094);
        let cov = 30f64;
        let nodes = vec![(26.0, 1), (30.0, 1000), (30.0, 1000), (23.0, 1)];
        let edges: Vec<_> = vec![
            (0, false, 1, true, 26.0),
            (0, false, 2, true, 20.0),
            (1, false, 3, true, 25.0),
            (2, false, 3, true, 10.0),
        ];
        let edge_cp = vec![1, 1, 1, 1];
        let mut graph = Graph::new(&nodes, &edges, cov);
        graph.optimize(&mut rng);
        let (pred, edge_pred) = graph.copy_numbers();
        assert_eq!(pred, nodes_cp);
        assert_eq!(edge_pred, edge_cp);
    }
}
