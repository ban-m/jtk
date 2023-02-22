pub mod copy_number_barrier;
pub mod copy_number_flow;
pub mod copy_number_gibbs;
pub mod copy_number_mrf;
pub mod copy_number_mst;

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
    pub fn calib(&self, observed: u64, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        let factor = self.cum_sum[idx] - gap_len * (self.lengths.len() - idx);
        observed as f64 * self.mean / (factor as f64 / self.lengths.len() as f64)
    }
    /// Calibrate the observed coverage in the gap-len bp region
    /// into the "actual" coverage, regarding that region as a single-bp region.
    /// Note that if the gap_len is longer than the longest reads in the dataset,
    /// it returns zero.
    pub fn calib_f64(&self, observed: f64, gap_len: usize) -> f64 {
        let idx = match self.lengths.binary_search(&gap_len) {
            Ok(idx) => idx,
            Err(idx) if idx == self.lengths.len() => return 0f64,
            Err(idx) => idx,
        };
        let factor = self.cum_sum[idx] - gap_len * (self.lengths.len() - idx);
        observed * self.mean / (factor as f64 / self.lengths.len() as f64)
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

// #[derive(Debug, Clone)]
// pub struct GraphWithCoverage {
//     nodes: Vec<NodeWithCoverage>,
//     edges: Vec<Vec<EdgeWithCoverage>>,
//     haploid_coverage: f64,
// }

// #[derive(Debug, Clone, Copy)]
// pub struct NodeWithCoverage {
//     length: usize,
//     coverage: f64,
// }

// #[derive(Debug, Clone, Copy)]
// pub struct EdgeWithCoverage {
//     from: usize,
//     from_position: bool,
//     to: usize,
//     to_position: bool,
// }

// #[derive(Debug, Clone)]
// pub struct GibbsSampler {
//     nodes: Vec<u64>,
//     edges: Vec<GbsEdge>,
//     // i->(indices of edges with (i,true), indices of edges with (i,false));
//     node_terminals: Vec<(Vec<usize>, Vec<usize>)>,
//     haploid_coverage: f64,
// }

// #[derive(Debug, Clone)]
// pub struct Graph {
//     // Node->Node Index.
//     node_to_idx: HashMap<Node, usize>,
//     // Edge to Edge Index.
//     edge_to_idx: HashMap<Edge, usize>,
//     // weight of the nodes.
//     nodes: Vec<u64>,
//     // weight of the edges.
//     edges: Vec<u64>,
//     // edge index i -> (u,v). No direction needed.
//     edge_terminals: Vec<(usize, usize)>,
//     // node index i -> set of edge index j such that the j-th edge is connected to the "downstream" direction of nodes[i].
//     downstream_edges: Vec<Vec<usize>>,
//     // node index i -> set of edge index j such that the j-th edge is connected to the "upstream" direction of nodes[i].
//     upstream_edges: Vec<Vec<usize>>,
// }

// #[derive(Debug, Clone)]
// pub struct Graph {
//     nodes: Vec<RawNode>,
//     node_copy_numbers: Vec<usize>,
//     edges: Vec<RawEdge>,
//     edge_copy_numbers: Vec<usize>,
//     residual_graph: ResGraph,
//     hap_cov: f64,
// }

// #[derive(Debug, Clone)]
// pub struct Graph {
//     // (u,u_pos, v, v_pos).
//     edges: Vec<(usize, bool, usize, bool)>,
//     // coverage of nodes, length(number of chunks) in nodes.
//     coverages: Vec<(u64, usize)>,
//     // for each node, for each position, return the set of edge indices.
//     edge_lists: Vec<[Vec<usize>; 2]>,
// }

// #[derive(Debug, Clone)]
// pub struct Graph {
//     hap_coverage: f64,
//     edges: Vec<FatEdge>,
//     graph: Vec<Vec<LightEdge>>,
//     self_loops: Vec<FatEdge>,
//     one_degree_nodes: Vec<usize>,
// }
