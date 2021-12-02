use definitions::*;
use rand::{seq::SliceRandom, Rng, SeedableRng};
use rand_xoshiro::Xoroshiro128StarStar;
use std::collections::HashMap;
use std::io::BufReader;
const PROB_CORRECT: f64 = 0.9;
fn main() {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let unit_len = 2_000;
    let graph = Graph::new(&ds.encoded_reads, unit_len);
    let hap_cov = 20f64;
    let (node_cp, edge_cp) = graph.estimate_copy_number(hap_cov, 3903498039);
    println!("NODE\tunit\tcluster\tcopy.number");
    for ((unit, cluster), cp) in node_cp {
        println!("NODE\t{}\t{}\t{}", unit, cluster, cp);
    }
    println!("EDGE\tfrom.unit\tfrom.cluster\tto.unit\tto.cluster\tcopy.number");
    for (((fu, fc), (tu, tc)), cp) in edge_cp {
        println!("EDGE\t{}\t{}\t{}\t{}\t{}", fu, fc, tu, tc, cp);
    }
}

type Node = (u64, u64);
type Edge = (Node, Node);
#[derive(Debug, Clone)]
pub struct Graph {
    // Node->Node Index.
    node_to_idx: HashMap<Node, usize>,
    // Edge to Edge Index.
    edge_to_idx: HashMap<Edge, usize>,
    // weight of the nodes.
    nodes: Vec<f64>,
    // weight of the edges.
    edges: Vec<f64>,
    // edge index i -> (u,v). No direction needed.
    edge_terminals: Vec<(usize, usize)>,
    // node index i -> set of edge index j such that the j-th edge is connected to the "downstream" direction of nodes[i].
    downstream_edges: Vec<Vec<usize>>,
    // node index i -> set of edge index j such that the j-th edge is connected to the "upstream" direction of nodes[i].
    upstream_edges: Vec<Vec<usize>>,
}

impl Graph {
    fn serialize_node<T: std::borrow::Borrow<EncodedRead>>(reads: &[T]) -> HashMap<Node, usize> {
        let mut nodes: HashMap<Node, usize> = HashMap::new();
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            let len = nodes.len();
            nodes.entry((node.unit, node.cluster)).or_insert(len);
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
        let u = (u.unit, u.cluster);
        let v = (v.unit, v.cluster);
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
            let from = node_to_idx[&(w[0].unit, w[0].cluster)];
            let to = node_to_idx[&(w[1].unit, w[1].cluster)];
            let edge_idx = edge_to_idx[&Self::normalize(&w[0], &w[1])];
            match w[0].is_forward {
                true => downstream[from].push(edge_idx),
                false => upstream[from].push(edge_idx),
            }
            match w[1].is_forward {
                true => upstream[to].push(edge_idx),
                false => downstream[to].push(edge_idx),
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
    fn new<T: std::borrow::Borrow<EncodedRead>>(reads: &[T], unit_len: usize) -> Self {
        let node_to_idx = Self::serialize_node(reads);
        let edge_to_idx = Self::serialize_edge(reads);
        let (downstream_edges, upstream_edges) =
            Self::get_connections(&node_to_idx, &edge_to_idx, reads);
        let edge_terminals = Self::get_terminals(&node_to_idx, &edge_to_idx);
        let mut nodes: Vec<_> = vec![0f64; node_to_idx.len()];
        let mut edges_len: Vec<_> = vec![0; edge_to_idx.len()];
        let mut edges: Vec<_> = vec![0f64; edge_to_idx.len()];
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            nodes[node_to_idx[&(node.unit, node.cluster)]] += 1f64;
        }
        for read in reads.iter().map(|r| r.borrow()) {
            for (i, edge) in read.edges.iter().enumerate() {
                assert_eq!(edge.from, read.nodes[i].unit);
                assert_eq!(edge.to, read.nodes[i + 1].unit);
            }
        }
        for read in reads.iter().map(|r| r.borrow()) {
            for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
                let idx = edge_to_idx[&Self::normalize(&w[0], &w[1])];
                edges[idx] += 1f64;
                edges_len[idx] += edge.offset;
            }
        }
        // Tune the count of the edges.
        let lens: Vec<_> = reads.iter().map(|r| r.borrow().original_length).collect();
        let calibrator = CoverageCalibrator::new(&lens);
        nodes
            .iter_mut()
            .for_each(|x| *x = calibrator.calib(*x, unit_len));
        edges.iter_mut().zip(edges_len.iter()).for_each(|(x, len)| {
            let count = x.ceil() as i64;
            let gap_len = (*len / count + 2 * unit_len as i64).max(0) as usize;
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
    // TODO:Maybe configurabale?
    fn estimate_copy_number(
        &self,
        hap_cov: f64,
        seed: u64,
    ) -> (Vec<(Node, usize)>, Vec<(Edge, usize)>) {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(seed);
        let mut node_cp: Vec<_> = self
            .nodes
            .iter()
            .map(|x| Self::round(&mut rng, x / hap_cov))
            .collect();
        let mut edge_cp: Vec<_> = self
            .edges
            .iter()
            .map(|x| Self::round(&mut rng, x / hap_cov))
            .collect();
        let outer_times = 1000;
        for _ in 0..outer_times {
            self.update_nodes(&mut node_cp, &edge_cp, hap_cov, &mut rng);
            self.update_edges(&node_cp, &mut edge_cp, hap_cov, &mut rng);
        }
        let node_cp: Vec<_> = self
            .node_to_idx
            .iter()
            .map(|(&node, &i)| (node, node_cp[i]))
            .collect();
        let edge_cp: Vec<_> = self
            .edge_to_idx
            .iter()
            .map(|(&edge, &i)| (edge, edge_cp[i]))
            .collect();
        (node_cp, edge_cp)
    }
    fn update_nodes<R: Rng>(
        &self,
        node_cp: &mut [usize],
        edge_cp: &[usize],
        hap_cov: f64,
        rng: &mut R,
    ) {
        let edges = self.downstream_edges.iter().zip(self.upstream_edges.iter());
        for (((down_edge, up_edge), &w), cp) in edges.zip(self.nodes.iter()).zip(node_cp.iter_mut())
        {
            let down_cp: usize = down_edge.iter().map(|&j| edge_cp[j]).sum();
            let up_cp: usize = up_edge.iter().map(|&j| edge_cp[j]).sum();
            let mut choises = vec![];
            for cp in [down_cp, up_cp] {
                choises.push((cp, PROB_CORRECT * poisson(w, cp, hap_cov)));
                if down_cp == 0 {
                    choises.push((cp + 1, (1f64 - PROB_CORRECT) * poisson(w, cp + 1, hap_cov)));
                } else {
                    let minus_one = (1f64 - PROB_CORRECT) / 2f64 * poisson(w, cp - 1, hap_cov);
                    choises.push((cp - 1, minus_one));
                    let plus_one = (1f64 - PROB_CORRECT) / 2f64 * poisson(w, cp + 1, hap_cov);
                    choises.push((cp + 1, plus_one));
                }
            }
            let sum: f64 = choises.iter().map(|x| x.1).sum();
            let (new_cp, _) = choises.choose_weighted(rng, |&(_, w)| w / sum).unwrap();
            *cp = *new_cp;
        }
    }
    fn update_edges<R: Rng>(
        &self,
        node_cp: &[usize],
        edge_cp: &mut [usize],
        hap_cov: f64,
        rng: &mut R,
    ) {
        let mut ordering: Vec<_> = (0..edge_cp.len()).collect();
        ordering.shuffle(rng);
        for edge in ordering {
            let (from, to) = self.edge_terminals[edge];
            let from_cp = node_cp[from] - 
        }
    }
}

// Poisson(obs|copy_num*coverage)~Norm(obs|mean=var=copy_num*coverage)
// If copy number is zero, It is Norm(obs|mean=0,var=coverage/2), this is rough heuristics.
fn poisson(obs: f64, copy_num: usize, coverage: f64) -> f64 {
    let (mean, var) = match copy_num {
        0 => (0f64, coverage / 2f64),
        _ => (copy_num as f64 * coverage, copy_num as f64 * coverage),
    };
    -(obs - mean).powi(2) / 2f64 / var - (2f64 * std::f64::consts::PI * var).ln() / 2f64
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
    pub fn calib(&self, observed: f64, gap_len: usize) -> f64 {
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
