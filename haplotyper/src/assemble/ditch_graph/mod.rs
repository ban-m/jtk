use super::AssembleConfig;
use definitions::{EncodedRead, Unit};
use rayon::prelude::*;
mod sequence_generation;
use std::collections::HashMap;
use std::collections::HashSet;
pub mod dg_test;
// const DRAFT_COVERAGE: usize = 20;
// const DRAFT_REP_NUM: usize = 7;

/// Position::Head is the up-stream position of a unit, and Tail is downstream by default.
#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash, PartialOrd, Ord)]
pub enum Position {
    Head,
    Tail,
}

impl std::fmt::Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let x = match self {
            Position::Head => 'H',
            Position::Tail => 'T',
        };
        write!(f, "{}", x)
    }
}

impl std::ops::Not for Position {
    type Output = Self;
    fn not(self) -> Self {
        match self {
            Position::Head => Position::Tail,
            Position::Tail => Position::Head,
        }
    }
}

type Node = (u64, u64);
type DitEdge = ((Node, Position), (Node, Position));
type EdgeConsensus = HashMap<DitEdge, EdgeLabel>;
type NodeConsensus = HashMap<Node, Vec<u8>>;

/// A focus from a node of a ditch graph.
/// Here, we tag a node with its position, i.e., Position.
#[derive(Debug, Clone)]
pub struct Focus {
    pub from: Node,
    pub from_position: Position,
    pub to: Node,
    pub to_position: Position,
    pub dist: usize,
    /// Log Likelihood ratio between the alt hypothesis/null hypothesis.
    pub log_likelihood_ratio: f64,
}

impl Focus {
    /// LogLilihood ratio
    pub fn llr(&self) -> f64 {
        self.log_likelihood_ratio
    }
    pub fn new(
        from: Node,
        from_position: Position,
        to: Node,
        to_position: Position,
        dist: usize,
        lk: f64,
    ) -> Self {
        Self {
            from,
            from_position,
            to,
            to_position,
            dist,
            log_likelihood_ratio: lk,
        }
    }
}

/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
/// It is bi-directed graph, in other words, if you see (from,from_pos, to, to_pos),
/// then, there is a edge (to, to_pos, from, from_pos).
#[derive(Clone)]
pub struct DitchGraph<'a> {
    nodes: HashMap<Node, DitchNode<'a>>,
}

impl<'a> std::fmt::Display for DitchGraph<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        use histgram_viz::Histgram;
        let nodes = self.nodes.len();
        let edges = self.nodes.values().map(|n| n.edges.len()).sum::<usize>();
        writeln!(f, "Node:{}, Edges:{}", nodes, edges)?;
        let occs: Vec<_> = self.nodes.values().map(|n| n.occ).collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Node Occs:{}", hist.format(20, 20))?;
        let occs: Vec<_> = self
            .nodes
            .values()
            .flat_map(|n| n.edges.iter().map(|e| e.occ))
            .collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Edge Occs:{}", hist.format(20, 20))?;
        let degrees = {
            let mut degs: HashMap<usize, usize> = HashMap::new();
            for node in self.nodes.values() {
                let degree = node.edges.len();
                *degs.entry(degree).or_default() += 1;
            }
            let mut degs = degs.into_iter().collect::<Vec<_>>();
            degs.sort_by_key(|x| x.0);
            degs.into_iter()
                .map(|(deg, count)| format!("{}:{}", deg, count))
                .collect::<Vec<_>>()
        };
        write!(f, "[{}]", degrees.join(","))
    }
}

impl<'a> std::fmt::Debug for DitchGraph<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edge = self.nodes.values().map(|e| e.edges.len()).sum::<usize>();
        writeln!(f, "Nodes:{}\tEdges:{}\n", self.nodes.len(), edge)?;
        let lines: Vec<_> = self
            .nodes
            .values()
            .map(|node| format!("{}", node))
            .collect();
        write!(f, "{}", lines.join("\n"))
    }
}

/// Ditch node.
/// It never allocates, or explicitly copy the nodes inside it,
/// rather, it borrows the reference to the `definitions::Node`s.
/// Note that even if a read contains (x,y,h) -> (z,w,t) connection,
/// the node of (z, w, t) contains an edge toward (x,y,h)-node.
/// Here, x and z are unit, y and w are cluster, and h and t are
/// head or tail.
#[derive(Debug, Clone)]
pub struct DitchNode<'a> {
    node: Node,
    occ: usize,
    seq: Vec<u8>,
    edges: Vec<DitchEdge>,
    // "Tip" of reads. In other words, as we tiling a read by units,
    // there is un-encoded regions at the both end of a read,
    // and we allocate memories for them.
    tips: Vec<DitchTip<'a>>,
    // Estimated copy number. If not able, None.
    copy_number: Option<usize>,
}

impl<'a> std::fmt::Display for DitchNode<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let cp = match self.copy_number {
            Some(res) => res as i64,
            None => -1,
        };
        writeln!(
            f,
            "---------{}:{}:{}---------",
            self.node.0, self.node.1, cp
        )?;
        writeln!(f, "Seq:{}", String::from_utf8_lossy(&self.seq))?;
        let lines: Vec<_> = self.edges.iter().map(|e| format!("{}", e)).collect();
        write!(f, "Edges:\n{}", lines.join("\n"))
    }
}

impl<'a> DitchNode<'a> {
    fn new(node: Node, seq: Vec<u8>) -> Self {
        Self {
            node,
            seq,
            edges: vec![],
            tips: vec![],
            occ: 0,
            copy_number: None,
        }
    }
    fn seq(&self) -> &[u8] {
        &self.seq
    }
    fn seq_as_string(&self) -> String {
        String::from_utf8(self.seq.clone()).unwrap()
    }
}

#[derive(Debug, Clone)]
pub struct DitchEdge {
    from: Node,
    to: Node,
    from_position: Position,
    to_position: Position,
    seq: EdgeLabel,
    occ: usize,
    // Estimated copy number. None if not available.
    copy_number: Option<usize>,
}

impl std::fmt::Display for DitchEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({:?}:{})->", self.from, self.from_position)?;
        write!(f, "({:?}:{})", self.to, self.to_position)?;
        write!(f, "{}:{}", self.occ, self.seq)
    }
}

impl std::cmp::PartialEq for DitchEdge {
    fn eq(&self, other: &DitchEdge) -> bool {
        self.from == other.from
            && self.to == other.to
            && self.from_position == other.from_position
            && self.to_position == other.to_position
    }
}

impl std::cmp::Eq for DitchEdge {}

impl DitchEdge {
    fn new(
        from: Node,
        from_position: Position,
        to: Node,
        to_position: Position,
        seq: EdgeLabel,
    ) -> Self {
        Self {
            from,
            to,
            from_position,
            to_position,
            occ: 0,
            seq,
            copy_number: None,
        }
    }
    // Reverse this edge.
    fn reverse(&self) -> Self {
        Self {
            from: self.to,
            from_position: self.to_position,
            to: self.from,
            to_position: self.from_position,
            seq: self.seq.reverse(),
            occ: 0,
            copy_number: self.copy_number,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DitchTip<'a> {
    seq: &'a [u8],
    position: Position,
    in_direction: bool,
}

impl<'a> DitchTip<'a> {
    fn new(seq: &'a [u8], position: Position, in_direction: bool) -> Self {
        Self {
            seq,
            position,
            in_direction,
        }
    }
}

fn revcmp_str(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
}

// I use this enumeration to
// tag a contig to nodes of a ditch graph.
// Sometimes, a node of a ditch graph
// has both end and start of the contig to
// each tip (head and tail). Thus,
// I explicitly express the position of each tags.
// A contig named `String`, the length of which is `usize`,
// ends at Position, or starts at `Position`, or both.
#[derive(Debug, Clone, Eq, PartialEq)]
enum ContigTag {
    Start(String, Position, usize),
    End(String, Position, usize),
    // Start position, end position.
    Both(String, Position, Position, usize),
    // None,
}

/// A summary of a contig. It tells us
/// the unit, the cluster, and the direction
/// it spelled.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ContigSummary {
    /// The ID of the focal contig.
    pub id: String,
    pub summary: Vec<ContigElement>,
}

impl ContigSummary {
    fn new(id: &str, summary: &[ContigElement]) -> Self {
        Self {
            id: id.to_string(),
            summary: summary.to_vec(),
        }
    }
}

impl std::fmt::Display for ContigSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line: Vec<_> = self
            .summary
            .iter()
            .map(|n| format!("{}-{}", n.unit, n.cluster))
            .collect();
        write!(f, "{}\t{}", self.id, line.join(":"))
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct ContigElement {
    pub unit: u64,
    pub cluster: u64,
    pub strand: bool,
    /// The "coverage."
    pub occ: usize,
    pub copy_number: Option<usize>,
}

// Take consensus on each node and each edge.
#[derive(Debug, Clone)]
enum EdgeLabel {
    // we should remove `-offset`bp from the end of the sequence.
    Ovlp(i64),
    Seq(Vec<u8>),
}

impl std::fmt::Display for EdgeLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            EdgeLabel::Ovlp(l) => write!(f, "{}", l),
            EdgeLabel::Seq(seq) => write!(f, "{}", String::from_utf8_lossy(seq)),
        }
    }
}

impl EdgeLabel {
    fn reverse(&self) -> Self {
        match self {
            EdgeLabel::Seq(seq) => EdgeLabel::Seq(bio_utils::revcmp(seq)),
            EdgeLabel::Ovlp(l) => EdgeLabel::Ovlp(*l),
        }
    }
    // Return the length of the edge. If the returned value is negative,
    // the two nodes connected by this edge is overlapping.
    pub fn len(&self) -> i64 {
        match self {
            EdgeLabel::Ovlp(l) => *l,
            EdgeLabel::Seq(s) => s.len() as i64,
        }
    }
}

fn take_consensus(
    reads: &[&EncodedRead],
    units: Option<&[Unit]>,
    c: &AssembleConfig,
) -> (NodeConsensus, EdgeConsensus) {
    let mut nodes: HashMap<_, Vec<_>> = HashMap::new();
    let mut edges: HashMap<_, Vec<_>> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            let tuple = (node.unit, node.cluster);
            nodes.entry(tuple).or_default().push(node);
        }
        for (i, edge) in read.edges.iter().enumerate() {
            let (from_node, to_node) = (&read.nodes[i], &read.nodes[i + 1]);
            let from = (from_node.unit, from_node.cluster);
            use Position::*;
            let from_pos = if from_node.is_forward { Tail } else { Head };
            let to = (to_node.unit, to_node.cluster);
            let to_pos = if to_node.is_forward { Head } else { Tail };
            let (tuple, value) = if from.0 < to.0 {
                (((from, from_pos), (to, to_pos)), (edge, true))
            } else {
                (((to, to_pos), (from, from_pos)), (edge, false))
            };
            edges.entry(tuple).or_default().push(value);
        }
    }
    let units: HashMap<_, _> = match units {
        Some(units) => units.iter().map(|u| (u.id, u.seq())).collect(),
        None => HashMap::new(),
    };
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let config = kiley::PolishConfig::with_model(100, 0, 0, 0, 0, hmm);
    let node_consensus: HashMap<_, _> = nodes
        .into_par_iter()
        .map(|(key, value)| {
            let draft = match units.get(&(key.0)) {
                Some(res) => res.to_vec(),
                None => value.iter().nth(0).map(|x| x.seq().to_vec()).unwrap(),
            };
            let cons = if c.to_polish {
                let seqs: Vec<_> = value.iter().map(|x| x.seq()).collect();
                kiley::polish_chunk_by_parts(&draft, &seqs, &config)
            } else {
                draft
            };
            (key, cons)
        })
        .collect();
    debug!("Consed nodes");
    let edge_consensus: HashMap<_, _> = edges
        .into_par_iter()
        .map(|(key, value)| {
            let offsets: i64 = value.iter().map(|x| x.0.offset).sum();
            let offset = offsets / value.len() as i64;
            let cons = match (c.to_polish, 0 < offset) {
                (_, false) => EdgeLabel::Ovlp(offset),
                (true, _) => {
                    let seqs: Vec<_> = value
                        .iter()
                        .map(|(ed, is_forward)| match is_forward {
                            true => ed.label().to_vec(),
                            false => bio_utils::revcmp(ed.label()),
                        })
                        .filter(|e| !e.is_empty())
                        .collect();
                    // let lens: Vec<_> = value.iter().map(|(ed, _)| ed.label().len()).collect();
                    // debug!("{:?}", lens);
                    // TODO:Maybe this is not so good...
                    let consensus = if seqs.len() > 3 {
                        kiley::consensus(&seqs, 0, 0, 100)
                    } else {
                        seqs.iter().max_by_key(|x| x.len()).unwrap().to_vec()
                    };
                    EdgeLabel::Seq(consensus)
                }
                (false, _) => {
                    let (edge, is_forward) = value.iter().max_by_key(|x| x.0.label.len()).unwrap();
                    match *is_forward {
                        true => EdgeLabel::Seq(edge.label().to_vec()),
                        false => EdgeLabel::Seq(bio_utils::revcmp(edge.label())),
                    }
                }
            };
            (key, cons)
        })
        .collect();
    (node_consensus, edge_consensus)
}

impl<'a> DitchGraph<'a> {
    pub fn sanity_check(&self) -> bool {
        self.nodes.values().all(|node| {
            node.edges.iter().all(|edge| {
                let rev = edge.reverse();
                let count = self.nodes[&edge.to]
                    .edges
                    .iter()
                    .filter(|e| e == &&rev)
                    .count();
                count == 1
            })
        })
    }
    pub fn new(reads: &[&'a EncodedRead], units: Option<&[Unit]>, c: &AssembleConfig) -> Self {
        // Take a consensus of nodes and edges.
        let (nodes_seq, edge_seq) = take_consensus(reads, units, c);
        // Allocate all nodes.
        let mut nodes: HashMap<_, _> = nodes_seq
            .into_iter()
            .map(|(node, seq)| (node, DitchNode::new(node, seq)))
            .collect();
        for (((from, from_pos), (to, to_pos)), seq) in edge_seq.into_iter() {
            let edge = DitchEdge::new(from, from_pos, to, to_pos, seq);
            nodes.get_mut(&to).unwrap().edges.push(edge.reverse());
            nodes.get_mut(&from).unwrap().edges.push(edge);
        }
        let mut graph = Self { nodes };
        // Append edges.
        for read in reads.iter() {
            graph.append_read(read, c);
            graph.append_tip(read, c);
        }
        graph
    }
    fn append_tip(&mut self, read: &'a EncodedRead, _c: &AssembleConfig) -> Option<()> {
        use Position::*;
        if let Some(first) = read.nodes.first() {
            let direction = true;
            let position = if first.is_forward { Head } else { Tail };
            self.nodes
                .get_mut(&(first.unit, first.cluster))
                .unwrap()
                .tips
                .push(DitchTip::new(&read.leading_gap, position, direction));
        }
        if let Some(last) = read.nodes.last() {
            let direction = true;
            let position = if last.is_forward { Head } else { Tail };
            self.nodes
                .get_mut(&(last.unit, last.cluster))
                .unwrap()
                .tips
                .push(DitchTip::new(&read.trailing_gap, position, direction));
        }
        Some(())
    }
    // Add weight of the graph.
    fn append_read(&mut self, read: &'a EncodedRead, _c: &AssembleConfig) -> Option<()> {
        for node in read.nodes.iter() {
            self.nodes.get_mut(&(node.unit, node.cluster)).unwrap().occ += 1;
        }
        for pairs in read.nodes.windows(2) {
            let (from, to) = (&pairs[0], &pairs[1]);
            let from_pos = if from.is_forward { Tail } else { Head };
            let to_pos = if to.is_forward { Head } else { Tail };
            let from = (from.unit, from.cluster);
            let to = (to.unit, to.cluster);
            use Position::*;
            let mock_edge = DitchEdge::new(from, from_pos, to, to_pos, EdgeLabel::Ovlp(0));
            self.nodes
                .get_mut(&from)
                .unwrap()
                .edges
                .iter_mut()
                .find(|e| e == &&mock_edge)
                .unwrap()
                .occ += 1;
            let mock_edge = mock_edge.reverse();
            self.nodes
                .get_mut(&to)
                .unwrap()
                .edges
                .iter_mut()
                .find(|e| e == &&mock_edge)
                .unwrap()
                .occ += 1;
        }
        Some(())
    }
}

impl<'a> DitchGraph<'a> {
    // Return tuples of node and their position from which
    // we can start simple-path-reduction.
    // In other words, it returns the list of the position which
    // does not has any edge in the opposite position.
    // For example, if there is no edge from Position::Tail at the (2,0) node,
    // then, ((2,0), Position::Head) would be included.
    // Also, it returns the list of the node whose parent has two or more children.
    fn enumerate_candidates(&self) -> Vec<(Node, Position)> {
        let mut selected = vec![false; self.nodes.len()];
        let mut primary_candidates = vec![];
        for (idx, (key, node)) in self.nodes.iter().enumerate() {
            for position in vec![Position::Head, Position::Tail] {
                let num_of_edge = node
                    .edges
                    .iter()
                    .filter(|e| e.from_position == position)
                    .count();
                if num_of_edge != 1 {
                    selected[idx] = true;
                    primary_candidates.push((*key, position));
                }
            }
        }
        // Secondary candidates.
        for ((key, node), _) in self
            .nodes
            .iter()
            .zip(selected)
            .filter(|(_, selected)| !selected)
        {
            for position in vec![Position::Head, Position::Tail] {
                let grand_child = node
                    .edges
                    .iter()
                    .find(|e| e.from_position == position)
                    .map(|e| {
                        let count = self.nodes[&e.to]
                            .edges
                            .iter()
                            .filter(|f| f.from_position == e.to_position)
                            .count();
                        if count == 0 {
                            debug!("{:?},{}", node, position);
                            debug!("Edge:{:?}", e);
                            for edge in &self.nodes[&e.to].edges {
                                debug!("Child:{:?}", edge)
                            }
                            panic!()
                        }
                        count
                    })
                    .unwrap();
                if 1 < grand_child {
                    primary_candidates.push((*key, position))
                }
            }
        }
        primary_candidates
    }
    // Only (from, _, to, _) with from <= to edges would be in these vectors.
    fn partition_edges_by_simple_path(&self) -> (Vec<&DitchEdge>, Vec<&DitchEdge>) {
        let (mut in_simple_path, mut between_simple_path) = (vec![], vec![]);
        for node in self.nodes.values() {
            for pos in vec![Position::Head, Position::Tail] {
                let count = node.edges.iter().filter(|e| e.from_position == pos).count();
                // if this is bidirectionaly unique, it is a edge in a simple path.
                let is_bi_unique = (count == 1) && {
                    let edge = node.edges.iter().find(|e| e.from_position == pos).unwrap();
                    let mut to_edge = self.nodes[&edge.to].edges.iter();
                    to_edge.all(|e| e.from_position != edge.to_position || e.to == edge.from)
                };
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                let edges = edges.filter(|e| e.from <= e.to);
                for edge in edges {
                    if is_bi_unique {
                        in_simple_path.push(edge);
                    } else {
                        between_simple_path.push(edge);
                    }
                }
            }
        }
        (in_simple_path, between_simple_path)
    }
    /// Estimoate copy number of nodes and edges.
    pub fn copy_number_estimation(
        &self,
        naive_cov: f64,
        lens: &[usize],
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        // Create simple-path reduced graph.
        let keys = self.nodes.keys().cloned().enumerate();
        let node_index: HashMap<_, _> = keys.map(|(i, k)| (k, i)).collect();
        let (edges_in_simple_path, edges_between_simple_path) =
            self.partition_edges_by_simple_path();
        use crate::find_union::FindUnion;
        let mut fu = FindUnion::new(node_index.len());
        for edge in edges_in_simple_path.iter() {
            fu.unite(node_index[&edge.from], node_index[&edge.to]);
        }
        let cluster_index: HashMap<_, _> = (0..node_index.len())
            .filter(|&n| fu.find(n) == Some(n))
            .enumerate()
            .map(|(idx, key)| (key, idx))
            .collect();
        use super::copy_number::CoverageCalibrator;
        let calibrator = CoverageCalibrator::new(&lens);
        let unit_len_sum: usize = self.nodes.values().map(|n| n.seq().len()).sum();
        let cov = calibrator.calib_f64(naive_cov, unit_len_sum / self.nodes.len());
        debug!("COVERAGE\t{:.3}\t{:.3}", naive_cov, cov,);
        // if (from, true, _, _) that read is in plugs[from][1].
        let mut plugs: Vec<_> = vec![vec![]; cluster_index.len()];
        let edges: Vec<_> = edges_between_simple_path
            .iter()
            .map(|edge| {
                let from_index = cluster_index[&fu.find(node_index[&edge.from]).unwrap()];
                let to_index = cluster_index[&fu.find(node_index[&edge.to]).unwrap()];
                // Decide the "position" of these nodes.
                // There should not be more than two plugs on each simple path.
                let from_node = (edge.from, edge.from_position);
                let fp = match plugs[from_index].iter().position(|n| n == &from_node) {
                    Some(idx) => idx == 1,
                    None => {
                        plugs[from_index].push(from_node);
                        plugs[from_index].len() == 2
                    }
                };
                assert!(plugs[from_index].len() <= 2);
                let to_node = (edge.to, edge.to_position);
                let tp = match plugs[to_index].iter().position(|n| n == &to_node) {
                    Some(idx) => idx == 1,
                    None => {
                        plugs[to_index].push(to_node);
                        plugs[to_index].len() == 2
                    }
                };
                assert!(plugs[to_index].len() <= 2);
                let unit_len =
                    self.nodes[&edge.from].seq().len() + self.nodes[&edge.to].seq().len();
                let gap_len = (unit_len as i64 + edge.seq.len()) as usize;
                let weight = calibrator.calib(edge.occ, gap_len) / cov;
                // debug!("CALIB\t{:.1}\t{}\t{}\tEDGE", edge.occ, gap_len, weight);
                (from_index, fp, to_index, tp, weight)
            })
            .collect();
        let mut node_weights = vec![(0f64, 0usize); cluster_index.len()];
        for node in self.nodes.values() {
            let cluster = fu.find(node_index[&node.node]).unwrap();
            let simple_path_index = cluster_index[&cluster];
            node_weights[simple_path_index].1 += 1;
            let unit_len = node.seq().len();
            let weight = calibrator.calib(node.occ, unit_len) / cov;
            // debug!("CALIB\t{:.1}\t{}\t{}\tNODE", node.occ, unit_len, weight);
            node_weights[simple_path_index].0 += weight;
        }
        let nodes: Vec<_> = node_weights
            .iter()
            .map(|&(sum, len)| sum / len as f64)
            .collect();
        let (node_cp, edge_cp) = super::copy_number::estimate_copy_number(&nodes, &edges);
        for (&(sum, len), cp) in node_weights.iter().zip(node_cp.iter()) {
            let weight = sum / len as f64;
            debug!("NODE\t{}\t{:.1}\t{}\t{:.2}", cp, sum, len, weight);
        }
        for ((from, fplus, to, tplus, w), cp) in edges.iter().zip(edge_cp.iter()) {
            debug!(
                "EDGE\t{}\t{}\t{}\t{}\t{:.1}\t{}",
                from, fplus, to, tplus, w, cp
            );
        }
        let node_copy_number: HashMap<_, _> = self
            .nodes
            .values()
            .map(|node| {
                let cluster = fu.find(node_index[&node.node]).unwrap();
                (node.node, node_cp[cluster_index[&cluster]])
            })
            .collect();
        let mut edge_copy_number = HashMap::new();
        for (&edge_cp, &(from, fplus, to, tplus, _)) in edge_cp.iter().zip(edges.iter()) {
            let from = plugs[from][fplus as usize];
            let to = plugs[to][tplus as usize];
            edge_copy_number.insert((from, to), edge_cp);
            edge_copy_number.insert((to, from), edge_cp);
        }
        for edge in self.nodes.values().flat_map(|node| node.edges.iter()) {
            let from = cluster_index[&fu.find(node_index[&edge.from]).unwrap()];
            let to = cluster_index[&fu.find(node_index[&edge.to]).unwrap()];
            if from == to {
                let cp = node_cp[from];
                let from = (edge.from, edge.from_position);
                let to = (edge.to, edge.to_position);
                edge_copy_number.insert((from, to), cp);
            }
        }
        (node_copy_number, edge_copy_number)
    }
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) = self.copy_number_estimation(naive_cov, lens);
        let get_edge_cn = |e: &DitchEdge| {
            let edge = ((e.from, e.from_position), (e.to, e.to_position));
            edge_copy_number.get(&edge).copied()
        };
        for (key, node) in self.nodes.iter_mut() {
            node.copy_number = node_copy_number.get(key).cloned();
            for edge in node.edges.iter_mut() {
                edge.copy_number = get_edge_cn(edge);
            }
        }
    }
    /// Even though the edge/node is zero copy number, we do not remove it if the conditions below hold:
    /// 1. If it is an edge, and all edge from the same position is ZCP, and it is the heviest edge among them.
    /// 2. If it is a node, and it has un-removable edge.
    /// 3. If it is a node, and the bounding constraint does not hold.
    ///    In other words, it connected to the non-ZCP edge.
    /// 4. If it is an edge, and the occ is more than `thr * max_out_dgree`
    pub fn remove_zero_copy_elements(&mut self, naive_cov: f64, lens: &[usize], thr: f64) {
        // let (node_copy_number, edge_copy_number) = self.copy_number_estimation(naive_cov, lens);
        // let get_edge_cn = |e: &DitchEdge| {
        //     let edge = ((e.from, e.from_position), (e.to, e.to_position));
        //     edge_copy_number[&edge]
        // };
        // Check the node violating for right-left condition.
        self.assign_copy_number(naive_cov, lens);
        let unsound_nodes: HashSet<_> = self
            .nodes
            .iter()
            .filter_map(|(key, node)| {
                let (plus, minus) = node.edges.iter().fold((0, 0), |(plus, minus), edge| {
                    let copy_num = match edge.copy_number {
                        Some(cp) => cp,
                        None => return (plus, minus),
                    };
                    match edge.from_position {
                        Position::Head => (plus + copy_num, minus),
                        Position::Tail => (plus, minus + copy_num),
                    }
                });
                match (plus == 0, minus == 0) {
                    (true, true) => Some(key),
                    (true, false) | (false, true) => None,
                    _ if plus != minus => Some(key),
                    _ => None,
                }
            })
            .collect();
        // let unsound_nodes: HashSet<_> = {
        // let mut node_count: HashMap<_, i64> =
        //     self.nodes.keys().copied().map(|k| (k, 0)).collect();
        // The Head position is "outcomming" edge and the Tail edge is "incoming".
        //     for edge in self.nodes.values().flat_map(|node| node.edges.iter()) {
        //         let copy_number = get_edge_cn(edge) as i64;
        //         *node_count.get_mut(&edge.from).unwrap() += match edge.from_position {
        //             Position::Head => -copy_number,
        //             Position::Tail => copy_number,
        //         };
        //         *node_count.get_mut(&edge.to).unwrap() += match edge.to_position {
        //             Position::Head => -copy_number,
        //             Position::Tail => copy_number,
        //         };
        //     }
        //     node_count
        //         .into_iter()
        //         .filter_map(|(key, value)| (value != 0).then(|| key))
        //         .collect()
        // };
        debug!("UNSOUND\t{}\t{}", unsound_nodes.len(), self.nodes.len());
        // (from,from_position,to,to_position) and from.0 <= to.0.
        let format_edge = |e: &DitchEdge| {
            if e.from.0 <= e.to.0 {
                (e.from, e.from_position, e.to, e.to_position)
            } else {
                (e.to, e.to_position, e.from, e.from_position)
            }
        };
        use super::copy_number::CoverageCalibrator;
        let calibrator = CoverageCalibrator::new(&lens);
        let sum_unit_len: usize = self.nodes.values().map(|node| node.seq().len()).sum();
        let ave_unit_len = sum_unit_len / self.nodes.len();
        let calib_occ = |edge: &DitchEdge| {
            let len = (edge.seq.len() + 2 * ave_unit_len as i64) as usize;
            calibrator.calib(edge.occ, len)
        };
        let mut is_ok_to_remove = HashSet::new();
        let mut retain_edges = HashSet::new();
        for (key, node) in self.nodes.iter() {
            if unsound_nodes.contains(key) {
                retain_edges.extend(node.edges.iter().map(format_edge));
                continue;
            }
            for position in vec![Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == position);
                // It is ok to take max by this way, as if the max is 0 there's no edge.
                let max = edges.clone().map(calib_occ).fold(0f64, |x, y| x.max(y));
                // It automatically drop the heviest edge from the selection.
                for edge in edges {
                    if edge.copy_number.map(|cp| cp == 0).unwrap_or(false)
                        && edge.occ as f64 / max < thr
                    {
                        is_ok_to_remove.insert(format_edge(edge));
                    } else {
                        retain_edges.insert(format_edge(edge));
                    }
                }
            }
        }
        for e in retain_edges.iter() {
            is_ok_to_remove.remove(e);
        }
        let (mut max_removed_weight, mut removed) = (0, 0);
        for node in self.nodes.values_mut() {
            node.edges.retain(|e| {
                if is_ok_to_remove.contains(&format_edge(e)) {
                    max_removed_weight = max_removed_weight.max(e.occ);
                    removed += 1;
                }
                !is_ok_to_remove.contains(&format_edge(e))
            });
        }
        debug!("MOREVED\t{}\t{}", removed, max_removed_weight);
        // If the condition 2. and 3. hold, then the degree should be 0.
        self.nodes.retain(|_, node| {
            node.copy_number.map(|cp| cp != 0).unwrap_or(true) || !node.edges.is_empty()
        });
    }
    #[allow(dead_code)]
    /// Enumerate bridges.
    /// Each element would be true if the node is an articulation point or each edge is bridge.
    fn enumerate_bridge(&self) -> HashMap<Node, (bool, Vec<bool>)> {
        // Serialize graph.
        let node_index: HashMap<_, _> = self
            .nodes
            .keys()
            .cloned()
            .enumerate()
            .map(|(i, n)| (n, i))
            .collect();
        // The i-th node would be 2*i and 2*i+1 node, representing Head/Tail position.
        let num_nodes = 2 * node_index.len();
        let mut edges: Vec<_> = (0..num_nodes)
            .map(|i| match i % 2 == 0 {
                true => vec![i + 1],
                false => vec![i - 1],
            })
            .collect();
        // Note that the `self` is already an undirected graph.
        for edge in self.nodes.values().flat_map(|n| n.edges.iter()) {
            let from = node_index[&edge.from];
            let to = node_index[&edge.to];
            let (i, j) = match (edge.from_position, edge.to_position) {
                (Position::Head, Position::Head) => (2 * from, 2 * to),
                (Position::Head, Position::Tail) => (2 * from, 2 * to + 1),
                (Position::Tail, Position::Head) => (2 * from, 2 * to),
                (Position::Tail, Position::Tail) => (2 * from + 1, 2 * to + 1),
            };
            edges[i].push(j);
        }
        let bridges = get_bridge(num_nodes, &edges);
        self.nodes
            .iter()
            .map(|(node, n)| {
                let index = 2 * node_index[node];
                let is_artic = bridges[index].contains(&(index + 1));
                // Linear search, little bit slow.
                let are_bridge: Vec<_> = n
                    .edges
                    .iter()
                    .map(|edge| {
                        let from = index + (edge.from_position == Position::Tail) as usize;
                        let to = 2 * node_index[&edge.to]
                            + (edge.to_position == Position::Tail) as usize;
                        edges[from].contains(&to)
                    })
                    .collect();
                (node.clone(), (is_artic, are_bridge))
            })
            .collect()
    }
    /// Remove small tips.
    /// TODO: More sophisticated algorithm is indeeed needed.
    pub fn remove_tips(&mut self) {
        let sum: usize = self.nodes.values().map(|e| e.occ).sum();
        let mean = sum / self.nodes.len();
        let thr = mean / 8;
        debug!("Removing nodes less than {} occ.", thr);
        let to_remove: HashSet<_> = self
            .nodes
            .iter()
            .filter_map(|(&key, node)| (node.occ < thr).then(|| key))
            .collect();
        self.remove_nodes(&to_remove);
    }
    // Removing nodes and re-map all of the indices.
    fn remove_nodes(&mut self, to_remove: &HashSet<Node>) {
        self.nodes.retain(|key, _| !to_remove.contains(key));
        for node in self.nodes.values_mut() {
            node.edges.retain(|edge| !to_remove.contains(&edge.to));
        }
    }
    /// Resolve repetitive units, selecting edges, simplify graph.
    /// This workflow is implemented by "foci resolving" algorithm.
    pub fn resolve_tangles(&self, reads: &[&EncodedRead], config: &AssembleConfig) {
        let _foci = self.get_foci(reads, config);
        unimplemented!()
    }
    /// Return a hash map containing all foci, thresholded by `config` parameter.
    pub fn get_foci(
        &self,
        reads: &[&EncodedRead],
        config: &AssembleConfig,
    ) -> HashMap<(Node, Position), Focus> {
        let mut foci = HashMap::new();
        for node in self.nodes.values() {
            for &pos in &[Position::Head, Position::Tail] {
                if let Some(focus) = self.examine_focus(node, pos, &reads, config) {
                    foci.insert((node.node, pos), focus);
                }
            }
        }
        foci
    }
    /// Return the focus if given node has a focus. Otherwise, return None.
    pub fn examine_focus(
        &self,
        node: &DitchNode,
        pos: Position,
        reads: &[&EncodedRead],
        _config: &AssembleConfig,
    ) -> Option<Focus> {
        // TODO: param.
        let dist_nodes = self.enumerate_node_upto(node, pos, 5);
        let reads: Vec<_> = reads.iter().filter(|r| r.contains(node.node)).collect();
        let mut focus: Option<Focus> = None;
        for (dist, nodes) in dist_nodes
            .iter()
            .enumerate()
            .filter(|(_, ns)| !ns.is_empty())
        {
            let dist = dist + 1;
            debug!("{}\t{:?}", dist, nodes);
            let mut occs: Vec<_> = vec![0; nodes.len()];
            for read in reads.iter() {
                // Safe.
                let start = read
                    .nodes
                    .iter()
                    .position(|n| (n.unit, n.cluster) == node.node)
                    .unwrap();
                let check_node = match (read.nodes[start].is_forward, pos == Position::Tail) {
                    (true, true) | (false, false) if start + dist < read.nodes.len() => {
                        &read.nodes[start + dist]
                    }
                    (false, true) | (true, false) if dist <= start => &read.nodes[start - dist],
                    _ => continue,
                };
                let check_node = (check_node.unit, check_node.cluster);
                if let Some(hit_node) = nodes.iter().position(|n| n.0 == check_node) {
                    occs[hit_node] += 1;
                }
            }
            // TODO: Assume that the copy numeber is 1. Thus, each element in nodes would
            // happen in 1/nodes.len() probability.
            // Maybe we need to correct this setting, by introducing copy number at each nodes/edges.
            let ith_ln = vec![-(nodes.len() as f64).ln(); nodes.len()];
            let null_prob: f64 = occs
                .iter()
                .zip(ith_ln.iter())
                .map(|(&occ, ln)| occ as f64 * ln)
                .sum();
            // TODO: Parametrize here.
            let err_prob = 0.1;
            let choice_num = nodes.len() as f64;
            let correct_lk = ((1f64 - err_prob).powi(2) + err_prob / choice_num).ln();
            let error_lk = match nodes.len() == 1 {
                true => ((1f64 - err_prob) * err_prob + err_prob / choice_num).ln(),
                _ => {
                    let correct_to_error = (1.0 - err_prob) * err_prob * (choice_num - 1.0).recip();
                    let error_to_error = err_prob / choice_num;
                    (correct_to_error + error_to_error).ln()
                }
            };
            let (alt_prob, alt_idx): (f64, _) = (0..nodes.len())
                .map(|k| {
                    let lk: f64 = occs
                        .iter()
                        .map(|&x| x as f64)
                        .enumerate()
                        .map(|(i, occ)| occ * (if i == k { correct_lk } else { error_lk }))
                        .sum();
                    (lk, k)
                })
                .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
                .unwrap();
            let lk_ratio = alt_prob - null_prob;
            if focus.as_ref().map(|f| f.llr()).unwrap_or(0f64) < lk_ratio {
                let (to, to_pos) = nodes[alt_idx];
                focus = Some(Focus::new(node.node, pos, to, to_pos, dist, lk_ratio));
            }
        }
        focus
    }
    /// Return nodes achievable from the given node. The i-th vector is
    /// the nodes reachable from `node` in i+1 hop.
    pub fn enumerate_node_upto(
        &self,
        node: &DitchNode,
        pos: Position,
        radius: usize,
    ) -> Vec<Vec<(Node, Position)>> {
        let initial_nodes: Vec<_> = node
            .edges
            .iter()
            .filter(|edge| edge.from_position == pos)
            .map(|edge| (edge.to, edge.to_position))
            .collect();
        let mut nodes_at: Vec<Vec<_>> = vec![initial_nodes];
        for i in 0..radius {
            let mut next_nodes = vec![];
            for (node, p) in nodes_at[i].iter() {
                let next_edges = match self.nodes.get(node) {
                    Some(n) => &n.edges,
                    None => continue,
                };
                let next_edges = next_edges
                    .iter()
                    .filter(|e| e.from_position == !*p)
                    .map(|e| (e.from, e.from_position));
                next_nodes.extend(next_edges);
            }
            next_nodes.sort();
            next_nodes.dedup();
            nodes_at.push(next_nodes);
        }
        nodes_at
    }
    /// Remove edge from given reads.
    pub fn remove_edges_from_short_reads(&mut self, reads: &[&EncodedRead]) {
        for read in reads.iter() {
            for w in read.nodes.windows(2) {
                let from = (w[0].unit, w[0].cluster);
                use Position::*;
                let from_pos = if w[0].is_forward { Tail } else { Head };
                let to = (w[1].unit, w[1].cluster);
                let to_pos = if w[1].is_forward { Tail } else { Head };
                let mock_edge = DitchEdge::new(from, from_pos, to, to_pos, EdgeLabel::Ovlp(0));
                let forward_edge = self
                    .nodes
                    .get_mut(&from)
                    .and_then(|n| n.edges.iter_mut().find(|e| e == &&mock_edge));
                if let Some(forward_edge) = forward_edge {
                    forward_edge.occ -= 1;
                }
                let mock_edge = mock_edge.reverse();
                let reverse_edge = self
                    .nodes
                    .get_mut(&to)
                    .and_then(|n| n.edges.iter_mut().find(|e| e == &&mock_edge));
                if let Some(reverse_edge) = reverse_edge {
                    reverse_edge.occ -= 1;
                }
            }
        }
    }
    /// Remove lightweight edges with occurence less than `thr`.
    /// Note that this function would never broke a connected component into two.
    pub fn remove_lightweight_edges(&mut self, thr: usize) {
        let mut removed_edges: HashMap<_, Vec<_>> =
            self.nodes.keys().map(|&k| (k, vec![])).collect();
        for (from, node) in self.nodes.iter() {
            for &pos in &[Position::Head, Position::Tail] {
                let to_cands = node
                    .edges
                    .iter()
                    .filter(|e| e.from_position == pos)
                    .map(|e| (e.to_position, e.to))
                    .count();
                if to_cands > 1 {
                    let max_occ = node.edges.iter().map(|x| x.occ).max().unwrap();
                    for e in node
                        .edges
                        .iter()
                        .filter(|e| e.from_position == pos && e.occ <= thr && e.occ < max_occ)
                    {
                        let to = e.to;
                        let t_pos = e.to_position;
                        removed_edges
                            .entry(*from)
                            .or_default()
                            .push((pos, to, t_pos));
                        removed_edges
                            .entry(to)
                            .or_default()
                            .push((t_pos, *from, pos));
                    }
                }
            }
        }
        self.nodes.iter_mut().for_each(|(key, node)| {
            let to_be_removed = &removed_edges[key];
            node.edges.retain(|x| {
                let probe = (x.from_position, x.to, x.to_position);
                !to_be_removed.contains(&probe)
            })
        });
    }
    /// Remove small connected components with the size less than `thr`.
    pub fn remove_small_component(&mut self, thr: usize) {
        let mut to_remove = HashSet::new();
        let mut arrived = HashSet::new();
        let mut cluster_num = 0;
        for node in self.nodes.keys() {
            if arrived.contains(node) {
                continue;
            }
            let mut actives = HashSet::new();
            let mut stack = vec![*node];
            while !stack.is_empty() {
                let last = stack.pop().unwrap();
                actives.insert(last);
                arrived.insert(last);
                for to in self.nodes[&last].edges.iter().map(|x| x.to) {
                    if !arrived.contains(&to) {
                        stack.push(to);
                    }
                }
            }
            debug!("CLUSTER\t{}\t{}", cluster_num, actives.len());
            cluster_num += 1;
            if actives.len() < thr {
                to_remove.extend(actives);
            }
        }
        self.remove_nodes(&to_remove);
    }
    /// Zip up bubble into a straight line.
    pub fn collapse_bubble(&mut self, c: &AssembleConfig) {
        let mut to_remove = HashSet::new();
        let mut queue = std::collections::VecDeque::new();
        for &key in self.nodes.keys() {
            queue.push_back((key, Position::Head));
            queue.push_back((key, Position::Tail));
        }
        while let Some((key, position)) = queue.pop_front() {
            let edges: Vec<_> = self.nodes[&key]
                .edges
                .iter()
                .filter(|e| e.from_position == position)
                .collect();
            if edges.len() <= 1 {
                continue;
            }
            let pos = edges[0].to_position;
            let unit = edges[0].to.0;
            let is_the_same_unit = edges.iter().all(|e| e.to_position == pos && e.to.0 == unit);
            let key_is_the_unique_parent = edges.iter().all(|e| {
                self.nodes[&e.to]
                    .edges
                    .iter()
                    .all(|f| (f.from_position != pos) | (f.to == key))
            });
            if is_the_same_unit && key_is_the_unique_parent {
                let (new_terminal, removed_nodes) = self.collapse_bubble_from(key, position, c);
                to_remove.extend(removed_nodes);
                queue.push_back(new_terminal)
            }
        }
        self.remove_nodes(&to_remove);
    }
    // Collapse bubble from the given root.
    // It is garanteed that the (root, position) node has at least two children,
    // all of their parent is the (root, position) itself.
    // Like,
    //
    //           /-------((to, _),_)----
    //  ----((u,c),p)
    //           \-------((to, _),_)---
    //
    // This.
    fn collapse_bubble_from(
        &mut self,
        root: Node,
        position: Position,
        _c: &AssembleConfig,
    ) -> ((Node, Position), Vec<Node>) {
        assert!(self.sanity_check());
        // Check the collapsing condition.
        let edges: Vec<_> = self.nodes[&root]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .collect();
        // Check.
        assert!(edges.len() > 1);
        // TODO: Maybe we need to consensus these bubbles.
        let (first_node, first_pos, seq, edgelabel) = {
            let first = edges.iter().max_by_key(|x| x.occ).unwrap();
            // let first = edges.first().unwrap();
            let first_node = first.to;
            let seq = self.nodes[&first.to].seq().to_vec();
            let edgelabel = first.seq.clone();
            assert_eq!(first.from, root);
            assert_eq!(first.from_position, position);
            (first_node, first.to_position, seq, edgelabel)
        };
        assert!(edges.iter().all(|e| e.to_position == first_pos));
        assert!(edges.iter().all(|e| e.to.0 == first_node.0));
        let merged_nodes: Vec<_> = edges.iter().map(|e| e.to).collect();
        let mut new_node = DitchNode::new(first_node, seq);
        let edge_to_root = DitchEdge::new(first_node, first_pos, root, position, edgelabel);
        new_node.edges.push(edge_to_root);
        let node_otherside: HashSet<_> = merged_nodes
            .iter()
            .flat_map(|node| self.nodes[node].edges.iter())
            .filter_map(|e| (e.from_position != first_pos).then(|| (e.to, e.to_position)))
            .collect();
        for node in merged_nodes.iter() {
            new_node.occ += self.nodes[node].occ;
            if let Some(node) = self.nodes.get_mut(node) {
                while let Some(tip) = node.tips.pop() {
                    new_node.tips.push(tip)
                }
            }
            while let Some(mut edge) = self.nodes.get_mut(&node).unwrap().edges.pop() {
                if let Some(existing) = new_node
                    .edges
                    .iter_mut()
                    .find(|e| e.to == edge.to && e.to_position == edge.to_position)
                {
                    existing.occ += edge.occ;
                } else {
                    edge.from = new_node.node;
                    new_node.edges.push(edge);
                }
            }
        }
        self.nodes
            .get_mut(&root)
            .unwrap()
            .edges
            .retain(|e| (e.from_position != position));
        let root_to_edge = match new_node
            .edges
            .iter()
            .find(|e| e.to == root && e.to_position == position)
        {
            Some(res) => res.reverse(),
            None => panic!(),
        };
        self.nodes.get_mut(&root).unwrap().edges.push(root_to_edge);
        self.nodes.insert(first_node, new_node);
        let num_edge = self.nodes[&root]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count();
        assert_eq!(num_edge, 1);
        for (other, other_position) in node_otherside {
            let (removed_occ, label) = self.nodes[&other]
                .edges
                .iter()
                .filter(|e| other_position == e.from_position && merged_nodes.contains(&e.to))
                .fold((0, EdgeLabel::Ovlp(0)), |(cum, _), e| {
                    (cum + e.occ, e.seq.clone())
                });
            self.nodes
                .get_mut(&other)
                .unwrap()
                .edges
                .retain(|e| !(other_position == e.from_position && merged_nodes.contains(&e.to)));
            let mut new_edge = DitchEdge::new(other, other_position, first_node, !first_pos, label);
            new_edge.occ += removed_occ;
            self.nodes.get_mut(&other).unwrap().edges.push(new_edge);
        }
        let mut merged_nodes = merged_nodes;
        merged_nodes.retain(|&x| x != first_node);
        ((first_node, !first_pos), merged_nodes)
    }
}

#[allow(dead_code)]
fn get_bridge(nodes: usize, edges: &[Vec<usize>]) -> Vec<Vec<usize>> {
    warn!("Have you tested this code? I think not.");
    let mut is_arrived = vec![false; nodes];
    let mut bridges = vec![vec![]; nodes];
    let mut order = vec![-1; nodes];
    let mut low = vec![-1; nodes];
    let mut parent = vec![0; nodes];
    let mut count = 0;
    for i in 0..nodes {
        if is_arrived[i] {
            continue;
        }
        let mut stack = vec![i];
        'dfs: while !stack.is_empty() {
            let last = *stack.last().unwrap();
            if !is_arrived[last] {
                is_arrived[last] = true;
                order[last] = count as i64;
                count += 1;
            }
            for &to in edges[last].iter() {
                if !is_arrived[to] {
                    parent[to] = last;
                    stack.push(to);
                    continue 'dfs;
                }
            }
            // Postorder, we have arrived all the nodes below this node.
            let last = stack.pop().unwrap();
            for &to in edges[last].iter() {
                if parent[to] == last {
                    low[last] = low[last].min(low[to]);
                    if order[last] < low[to] {
                        bridges[last].push(to);
                        bridges[to].push(last);
                    }
                } else {
                    // This is a back edge.
                    low[last] = low[last].min(order[to]);
                }
            }
        }
    }
    bridges
}
