use super::copy_number;
use super::AssembleConfig;
// use crate::find_union;
use definitions::{EncodedRead, Unit};
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
// const DRAFT_COVERAGE: usize = 20;
// const DRAFT_REP_NUM: usize = 7;

#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash)]
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

/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
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
}

impl<'a> std::fmt::Display for DitchNode<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "---------{}:{}---------", self.node.0, self.node.1)?;
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
    /// Reduce simple path of this graph and returns the edges and nodes of the reduced graph..
    pub fn spell(
        &self,
        c: &AssembleConfig,
        cl: usize,
    ) -> (
        Vec<gfa::Segment>,
        Vec<(gfa::Edge, Vec<gfa::SamTag>)>,
        gfa::Group,
        Vec<ContigSummary>,
    ) {
        let mut arrived = HashSet::new();
        let mut sids: HashMap<_, _> = HashMap::new();
        let (mut g_segs, mut g_edges, mut summaries) = (vec![], vec![], vec![]);
        let candidates = self.enumerate_candidates();
        for (node, p) in candidates {
            if arrived.contains(&node) {
                continue;
            }
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges, summary) =
                self.traverse_from(&mut arrived, &mut sids, node, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
        }
        for key in self.nodes.keys() {
            if arrived.contains(key) {
                continue;
            }
            let p = Position::Head;
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges, summary) =
                self.traverse_from(&mut arrived, &mut sids, *key, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
        }
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.0.eid.clone()))
            .collect();
        let uid = Some(format!("group-{}", cl));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group, summaries)
    }
    pub fn copy_number_estimation(
        &self,
        cv: f64,
        len: f64,
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        // Create simple-path reduced graph.
        let (reduced_nodes, nodes_in_sp) = self.simple_path_reduction();
        let is_tip: Vec<_> = nodes_in_sp
            .iter()
            .map(|nodes| self.outdeg_of_nodes(&nodes) < 2)
            .collect();
        let node_coverage: Vec<_> = nodes_in_sp
            .iter()
            .map(|nodes| self.coverage_of_nodes(nodes) / cv)
            .collect();
        let reduced_edges = self.edges_between_simple_path(&reduced_nodes, cv, len);
        let (nodes_cp, edge_cp) =
            copy_number::estimate_copy_number(&node_coverage, &is_tip, &reduced_edges);
        let nodes_cp: HashMap<_, _> = nodes_cp
            .iter()
            .zip(nodes_in_sp)
            .flat_map(|(&cp, nodes)| nodes.iter().map(|&node| (node, cp)).collect::<Vec<_>>())
            .collect();
        let mut edge_cp: HashMap<DitEdge, _> = edge_cp
            .iter()
            .zip(reduced_edges)
            .map(|(&cp, (from, from_first, to, to_first, _))| {
                let from = match reduced_nodes[from] {
                    (node, pos, _, _) if from_first => (node, pos),
                    (_, _, node, pos) => (node, pos),
                };
                let to = match reduced_nodes[to] {
                    (node, pos, _, _) if to_first => (node, pos),
                    (_, _, node, pos) => (node, pos),
                };
                ((from, to), cp)
            })
            .collect();
        for node in self.nodes.values() {
            for edge in node.edges.iter() {
                let from = (edge.from, edge.from_position);
                let to = (edge.to, edge.to_position);
                if edge_cp.contains_key(&(from, to)) {
                    // It is OK to insert either copy number.
                    let copy_number = nodes_cp[&edge.from];
                    edge_cp.insert((from, to), copy_number);
                }
            }
        }
        (nodes_cp, edge_cp)
    }
    fn simple_path_reduction(&self) -> (Vec<(Node, Position, Node, Position)>, Vec<Vec<Node>>) {
        unimplemented!()
    }
    fn outdeg_of_nodes(&self, nodes: &[Node]) -> usize {
        unimplemented!()
    }
    fn coverage_of_nodes(&self, nodes: &[Node]) -> f64 {
        unimplemented!()
    }
    fn edges_between_simple_path(
        &self,
        nodes: &[(Node, Position, Node, Position)],
        cv: f64,
        len: f64,
    ) -> Vec<(usize, bool, usize, bool, f64)> {
        unimplemented!()
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
    /// Resolve repetitive units. (TODO.)
    pub fn resolve_repeats(&self) {}
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
        // use find_union::FindUnion;
        // let mut cluster = FindUnion::new(self.nodes.len());
        // for node in self.nodes.values() {
        //     for edge in node.edges.iter() {
        //         cluster.unite(edge.from, edge.to);
        //     }
        // }
        // for node in 0..self.nodes.len() {
        //     if cluster.find(node).unwrap() != node {
        //         continue;
        //     }
        //     let count = (0..self.nodes.len())
        //         .filter(|&x| cluster.find(x).unwrap() == node)
        //         .count();
        //     if count < thr {
        //         for (component, to_remove) in to_remove.iter_mut().enumerate() {
        //             if cluster.find(component).unwrap() == node {
        //                 *to_remove = true;
        //             }
        //         }
        //     }
        // }
        // self.remove_nodes(&to_remove);
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
    // Traverse from the given `start` node of `start_position` Position.
    fn traverse_from(
        &self,
        arrived: &mut HashSet<Node>,
        sids: &mut HashMap<Node, ContigTag>,
        start: Node,
        start_position: Position,
        seqname: String,
        _c: &AssembleConfig,
    ) -> (
        gfa::Segment,
        Vec<(gfa::Edge, Vec<gfa::SamTag>)>,
        ContigSummary,
    ) {
        // Find edges.
        let mut edges: Vec<_> = self.nodes[&start]
            .edges
            .iter()
            .filter(|e| e.from_position == start_position)
            .filter_map(|e| {
                assert!(e.from == start);
                let (sid2, beg2) = match sids.get(&e.to)? {
                    &ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    &ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    &ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    &ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(0, false);
                let a = None;
                let edge = gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a);
                let len = match &e.seq {
                    EdgeLabel::Ovlp(_) => 0,
                    EdgeLabel::Seq(l) => l.len(),
                };
                let samtag = vec![
                    gfa::SamTag::new(format!("cv:i:{}", e.occ)),
                    gfa::SamTag::new(format!("ln:i:{}", len)),
                ];
                Some((edge, samtag))
            })
            .collect();
        let (mut node, mut position) = (start, start_position);
        let mut seq = self.initial_sequence(start, start_position);
        // Start traveresing.
        let mut unit_names = vec![];
        loop {
            arrived.insert(node);
            // Move forward.
            let cons = match position {
                Position::Head => self.nodes[&node].seq_as_string(),
                Position::Tail => revcmp_str(&self.nodes[&node].seq_as_string()),
            };
            {
                let elm = ContigElement {
                    unit: self.nodes[&node].node.0,
                    cluster: self.nodes[&node].node.1,
                    strand: position == Position::Head,
                    occ: self.nodes[&node].occ,
                };
                unit_names.push(elm);
            }
            seq += &cons;
            position = !position;
            // Check.
            let num_edges = self.nodes[&node]
                .edges
                .iter()
                .filter(|e| e.from_position == position)
                .count();
            if num_edges == 0 || num_edges > 1 {
                break;
            }
            // There is only one child.
            let selected_edge = self.nodes[&node]
                .edges
                .iter()
                .find(|e| e.from_position == position)
                .unwrap();
            assert_eq!(selected_edge.from, node);
            assert_eq!(selected_edge.from_position, position);
            // Succeed along with the edge.
            match &selected_edge.seq {
                EdgeLabel::Ovlp(l) => assert!((0..(-l)).all(|_| seq.pop().is_some())),
                EdgeLabel::Seq(label) => seq.extend(label.iter().map(|&x| x as char)),
            };
            let (next, next_position) = (selected_edge.to, selected_edge.to_position);
            // Check the number of child.
            let num_children = self.nodes[&next]
                .edges
                .iter()
                .filter(|e| e.from_position == next_position)
                .count();
            if num_children >= 2 || arrived.contains(&next) {
                break;
            }
            // Jump to the next node.
            assert!(num_children == 1);
            node = next;
            position = next_position;
        }
        seq += &self.trailing_sequence(node, position);
        let summary = ContigSummary::new(&seqname, &unit_names);
        // Register start and tail node.
        // This if statement is no needed?
        if start == node {
            let tag = ContigTag::Both(seqname.clone(), start_position, position, seq.len());
            sids.insert(node, tag);
        } else {
            let tag = ContigTag::Start(seqname.clone(), start_position, seq.len());
            sids.insert(start, tag);
            let tag = ContigTag::End(seqname.clone(), position, seq.len());
            sids.insert(node, tag);
        }
        let seg = gfa::Segment::from(seqname.clone(), seq.len(), Some(seq.clone()));
        // Add gfa edges.
        let tail_edges = self.nodes[&node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .filter_map(|e| {
                assert!(e.from == node);
                let (sid2, beg2) = match sids.get(&e.to)? {
                    &ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    &ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    &ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    &ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(seq.len(), true);
                let a = None;
                let edge = gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a);
                let len = match &e.seq {
                    EdgeLabel::Ovlp(_) => 0,
                    EdgeLabel::Seq(l) => l.len(),
                };
                let samtag = vec![
                    gfa::SamTag::new(format!("cv:i:{}", e.occ)),
                    gfa::SamTag::new(format!("ln:i:{}", len)),
                ];
                Some((edge, samtag))
            });
        edges.extend(tail_edges);
        (seg, edges, summary)
    }
    fn initial_sequence(&self, node: Node, position: Position) -> String {
        let num_edges = self.nodes[&node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count();
        if num_edges > 0 {
            String::new()
        } else {
            self.nodes[&node]
                .tips
                .iter()
                .filter(|tip| tip.position == position && !tip.seq.is_empty())
                .map(|tip| match tip.in_direction {
                    true => tip.seq.to_vec(),
                    false => bio_utils::revcmp(tip.seq),
                })
                .max_by_key(|x| x.len())
                .and_then(|x| String::from_utf8(x).ok())
                .unwrap_or_else(String::new)
        }
    }
    fn trailing_sequence(&self, node: Node, position: Position) -> String {
        let num_edges = self.nodes[&node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count();
        if num_edges > 0 {
            String::new()
        } else {
            self.nodes[&node]
                .tips
                .iter()
                .filter(|tip| tip.position == position && !tip.seq.is_empty())
                .map(|tip| match tip.in_direction {
                    true => bio_utils::revcmp(tip.seq),
                    false => tip.seq.to_vec(),
                })
                .max_by_key(|x| x.len())
                .and_then(|x| String::from_utf8(x).ok())
                .unwrap_or_else(String::new)
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(dead_code)]
    use definitions::{Edge, Node};
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128StarStar;
    // Raw data, expected result.
    type DataPack = (Vec<EncodedRead>, Vec<String>);
    fn gen_node(unit: u64, cluster: u64, is_forward: bool, seq: String) -> Node {
        let mut node = Node::default();
        node.unit = unit;
        node.cluster = cluster;
        node.is_forward = is_forward;
        node.seq = seq;
        node
    }
    fn gen_edge(from: u64, to: u64, offset: i64, label: String) -> Edge {
        Edge {
            from,
            to,
            offset,
            label,
        }
    }
    const BASES: &'static [char] = &['A', 'C', 'G', 'T'];
    use rand::seq::SliceRandom;
    fn revcmp(seq: &[char]) -> Vec<char> {
        seq.iter()
            .rev()
            .map(|&x| match x {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                _ => unreachable!(),
            })
            .collect()
    }
    fn gen_mock1() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        let seq: Vec<_> = (0..30)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let mut read = EncodedRead::default();
        read.nodes
            .push(gen_node(0, 1, true, seq[..10].iter().collect()));
        read.nodes
            .push(gen_node(1, 1, true, seq[11..26].iter().collect()));
        read.nodes
            .push(gen_node(2, 1, false, revcmp(&seq[20..27]).iter().collect()));
        read.nodes
            .push(gen_node(3, 1, true, seq[29..30].iter().collect()));
        read.edges
            .push(gen_edge(0, 1, 1, seq[10..11].iter().collect()));
        read.edges.push(gen_edge(1, 2, -6, String::new()));
        read.edges
            .push(gen_edge(2, 3, 2, seq[27..29].iter().collect()));
        let reads = vec![read];
        (reads, vec![seq.iter().collect()])
    }
    fn gen_mock_large() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        let seq_arms: Vec<Vec<_>> = (0..4)
            .map(|_| {
                (0..550)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let seq_center: Vec<_> = (0..100)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..5 {
            let seq: String = seq_arms[0][unit * 110..unit * 110 + 100].iter().collect();
            read0.nodes.push(gen_node(unit as u64, 0, true, seq));
            let seq: String = seq_arms[1][unit * 110..unit * 110 + 100].iter().collect();
            read1.nodes.push(gen_node(unit as u64, 1, true, seq));
            let lab: String = seq_arms[0][unit * 110 + 100..(unit + 1) * 110]
                .iter()
                .collect();
            read0
                .edges
                .push(gen_edge(unit as u64, unit as u64 + 1, 10, lab));
            let lab: String = seq_arms[1][unit * 110 + 100..(unit + 1) * 110]
                .iter()
                .collect();
            read1
                .edges
                .push(gen_edge(unit as u64, unit as u64 + 1, 10, lab));
        }
        read0
            .nodes
            .push(gen_node(5, 0, true, seq_center.iter().collect()));
        read1
            .nodes
            .push(gen_node(5, 0, true, seq_center.iter().collect()));
        let (mut read2, mut read3) = (read.clone(), read.clone());
        read2
            .nodes
            .push(gen_node(5, 0, true, seq_center.iter().collect()));
        read3
            .nodes
            .push(gen_node(5, 0, true, seq_center.iter().collect()));
        for unit in 0..5 {
            let lab: String = seq_arms[2][unit * 110..unit * 110 + 10].iter().collect();
            read2
                .edges
                .push(gen_edge(unit as u64 + 5, unit as u64 + 6, 10, lab));
            let lab: String = seq_arms[3][unit * 110..unit * 110 + 10].iter().collect();
            read3
                .edges
                .push(gen_edge(unit as u64 + 5, unit as u64 + 6, 10, lab));
            let seq: String = seq_arms[2][10 + unit * 110..(unit + 1) * 110]
                .iter()
                .collect();
            read2.nodes.push(gen_node(unit as u64 + 6, 0, true, seq));
            let seq: String = seq_arms[3][10 + unit * 110..(unit + 1) * 110]
                .iter()
                .collect();
            read3.nodes.push(gen_node(unit as u64 + 6, 1, true, seq));
        }
        let answer: Vec<_> = seq_arms
            .into_iter()
            .chain(std::iter::once(seq_center))
            .map(|cs| cs.into_iter().collect::<String>())
            .collect();
        (vec![read0, read1, read2, read3], answer)
    }
    fn gen_mock_large_2() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(11231312);
        let seq: Vec<Vec<_>> = (0..=8)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        let edg = String::new();
        for unit in 0..7 {
            let (subseq, dir): (String, _) = if [0, 1, 4, 5, 7].contains(&unit) {
                (seq[unit].iter().collect(), true)
            } else {
                (revcmp(&seq[unit]).iter().collect(), false)
            };
            let (cl0, cl1) = if [0, 1, 6, 7].contains(&unit) {
                (0, 0)
            } else {
                (0, 1)
            };
            let unit = unit as u64;
            read0.nodes.push(gen_node(unit, cl0, dir, subseq.clone()));
            read1.nodes.push(gen_node(unit, cl1, dir, subseq));
            read0.edges.push(gen_edge(unit, unit + 1, 0, edg.clone()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, edg.clone()));
        }
        read0
            .nodes
            .push(gen_node(7, 0, true, seq[7].iter().collect()));
        read1
            .nodes
            .push(gen_node(7, 0, true, seq[7].iter().collect()));
        let (read0, read1) = (read0, read1);
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (1..=7).rev() {
            let (subseq, dir): (String, _) = if [0, 1, 4, 5, 7].contains(&unit) {
                (revcmp(&seq[unit]).iter().collect(), false)
            } else {
                (seq[unit].iter().collect(), true)
            };
            let (cl0, cl1) = if [0, 1, 6, 7].contains(&unit) {
                (0, 0)
            } else {
                (0, 1)
            };
            let unit = unit as u64;
            read2.nodes.push(gen_node(unit, cl0, dir, subseq.clone()));
            read3.nodes.push(gen_node(unit, cl1, dir, subseq));
            read2.edges.push(gen_edge(unit - 1, unit, 0, edg.clone()));
            read3.edges.push(gen_edge(unit - 1, unit, 0, edg.clone()));
        }
        let seq0: String = revcmp(&seq[0]).iter().collect();
        read2.nodes.push(gen_node(0, 0, false, seq0.clone()));
        read3.nodes.push(gen_node(0, 0, false, seq0.clone()));
        let reads = vec![read0, read1, read2, read3];
        let middle = (2..=5).fold(String::new(), |mut x, y| {
            x.extend(&seq[y]);
            x
        });
        let answer: Vec<String> = vec![
            seq[0].iter().chain(seq[1].iter()).collect(),
            middle.clone(),
            middle,
            seq[6].iter().chain(seq[7].iter()).collect(),
        ];
        (reads, answer)
    }
    fn gen_circle() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(111);
        let seq: Vec<_> = (0..880)
            .filter_map(|_| BASES.choose(&mut rng))
            .copied()
            .collect();
        let mut read = EncodedRead::default();
        for unit in 0..8 {
            let subseq: String = seq[110 * unit..110 * unit + 100].iter().collect();
            read.nodes.push(gen_node(unit as u64, 0, true, subseq));
            let lab: String = seq[110 * unit + 100..110 * (unit + 1)].iter().collect();
            read.edges
                .push(gen_edge(unit as u64, (unit as u64 + 1) % 7, 10, lab));
        }
        let subseq: String = seq[..100].iter().collect();
        read.nodes.push(gen_node(0, 0, true, subseq));
        (vec![read], vec![seq.into_iter().collect()])
    }
    fn gen_complex() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(11212);
        let units: Vec<String> = (0..18)
            .map(|_| {
                (0..100)
                    .filter_map(|_| BASES.choose(&mut rng))
                    .copied()
                    .collect()
            })
            .collect();
        let mut reads = vec![];
        let (mut read0, mut read1) = (EncodedRead::default(), EncodedRead::default());
        for unit in 8..10 {
            let subseq = units[unit].clone();
            let unit = unit as u64;
            read0.nodes.push(gen_node(unit, 0, true, subseq.clone()));
            read1.nodes.push(gen_node(unit, 1, true, subseq.clone()));
            read0.edges.push(gen_edge(unit, unit + 1, 0, String::new()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, String::new()));
        }
        let unit = 10u64;
        let subseq = units[unit as usize].clone();
        read0.nodes.push(gen_node(unit, 0, true, subseq.clone()));
        read1.nodes.push(gen_node(unit, 0, true, subseq.clone()));
        reads.push(read0);
        reads.push(read1);
        let mut read = EncodedRead::default();
        for &unit in &[10u64, 11u64] {
            let subseq = units[unit as usize].clone();
            read.nodes.push(gen_node(unit, 0, true, subseq));
            read.edges.push(gen_edge(unit, unit + 1, 0, String::new()));
        }
        read.nodes.push(gen_node(12, 0, true, units[12].clone()));
        reads.push(read);
        let mut read = EncodedRead::default();
        read.nodes.push(gen_node(10, 0, true, units[10].clone()));
        read.edges.push(gen_edge(unit as u64, 13, 0, String::new()));
        for &unit in &[13, 14, 15, 16] {
            read.nodes
                .push(gen_node(unit, 0, true, units[unit as usize].clone()));
            read.edges.push(gen_edge(unit, unit + 1, 0, String::new()));
        }
        read.nodes.push(gen_node(17, 0, true, units[17].clone()));
        read.edges.push(gen_edge(17, 11, 0, String::new()));
        read.nodes.push(gen_node(11, 0, true, units[11].clone()));
        reads.push(read);
        let answer = vec![
            units[8].clone() + &units[9],
            units[8].clone() + &units[9],
            units[10].clone(),
            units[11].clone() + &units[12],
            units[13..=17].iter().fold(String::new(), |x, y| x + y),
        ];
        (reads, answer)
    }
    fn gen_mock_complex() -> DataPack {
        let (mut circle_reads, mut circle) = gen_circle();
        let (complex_reads, complex) = gen_complex();
        circle_reads.extend(complex_reads);
        circle.extend(complex);
        (circle_reads, circle)
    }
    fn gen_remove_test_1() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..3)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let edge = String::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..3 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1) = if unit == 1 { (0, 1) } else { (0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq));
        }
        for unit in 0..2 {
            read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
        }
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (0..3).rev() {
            let seq: Vec<_> = units[unit].chars().collect();
            let seq: String = revcmp(&seq).iter().collect();
            let unit = unit as u64;
            let (cl2, cl3) = if unit == 1 { (0, 1) } else { (0, 0) };
            read2.nodes.push(gen_node(unit, cl2, false, seq.clone()));
            read3.nodes.push(gen_node(unit, cl3, false, seq));
        }
        for unit in (0..2).rev() {
            read2.edges.push(gen_edge(unit + 1, unit, 0, edge.clone()));
            read3.edges.push(gen_edge(unit + 1, unit, 0, edge.clone()));
        }
        let answer = units.iter().fold(String::new(), |x, y| x + y);
        (vec![read0, read1, read2, read3], vec![answer])
    }
    fn gen_by_units(units: &[String], read: &[(u64, u64, bool)]) -> EncodedRead {
        let edge = String::new();
        let mut eread = EncodedRead::default();
        for &(u, c, b) in read.iter() {
            let seq = units[u as usize].clone();
            if b {
                eread.nodes.push(gen_node(u, c, b, seq));
            } else {
                let seq: Vec<char> = seq.chars().collect();
                let seq: String = seq.into_iter().collect();
                eread.nodes.push(gen_node(u, c, b, seq));
            }
        }
        for w in read.windows(2) {
            eread.edges.push(gen_edge(w[0].0, w[1].0, 0, edge.clone()));
        }
        eread
    }
    fn gen_remove_test_2() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..5)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let reads = vec![
            vec![(0, 0, true), (1, 0, true), (2, 0, true)],
            vec![(2, 0, false), (1, 0, false), (0, 0, false)],
            vec![(0, 0, true), (1, 1, true), (2, 0, true)],
            vec![(2, 0, false), (1, 1, false), (0, 0, false)],
            vec![(1, 1, true), (4, 0, true), (3, 0, true)],
            vec![(3, 0, false), (4, 0, false), (1, 1, false)],
        ];
        let reads: Vec<_> = reads
            .into_iter()
            .map(|r| gen_by_units(&units, &r))
            .collect();
        let answers = vec![
            units[0..2].iter().flat_map(|e| e.chars()).collect(),
            units[2].clone(),
            units[4].clone() + &units[3],
        ];
        (reads, answers)
    }
    fn gen_remove_test_3() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..3)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let edge = String::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1, mut read2) = (read.clone(), read.clone(), read.clone());
        for unit in 0..3 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1, cl2) = if unit == 1 { (0, 1, 2) } else { (0, 0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq.clone()));
            read2.nodes.push(gen_node(unit, cl2, true, seq));
            if unit < 2 {
                read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read2.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            }
        }
        let (mut read3, mut read4, mut read5) = (read.clone(), read.clone(), read.clone());
        for unit in (0..3).rev() {
            let seq: Vec<_> = units[unit].chars().collect();
            let seq: String = revcmp(&seq).iter().collect();
            let unit = unit as u64;
            let (cl3, cl4, cl5) = if unit == 1 { (0, 1, 2) } else { (0, 0, 0) };
            read3.nodes.push(gen_node(unit, cl3, false, seq.clone()));
            read4.nodes.push(gen_node(unit, cl4, false, seq.clone()));
            read5.nodes.push(gen_node(unit, cl5, false, seq.clone()));
            if 0 < unit {
                read3.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read4.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read5.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
            }
        }
        let reads = vec![read0, read1, read2, read3, read4, read5];
        let answers = units.into_iter().fold(String::new(), |x, y| x + &y);
        (reads, vec![answers])
    }
    fn gen_remove_test_4() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..5)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let reads = vec![
            vec![(0, 0, true), (1, 0, true), (2, 0, true)],
            vec![(0, 0, true), (1, 1, true), (2, 0, true)],
            vec![(0, 0, true), (3, 0, true), (4, 0, true)],
            vec![(4, 0, false), (3, 0, false), (0, 0, false)],
            vec![(2, 0, false), (1, 0, false), (0, 0, false)],
        ];
        let reads: Vec<_> = reads.iter().map(|r| gen_by_units(&units, r)).collect();
        let answers = vec![
            units[0].clone(),
            units[1].clone() + &units[2],
            units[3].clone() + &units[4],
        ];
        (reads, answers)
    }
    fn gen_remove_test_5() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..5)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let edge = String::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1) = (read.clone(), read.clone());
        for unit in 0..5 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1) = if unit == 1 { (0, 1) } else { (0, 0) };
            read0.nodes.push(gen_node(unit, cl0, true, seq.clone()));
            read1.nodes.push(gen_node(unit, cl1, true, seq));
            if unit < 4 {
                read0.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
                read1.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            }
        }
        let (mut read2, mut read3) = (read.clone(), read.clone());
        for unit in (0..5).rev() {
            let seq: Vec<_> = units[unit].chars().collect();
            let seq: String = revcmp(&seq).iter().collect();
            let unit = unit as u64;
            let (cl2, cl3) = if unit == 1 { (0, 1) } else { (0, 0) };
            read2.nodes.push(gen_node(unit, cl2, false, seq.clone()));
            read3.nodes.push(gen_node(unit, cl3, false, seq));
            if 0 < unit {
                read2.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
                read3.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
            }
        }
        let answer = units.iter().fold(String::new(), |x, y| x + y);
        (vec![read0, read1, read2, read3], vec![answer])
    }

    use super::*;
    #[test]
    fn create() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let _graph = DitchGraph::new(&reads, None, &c);
    }
    #[test]
    fn generate() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, &c);
        let _ = graph.spell(&c, 0);
    }
    #[test]
    fn validate() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        let (segment, _edge, _group, _) = graph.spell(&c, 0);
        eprintln!("{:?}", segment);
        eprintln!("{:?}", answer);
        assert_eq!(segment.len(), 1);
        let seq = segment[0].sequence.as_ref().unwrap();
        let revcmp = revcmp_str(segment[0].sequence.as_ref().unwrap());
        let is_the_same = (seq == answer[0].as_str()) || (answer[0].as_str() == &revcmp);
        assert!(is_the_same);
    }
    #[test]
    fn validate_large() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_large();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 4);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_large_2() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_large_2();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 4);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_complex() {
        let c = AssembleConfig::default();
        let (reads, _answer) = gen_mock_complex();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 6);
        assert_eq!(segments.len(), 6);
    }
    #[test]
    fn validate_remove_1() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_1();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_remove_2() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_2();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 2);
        assert_eq!(segments.len(), 3);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_remove_3() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_3();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_remove_4() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_4();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 2);
        assert_eq!(segments.len(), 3);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
    }
    #[test]
    fn validate_remove_5() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_remove_test_5();
        let reads: Vec<_> = reads.iter().collect();
        let mut graph = DitchGraph::new(&reads, None, &c);
        eprintln!("{:?}", graph);
        graph.collapse_bubble(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group, _) = graph.spell(&c, 0);
        let mut records = vec![];
        let nodes = segments
            .clone()
            .into_iter()
            .map(gfa::Content::Seg)
            .map(|n| gfa::Record::from_contents(n, vec![]));
        records.extend(nodes);
        {
            let edges = edges
                .clone()
                .into_iter()
                .map(|(e, _)| gfa::Content::Edge(e))
                .map(|n| gfa::Record::from_contents(n, vec![]));
            records.extend(edges);
        }
        let group = gfa::Record::from_contents(gfa::Content::Group(group.clone()), vec![]);
        records.push(group);
        eprintln!("=============Assembly================");
        eprintln!("{}", gfa::GFA::from_records(records));
        // Assertion.
        assert_eq!(edges.len(), 0);
        assert_eq!(segments.len(), 1);
        for ans in answer {
            let forward = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &ans);
            let rev_ans: Vec<_> = ans.chars().collect();
            let rev_ans: String = revcmp(&rev_ans).into_iter().collect();
            let reverse = segments
                .iter()
                .any(|seg| seg.sequence.as_ref().unwrap() == &rev_ans);
            assert!(forward || reverse, "Couldn't find {}", ans);
        }
        assert!(graph.sanity_check());
    }
}
