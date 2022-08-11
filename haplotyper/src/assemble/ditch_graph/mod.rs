use super::AssembleConfig;
use definitions::DNASeq;
use definitions::{EncodedRead, Unit};
mod copy_num_by_flow;
mod copy_num_by_mst;
pub mod sequence_generation;
pub use sequence_generation::*;
use std::collections::HashMap;
use std::collections::HashSet;
mod iterators;
mod position;
mod repeat_resolve_by_focus;
mod squish_graph;
mod update_copy_numbers;
use position::Position;
pub mod dg_test;

type Node = (u64, u64);
type DitEdge = ((NodeIndex, Position), (NodeIndex, Position));
type EdgeBetweenSimplePath = (usize, bool, usize, bool, f64);
type GraphNode = (NodeIndex, Position);
type GraphBoundary = Vec<GraphNode>;

/// Type to index nodes in the graph.
#[derive(Debug, Copy, Clone, Eq, PartialEq, Hash, PartialOrd, Ord)]
pub struct NodeIndex(usize);

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::convert::From<NodeIndex> for usize {
    fn from(NodeIndex(i): NodeIndex) -> Self {
        i
    }
}

/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
/// It is bi-directed graph, in other words, if you see (from,from_pos, to, to_pos),
/// then, there is a edge (to, to_pos, from, from_pos).
/// My implementation notes:
/// - Pay attention to the `is_deleted` flag on each node!
/// - All edges *should not* connected to/from deleted nodes.
#[derive(Clone)]
pub struct DitchGraph<'a> {
    /// Node and the next node representing the same (unit,cluster).
    nodes: Vec<DitchNode<'a>>,
    nodes_index: HashMap<Node, NodeIndex>,
}

impl<'a> std::fmt::Display for DitchGraph<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        use histgram_viz::Histgram;
        let nodes = self.active_nodes();
        let edges = self.edges().count();
        writeln!(f, "Node:{}, Edges:{}", nodes, edges)?;
        let occs: Vec<_> = self.nodes().map(|(_, n)| n.occ).collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Node Occs:{}", hist.format(20, 20))?;
        let occs: Vec<_> = self
            .nodes()
            .flat_map(|(_, n)| n.edges.iter().map(|e| e.occ))
            .collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Edge Occs:{}", hist.format(20, 20))?;
        let degrees = {
            let mut degs: HashMap<usize, usize> = HashMap::new();
            for (_, node) in self.nodes() {
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
        let edge = self.edges().count();
        writeln!(f, "Nodes:{}\tEdges:{}\n", self.active_nodes(), edge)?;
        let lines: Vec<_> = self.nodes().map(|(_, node)| format!("{}", node)).collect();
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
    pub node: Node,
    pub occ: usize,
    pub seq: Vec<u8>,
    pub edges: Vec<DitchEdge>,
    // "Tip" of reads. In other words, as we tiling a read by units,
    // there is un-encoded regions at the both end of a read,
    // and we allocate memories for them.
    pub tips: Vec<DitchTip<'a>>,
    // Estimated copy number. If not able, None.
    pub copy_number: Option<usize>,
    // If true, this node is deleted.
    pub is_deleted: bool,
    // If Some(idx), nodes[idx] has the same `node` member.
    pub next_index: Option<NodeIndex>,
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
            is_deleted: false,
            next_index: None,
        }
    }
    // pub fn delete(&mut self) {
    //     self.is_deleted = true;
    // }
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
    fn seq_as_string(&self) -> String {
        String::from_utf8(self.seq.clone()).unwrap()
    }
    fn decrement_copy_number(&mut self) -> usize {
        match self.copy_number {
            Some(cp) if 0 < cp => {
                let occ = self.occ / cp;
                self.occ -= occ;
                self.copy_number = Some(cp - 1);
                occ
            }
            _ => 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DitchEdge {
    pub from: NodeIndex,
    pub to: NodeIndex,
    from_node: Node,
    to_node: Node,
    from_position: Position,
    to_position: Position,
    seq: EdgeLabel,
    pub occ: usize,
    // Estimated copy number. None if not available.
    copy_number: Option<usize>,
}

impl std::fmt::Display for DitchEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({:?}:{})->", self.from_node, self.from_position)?;
        write!(f, "({:?}:{})", self.to_node, self.to_position)?;
        let cp = self.copy_number.unwrap_or(0);
        write!(f, "{}:{}:{}", self.occ, self.seq, cp)
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
    pub fn key(&self) -> DitEdge {
        ((self.from, self.from_position), (self.to, self.to_position))
    }
    /// Normalized key (f,t) where f <= t.
    pub fn norm_key(&self) -> DitEdge {
        let (f, t) = self.key();
        (f.min(t), f.max(t))
    }
    pub fn occ(&self) -> usize {
        self.occ
    }
    fn new(
        (from, from_node, from_position): (NodeIndex, Node, Position),
        (to, to_node, to_position): (NodeIndex, Node, Position),
        seq: EdgeLabel,
        occ: i64,
    ) -> Self {
        Self {
            from,
            to,
            from_node,
            from_position,
            to_node,
            to_position,
            occ: occ.max(0) as usize,
            seq,
            copy_number: None,
        }
    }
    pub fn label(&self) -> Option<&[u8]> {
        match &self.seq {
            EdgeLabel::Ovlp(_) => None,
            EdgeLabel::Seq(seq) => Some(seq.as_slice()),
        }
    }
    // Reverse this edge.
    fn reverse(&self) -> Self {
        Self {
            from_node: self.to_node,
            from: self.to,
            from_position: self.to_position,
            to: self.from,
            to_node: self.from_node,
            to_position: self.from_position,
            seq: self.seq.reverse(),
            occ: self.occ,
            copy_number: self.copy_number,
        }
    }
    // reduce the copy number by one,
    // and decrease the .occ parameter appropriately.
    fn reduce_one_copy_number(&mut self) -> usize {
        match self.copy_number.as_mut() {
            Some(cp) if *cp == 0 => 0,
            Some(cp) if *cp == 1 => {
                let occ = self.occ;
                self.occ = 0;
                *cp = 0;
                occ
            }
            Some(cp) => {
                let occ_reduce = self.occ / *cp;
                self.occ -= occ_reduce;
                *cp -= 1;
                occ_reduce
            }
            None => 0,
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
    fn new(seq: &'a DNASeq, position: Position, in_direction: bool) -> Self {
        Self {
            seq: seq.as_slice(),
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

fn take_representative<R: std::borrow::Borrow<EncodedRead>>(
    reads: &[R],
    node_to_idx: &HashMap<Node, NodeIndex>,
) -> Vec<DitchEdge> {
    use std::collections::BTreeMap;
    let mut edges: BTreeMap<_, (i64, i64, Vec<u8>)> = BTreeMap::new();
    use Position::*;
    for read in reads.iter().map(|r| r.borrow()) {
        for (i, edge) in read.edges.iter().enumerate() {
            let from_node = read.nodes.get(i).map(|n| (n.unit, n.cluster)).unwrap();
            let to_node = read.nodes.get(i + 1).map(|n| (n.unit, n.cluster)).unwrap();
            let from = node_to_idx[&from_node];
            let to = node_to_idx[&to_node];
            let from_pos = if read.nodes[i].is_forward { Tail } else { Head };
            let to_pos = match read.nodes[i + 1].is_forward {
                true => Head,
                false => Tail,
            };
            // No bigger element...
            let (key, is_forward) = match (from, from_pos) <= (to, to_pos) {
                true => (((from, from_node, from_pos), (to, to_node, to_pos)), true),
                false => (((to, to_node, to_pos), (from, from_node, from_pos)), false),
            };
            let lab = edges.entry(key).or_default();
            lab.0 += edge.offset;
            lab.1 += 1;
            if lab.2.len() < edge.label().len() {
                lab.2 = match is_forward {
                    true => edge.label().to_vec(),
                    false => bio_utils::revcmp(edge.label()),
                };
            }
        }
    }
    edges
        .into_iter()
        .map(|(key, (offset, edge_num, seq))| {
            let seq = match seq.is_empty() {
                true => EdgeLabel::Ovlp(offset / edge_num),
                false => EdgeLabel::Seq(seq),
            };
            (key, seq, edge_num)
        })
        .map(|((from_elm, to_elm), seq, occ)| DitchEdge::new(from_elm, to_elm, seq, occ))
        .collect()
}

impl<'b, 'a: 'b> DitchGraph<'a> {
    pub fn new<R: std::borrow::Borrow<EncodedRead>>(
        reads: &'a [R],
        units: &[Unit],
        _read_type: definitions::ReadType,
        c: &AssembleConfig,
    ) -> Self {
        let nodes_seq: HashMap<_, _> = units.iter().map(|c| (c.id, c.seq())).collect();
        use std::collections::BTreeMap;
        let mut nodes_counts: BTreeMap<_, usize> = BTreeMap::new();
        for node in reads.iter().flat_map(|r| r.borrow().nodes.iter()) {
            *nodes_counts.entry((node.unit, node.cluster)).or_default() += 1;
        }
        let (nodes, nodes_index): (Vec<_>, HashMap<Node, _>) = nodes_counts
            .into_iter()
            .enumerate()
            .map(|(idx, (node, occ))| {
                let seq = nodes_seq[&node.0].to_vec();
                let mut d_node = DitchNode::new(node, seq);
                d_node.occ = occ;
                (d_node, (node, NodeIndex(idx)))
            })
            .unzip();
        for (&node, &idx) in nodes_index.iter() {
            assert_eq!(nodes[idx.0].node, node);
        }
        let edge_seq: Vec<_> = take_representative(reads, &nodes_index);
        let mut graph = Self { nodes, nodes_index };
        assert!(graph.sanity_check(), "{}", line!());
        for edge in edge_seq.into_iter() {
            graph.add_edge(edge);
        }
        for read in reads.iter().map(|r| r.borrow()) {
            graph.append_tip(read, c);
        }
        assert!(graph.sanity_check(), "{}", line!());
        graph
    }
    fn add_edge(&'b mut self, edge: DitchEdge) {
        let (f, t) = edge.key();
        if f != t {
            self.node_mut(edge.to).unwrap().edges.push(edge.reverse());
        }
        self.node_mut(edge.from).unwrap().edges.push(edge);
    }
    // Currently, just update the occ of the edge. Make sure that it truly has the edge.
    fn merge_edge(&'b mut self, edge: &DitchEdge) {
        let forward = self
            .node_mut(edge.from)
            .unwrap()
            .edges
            .iter_mut()
            .find(|e| {
                e.from_position == edge.from_position
                    && e.to == edge.to
                    && e.to_position == edge.to_position
            })
            .unwrap();
        forward.occ += edge.occ;
        if let Some(cp) = forward.copy_number.as_mut() {
            *cp += edge.copy_number.unwrap_or(0);
        }
        let reverse = self
            .node_mut(edge.to)
            .unwrap()
            .edges
            .iter_mut()
            .find(|e| {
                e.from_position == edge.to_position
                    && e.to == edge.from
                    && e.to_position == edge.from_position
            })
            .unwrap();
        reverse.occ += edge.occ;
        if let Some(cp) = reverse.copy_number.as_mut() {
            *cp += edge.copy_number.unwrap_or(0);
        }
    }
    fn has_edge(&self, ((f, fpos), (t, tpos)): DitEdge) -> bool {
        let f_has = self
            .node(f)
            .map(|node| {
                node.edges
                    .iter()
                    .any(|e| e.from_position == fpos && e.to == t && e.to_position == tpos)
            })
            .unwrap_or(false);
        let t_has = self
            .node(t)
            .map(|node| {
                node.edges
                    .iter()
                    .any(|e| e.from_position == tpos && e.to == f && e.to_position == fpos)
            })
            .unwrap_or(false);
        assert_eq!(f_has, t_has);
        f_has
    }
    // Add weight of the graph.
    fn append_tip(&'b mut self, read: &'a EncodedRead, _c: &AssembleConfig) -> Option<()> {
        use Position::*;
        if let Some(first) = read.nodes.first() {
            let direction = true;
            let position = if first.is_forward { Head } else { Tail };
            let index = self.nodes_index[&(first.unit, first.cluster)];
            let tip = DitchTip::new(&read.leading_gap, position, direction);
            self.node_mut(index).unwrap().tips.push(tip);
        }
        if let Some(last) = read.nodes.last() {
            let direction = true;
            let position = if last.is_forward { Head } else { Tail };
            let index = self.nodes_index[&(last.unit, last.cluster)];
            let tip = DitchTip::new(&read.trailing_gap, position, direction);
            self.node_mut(index).unwrap().tips.push(tip);
        }
        Some(())
    }
    pub fn sanity_check(&self) -> bool {
        for (_, node) in self.nodes() {
            assert!(!node.is_deleted);
            let mut edges = HashSet::new();
            for e in node.edges.iter() {
                if edges.contains(&e.key()) {
                    error!("DUPLCATE EDGE\t{}", node);
                    return false;
                }
                edges.insert(e.key());
            }
        }
        for (i, node) in self.nodes() {
            for edge in node.edges.iter() {
                if edge.from != i {
                    error!("EDGE INCONSIS\t{node}\t{edge}\t{i}");
                    return false;
                }
                if edge.from_node != node.node {
                    error!("EDGE NOMATCH\t{node}\t{edge}");
                    return false;
                }
                let rev = edge.reverse();
                assert!(!self.node(edge.to).unwrap().is_deleted);
                let mut match_edge = self
                    .node(edge.to)
                    .unwrap()
                    .edges
                    .iter()
                    .filter(|e| e == &&rev);
                let rev_edge = match match_edge.next() {
                    Some(e) => e,
                    None => {
                        error!("NO REV EDGE\t{}\t{}", node, edge);
                        return false;
                    }
                };
                if rev_edge.copy_number != edge.copy_number {
                    error!("MULTIPNOTAGREE\n{rev_edge}\n{edge}");
                    return false;
                }
                if match_edge.next().is_some() {
                    error!("REV EDGE MULTIP\t{}", node);
                    error!("{edge}");
                    return false;
                }
            }
        }
        true
    }
    pub fn active_nodes(&self) -> usize {
        self.nodes().count()
    }
    pub fn is_deleted(&self, NodeIndex(i): NodeIndex) -> bool {
        self.nodes.get(i).map(|n| n.is_deleted).unwrap()
    }
    pub fn delete(&'b mut self, node: NodeIndex) {
        let removed: Vec<_> = self
            .node(node)
            .unwrap()
            .edges
            .iter()
            .map(|e| e.key())
            .collect();
        for (f, t) in removed {
            self.remove_edge_between(f, t);
        }
        self.node_mut(node).unwrap().is_deleted = true;
    }
    /// Get node at the `i`-th index. If it is deleted, panic.
    pub fn node(&'b self, NodeIndex(i): NodeIndex) -> Option<&'b DitchNode<'a>> {
        match self.nodes.get(i) {
            Some(n) if n.is_deleted => {
                error!("{:?}", n.node);
                panic!("Accessing a deleted node.")
            }
            x => x,
        }
    }
    pub fn edges_from(&self, node: NodeIndex, position: Position) -> Vec<&DitchEdge> {
        self.node(node)
            .unwrap()
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .collect()
    }
    pub fn node_mut(&'b mut self, NodeIndex(i): NodeIndex) -> Option<&'b mut DitchNode<'a>> {
        match self.nodes.get_mut(i) {
            Some(n) if n.is_deleted => panic!("Accessing a deleted node."),
            x => x,
        }
    }
    /// Duplicate given node, return newly allocated index.
    /// The `self.node(index).unwrap().next_index` might not be the same as the
    /// returned value.
    /// This function only duplicate node, not edges with it.
    pub fn duplicate(&'b mut self, index: NodeIndex) -> NodeIndex {
        let mut node = self.node(index).unwrap().clone();
        node.is_deleted = false;
        node.next_index = None;
        node.edges.clear();
        self.nodes.push(node);
        let len = self.nodes.len() - 1;
        let mut par_index = index;
        loop {
            // Here, self.node(par_index) is not appropriate, as
            // sometimes the child is removed...
            par_index = match self.nodes[par_index.0].next_index {
                Some(next) => next,
                None => break,
            };
        }
        let par_node = self.node_mut(par_index).unwrap();
        par_node.next_index = Some(NodeIndex(len));
        NodeIndex(len)
    }
    /// Check if there is a node labeled with `node`.
    pub fn contains_node(&self, node: Node) -> bool {
        self.nodes_index.contains_key(&node)
    }
    // POINTER: ASSEMBLEIMPL
    pub fn clean_up_graph_for_assemble(
        &'b mut self,
        cov: f64,
        reads: &[&EncodedRead],
        c: &super::AssembleConfig,
        _read_type: definitions::ReadType,
    ) {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256PlusPlus;
        let seed = self.nodes.len() as u64 * 7329;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        debug!("CC\tBFASN\t{}", self.cc());
        if log_enabled!(log::Level::Trace) {
            dump(self, 0, c);
        }
        self.assign_copy_number(cov, &mut rng);
        self.remove_tips(0.8, 4);
        assert!(self.sanity_check());
        debug!("CC\tRMZERO\t{}", self.cc());
        self.remove_tips(0.8, 4);
        // From good Likelihood ratio focus, to weaker ones.
        let min_llr = c.span_likelihood_ratio;
        let llr_stream = ((min_llr + 0.0).floor() as usize..(10.0 + min_llr).floor() as usize)
            .rev()
            .map(|i| i as f64 + 0.00001)
            .take_while(|&x| min_llr < x)
            .enumerate();
        for (i, llr) in llr_stream {
            self.assign_copy_number(cov, &mut rng);
            self.remove_zero_copy_elements(0.8);
            debug!("REPEATRESOLVE\t{}", i);
            self.remove_zero_copy_path(0.1);
            self.resolve_repeats(reads, c, llr as f64, false);
            debug!("CC\tSOLVEREP\t{}\t{i}", self.cc());
            self.zip_up_overclustering(2);
            if i == 5 {
                self.assign_copy_number(cov, &mut rng);
                self.remove_zero_copy_elements(0.9);
                self.remove_zero_copy_path(0.3);
                self.remove_lightweight_edges(0, true);
                self.remove_tips(0.8, 4);
                self.squish_small_net(3);
                assert!(self.sanity_check(), "{}", line!());
            }
            if log_enabled!(log::Level::Trace) {
                dump(self, i + 1, c);
            }
        }
        self.assign_copy_number(cov, &mut rng);
        self.remove_zero_copy_elements(0.9);
        self.remove_zero_copy_path(0.3);
        self.remove_lightweight_edges(0, true);
        self.remove_tips(0.8, 4);
        self.squish_small_net(3);
        self.assign_copy_number(cov, &mut rng);
        self.zip_up_overclustering_dev();
        self.resolve_repeats(reads, c, min_llr, true);
        self.remove_zero_copy_elements(100f64);
        // self.z_edge_selection();
        // self.remove_zero_copy_path(0.2);
    }
    /// Retun the number of the connected components
    pub fn cc(&self) -> usize {
        // Connected component.
        let mut fu = crate::find_union::FindUnion::new(self.nodes.len());
        for (_, node) in self.nodes() {
            for edge in node.edges.iter() {
                fu.unite(edge.from.0, edge.to.0);
            }
        }
        self.nodes()
            .filter(|(idx, _)| fu.find(idx.0).unwrap() == idx.0)
            .count()
    }
    pub fn connected_components(&self) -> Vec<Vec<Node>> {
        let mut fu = crate::find_union::FindUnion::new(self.nodes.len());
        for (_, node) in self.nodes() {
            for edge in node.edges.iter() {
                fu.unite(edge.from.0, edge.to.0);
            }
        }
        let mut clusters: HashMap<_, Vec<_>> = HashMap::new();
        for (idx, node) in self.nodes() {
            let elm = node.node;
            clusters
                .entry(fu.find(idx.0).unwrap())
                .or_default()
                .push(elm);
        }
        clusters.into_values().collect()
    }
}

fn dump(graph: &DitchGraph, i: usize, c: &AssembleConfig) {
    let (segments, edge, _group, summaries, _) = graph.spell(c);
    let mut groups: HashMap<_, Vec<_>> = HashMap::new();
    let nodes: Vec<_> = segments
        .into_iter()
        .map(|mut node| {
            let tags = match summaries.iter().find(|x| x.id == node.sid) {
                Some(contigsummary) => {
                    let ids: Vec<_> = contigsummary
                        .summary
                        .iter()
                        .map(|elm| format!("{}-{}", elm.unit, elm.cluster))
                        .collect();
                    let total: usize = contigsummary.summary.iter().map(|n| n.occ).sum();
                    let coverage =
                        gfa::SamTag::new(format!("cv:i:{}", total / contigsummary.summary.len()));
                    let (cp, cpnum) = contigsummary
                        .summary
                        .iter()
                        .filter_map(|elm| elm.copy_number)
                        .fold((0, 0), |(cp, num), x| (cp + x, num + 1));
                    let copynum = (cp as f64 / cpnum.max(1) as f64).round() as usize;
                    log::debug!(
                        "ASMDUMP\t{i}\t{}\t{copynum}\t{}",
                        contigsummary.id,
                        ids.join("\t")
                    );
                    groups.entry(copynum).or_default().push(node.sid.clone());
                    let cp = gfa::SamTag::new(format!("cp:i:{copynum}"));
                    vec![coverage, cp]
                }
                None => Vec::new(),
            };
            node.sequence = Some("A".to_string());
            gfa::Record::from_contents(gfa::Content::Seg(node), tags.into())
        })
        .collect();
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags.into()));
    let groups = groups.into_iter().map(|(cp, ids)| {
        let group = gfa::UnorderedGroup {
            uid: Some(format!("cp:i:{}", cp)),
            ids,
        };
        let group = gfa::Content::Group(gfa::Group::Set(group));
        gfa::Record::from_contents(group, vec![].into())
    });
    // let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![].into());
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![].into());
    let records = std::iter::once(header)
        .chain(groups)
        .chain(nodes)
        .chain(edges)
        .collect();
    let gfa = gfa::GFA::from_records(records);
    if let Ok(mut wtr) = std::fs::File::create(format!("{}.gfa", i)).map(std::io::BufWriter::new) {
        use std::io::Write;
        if let Err(why) = writeln!(wtr, "{}", gfa) {
            trace!("{:?}", why);
        }
    }
}

impl<'b, 'a: 'b> DitchGraph<'a> {
    // Return primary node.
    // fn get_primary_node(&'b self, node: Node) -> &'b DitchNode<'a> {
    //     self.nodes_index
    //         .get(&node)
    //         .and_then(|&i| self.node(i))
    //         .unwrap_or_else(|| panic!("{:?} does not exists.", node))
    // }
    // // Return iterator yeilding the edge from a specified node and position.
    // fn get_primary_edges(
    //     &self,
    //     node: Node,
    //     pos: Position,
    // ) -> impl std::iter::Iterator<Item = &DitchEdge> {
    //     self.get_primary_node(node)
    //         .edges
    //         .iter()
    //         .filter(move |e| e.from_position == pos)
    // }
    fn count_edges(&self, node_idx: NodeIndex, position: Position) -> usize {
        self.node(node_idx)
            .unwrap()
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count()
    }
    // Decrease the copy number of the given edge by one,
    // tune the occurence of these edges appropriately.
    fn decrement_edge_copy_number(
        &'b mut self,
        (from_idx, from_pos): (NodeIndex, Position),
        (to_idx, to_pos): (NodeIndex, Position),
    ) -> usize {
        let (mut num, mut total) = (0, 0);
        for edge in self.node_mut(from_idx).unwrap().edges.iter_mut() {
            if edge.from_position == from_pos && edge.to == to_idx && edge.to_position == to_pos {
                total += edge.reduce_one_copy_number();
                num += 1;
            }
        }
        for edge in self.node_mut(to_idx).unwrap().edges.iter_mut() {
            if edge.from_position == to_pos && edge.to == from_idx && edge.to_position == from_pos {
                total += edge.reduce_one_copy_number();
                num += 1;
            }
        }
        total / num.max(1)
    }
    // Return tuples of node and their position from which
    // we can start simple-path-reduction.
    // In other words, it returns the list of the position which
    // does not has any edge in the opposite position.
    // For example, if there is no edge from Position::Tail at the first node,
    // then, (0, Position::Head) would be included.
    // Also, it returns the list of the node whose parent has two or more children.
    fn enumerate_candidates(&self) -> GraphBoundary {
        let mut selected = vec![false; self.nodes.len()];
        let mut primary_candidates = vec![];
        for (idx, node) in self.nodes() {
            for position in [Position::Head, Position::Tail] {
                let num_of_edge = node
                    .edges
                    .iter()
                    .filter(|e| e.from_position == position)
                    .count();
                if num_of_edge != 1 {
                    selected[idx.0] = true;
                    primary_candidates.push((idx, position));
                }
            }
        }
        // Secondary candidates.
        for (index, node) in self.nodes().filter(|&(idx, _)| !selected[idx.0]) {
            for position in [Position::Head, Position::Tail] {
                let grand_child = node
                    .edges
                    .iter()
                    .find(|e| e.from_position == position)
                    .map(|e| {
                        let count = self.count_edges(e.to, e.to_position);
                        if count == 0 {
                            debug!("{:?},{}", node, position);
                            debug!("Edge:{:?}", e);
                            for edge in &self.node(e.to).unwrap().edges {
                                debug!("Child:{:?}", edge)
                            }
                            panic!()
                        }
                        count
                    })
                    .unwrap();
                if 1 < grand_child {
                    primary_candidates.push((index, position))
                }
            }
        }
        primary_candidates
    }

    // fn generate_coverage_calib(&self, naive_cov: f64, lens: &[usize]) -> (CoverageCalibrator, f64) {
    //     let calibrator = CoverageCalibrator::new(lens);
    //     let unit_len_sum: usize = self.nodes().map(|n| n.seq().len()).sum();
    //     let cov = calibrator.calib_f64(naive_cov, unit_len_sum / self.nodes.len());
    //     debug!("COPYNUM\tCOVERAGE\t{:.3}\t{:.3}", naive_cov, cov,);
    //     (calibrator, cov)
    // }

    /// Even though the edge/node is zero copy number, we do not remove it if the conditions below hold:
    /// 1. If it is an edge, and all edge from the same position is ZCP, and it is the heviest edge among them.
    /// 2. If it is a node, and it has un-removable edge.
    /// 3. If it is a node, and the bounding constraint does not hold.
    ///    In other words, it connected to the non-ZCP edge.
    /// 4. If it is an edge, and the occ is more than `thr * max_out_dgree`
    pub fn remove_zero_copy_elements(&mut self, thr: f64) {
        // Check the node violating for right-left condition.
        let unsound_nodes: HashSet<NodeIndex> = self
            .nodes()
            .filter_map(|(index, node)| {
                let degrees: [usize; 2] = node.edges.iter().fold([0, 0], |mut copy_nums, edge| {
                    if let Some(cp) = edge.copy_number {
                        copy_nums[(edge.from_position == Position::Head) as usize] += cp;
                    }
                    copy_nums
                });
                match degrees {
                    [0, _] | [_, 0] => None,
                    [x, y] if x != y => Some(index),
                    _ => None,
                }
            })
            .collect();
        debug!("UNSOUND\t{}\t{}", unsound_nodes.len(), self.nodes.len());
        let mut is_ok_to_remove = HashSet::new();
        let mut retain_edges = HashSet::new();
        for (index, node) in self.nodes() {
            // If the estimation of the copy number is is poor, do not remove edges.
            if unsound_nodes.contains(&index) {
                retain_edges.extend(node.edges.iter().map(|e| e.norm_key()));
                continue;
            }
            for position in [Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == position);
                let max = edges.clone().map(|e| e.occ).max().unwrap_or(0);
                for edge in edges {
                    if matches!(edge.copy_number, Some(0)) && edge.occ as f64 / (max as f64) < thr {
                        is_ok_to_remove.insert(edge.norm_key());
                    } else {
                        retain_edges.insert(edge.norm_key());
                    }
                }
            }
        }
        let (mut max_removed_weight, mut argmax, mut removed_num) = (0, None, 0);
        for &(f, t) in is_ok_to_remove.difference(&retain_edges) {
            if let Some(removed) = self.remove_edge_between(f, t) {
                assert_eq!(removed.copy_number, Some(0), "{}", removed);
                removed_num += 1;
                if max_removed_weight < removed.occ {
                    max_removed_weight = removed.occ;
                    argmax = Some(removed.clone());
                }
            }
        }
        if let Some(e) = argmax {
            debug!("REMOVED\t{removed_num}\t{max_removed_weight}");
            trace!("MAX_EDGE\t{e}");
        }
        self.modify_by(|node| {
            if matches!(node.copy_number, Some(0)) && node.edges.is_empty() {
                node.is_deleted = true;
            }
        });
    }
    /// Remove zero-copy path.
    /// 1. Find a zero-copy path.
    /// 2. Check whether the destination of the zero-copy path are reachable from other paths.
    /// 3. If so, remove that path.
    pub fn remove_zero_copy_path(&'b mut self, thr: f64) {
        // enumerate parent nodes of the single-copy-path.
        let mut parents_of_zc = vec![];
        for (index, node) in self.nodes() {
            // Parents should have non-zero-copy.
            if matches!(node.copy_number, Some(0) | None) {
                continue;
            }
            for pos in [Position::Head, Position::Tail] {
                let eds = node.edges.iter().filter(|e| e.from_position == pos);
                if 2 <= eds.clone().count() {
                    let has_zc = eds
                        .filter_map(|e| self.node(e.to))
                        .any(|n| matches!(n.copy_number, Some(0)));
                    if has_zc {
                        parents_of_zc.push((index, pos));
                    }
                }
            }
        }
        // Check reachability.
        for (parent_index, par_pos) in parents_of_zc {
            if !self.is_deleted(parent_index) {
                continue;
            }

            if self.count_edges(parent_index, par_pos) <= 1 {
                continue;
            }
            let edges = self.node(parent_index).unwrap().edges.iter();
            let has_zc = edges
                .filter(|e| e.from_position == par_pos)
                .filter_map(|e| self.node(e.to).and_then(|n| n.copy_number))
                .any(|cp| cp == 0);
            if has_zc {
                // Check the destination of the zero-copy path and non-zero-copy paths.
                let (zc_edges, nzc_edges): (Vec<_>, Vec<_>) = self
                    .node(parent_index)
                    .unwrap()
                    .edges
                    .iter()
                    .filter(|e| e.from_position == par_pos)
                    .partition(|e| matches!(self.node(e.to).and_then(|n| n.copy_number), Some(0)));
                let zc_dests: HashSet<_> = zc_edges
                    .iter()
                    .flat_map(|e| self.simple_path_and_dest(e.to, e.to_position).1)
                    .collect();
                let zc_max = zc_edges.iter().map(|e| self.node(e.to).unwrap().occ).max();
                let nzc_dests: HashSet<_> = nzc_edges
                    .iter()
                    .flat_map(|e| self.simple_path_and_dest(e.to, e.to_position).1)
                    .collect();
                let nzc_max = nzc_edges.iter().map(|e| self.node(e.to).unwrap().occ).max();
                let zc_nzc_ratio = match (zc_max, nzc_max) {
                    (Some(x), Some(y)) => x as f64 / y as f64,
                    (Some(_), _) => 1f64,
                    _ => panic!(),
                };
                if zc_dests.is_subset(&nzc_dests) && zc_nzc_ratio < thr {
                    debug!(
                        "RZCP\t{:?}\t{}\t{}\t{:.3}",
                        parent_index,
                        par_pos,
                        zc_dests.is_subset(&nzc_dests),
                        zc_nzc_ratio
                    );
                    let targets: Vec<_> = zc_edges.iter().map(|e| (e.to, e.to_position)).collect();
                    for (node, pos) in targets {
                        self.remove_edge_and_pruning((parent_index, par_pos), (node, pos));
                    }
                }
            }
        }
    }
    /// Remove transitive edge if the copy_number of that edge is zero.
    pub fn transitive_edge_reduction(&mut self) {
        let mut removed_edges: HashMap<_, Vec<_>> =
            self.nodes().map(|(index, _)| (index, vec![])).collect();
        for (from, node) in self.nodes() {
            for &pos in &[Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                let to_cands = edges.clone().count();
                if to_cands > 1 {
                    for e in edges.clone().filter(|e| !self.is_deleted(e.to)) {
                        // Check whether or not this edge is transitive edge.
                        let is_transitive = edges.clone().any(|edge| {
                            self.node(edge.to)
                                .unwrap()
                                .edges
                                .iter()
                                .filter(|gc| gc.from_position == !edge.to_position)
                                .any(|gc| gc.to == e.to && gc.to_position == e.to_position)
                        });
                        if is_transitive && matches!(e.copy_number, Some(0)) {
                            let (to, t_pos) = (e.to, e.to_position);
                            removed_edges
                                .entry(from)
                                .or_default()
                                .push((pos, to, t_pos));
                            removed_edges
                                .entry(to)
                                .or_default()
                                .push((t_pos, from, pos));
                        }
                    }
                }
            }
        }
        self.modify_by_with_index(|(index, node)| {
            let to_be_removed = &removed_edges[&index];
            node.edges.retain(|x| {
                let probe = (x.from_position, x.to, x.to_position);
                !to_be_removed.contains(&probe)
            })
        });
    }

    /// Zip up overclustered regions.
    /// The graph should have estimeted repeat numbers on nodes.
    pub fn zip_up_overclustering(&mut self, len: usize) {
        let mut to_remove = HashSet::new();
        let nodes = self.nodes().filter(|n| matches!(n.1.copy_number, Some(1)));
        for (_index, node) in nodes {
            for &pos in &[Position::Head, Position::Tail] {
                let edges = node
                    .edges
                    .iter()
                    .filter(|e| e.from_position == pos && !to_remove.contains(&e.to));
                if edges.clone().count() != 2 {
                    continue;
                }
                // This is a "fork" branch.
                let mut dests = edges.clone().map(|edge| self.destination(edge));
                let first_dest = dests.next().unwrap();
                let second_dest = dests.next().unwrap();
                assert_eq!(dests.next(), None);
                if first_dest == second_dest {
                    let from_edge = edges.clone().max_by_key(|e| e.occ).unwrap();
                    let path = self.simple_path_from(from_edge);
                    debug!("ZIPPINGUP\t{:?}\t{}\t{}", node.node, pos, path.len());
                    if path.len() <= len {
                        to_remove.extend(path);
                    }
                }
            }
        }
        self.remove_nodes(&to_remove);
    }
    fn has_self_loop(&self, index: NodeIndex) -> bool {
        let node = self.node(index).unwrap();
        let unit_cl = node.node;
        node.edges.iter().any(|e| e.to_node == unit_cl)
    }
    pub fn zip_up_overclustering_dev(&'b mut self) {
        let mut keys: Vec<_> = self.nodes().map(|n| n.0).collect();
        keys.sort_unstable();
        for node in keys {
            if self.is_deleted(node) || self.has_self_loop(node) {
                continue;
            }
            let (retain, sibs) = match self.zippable(node) {
                Some(res) => res,
                None => continue,
            };
            // Make all the edges into sibs to retain.
            let (edges, increase_occ, increase_copy_num) = {
                let (mut edges, mut occ, mut cp) = (vec![], 0, 0);
                for &node in sibs.iter() {
                    let removed = self.node(node).unwrap();
                    edges.extend(removed.edges.clone());
                    occ += removed.occ;
                    cp += removed.copy_number.unwrap_or(0);
                    self.delete(node);
                }
                let retain_node = self.node(retain).unwrap().node;
                for edge in edges.iter_mut() {
                    edge.from = retain;
                    edge.from_node = retain_node;
                }
                (edges, occ, cp)
            };
            {
                let retain_node = self.node_mut(retain).unwrap();
                retain_node.occ += increase_occ;
                if let Some(cp) = retain_node.copy_number.as_mut() {
                    *cp += increase_copy_num;
                }
            };
            for edge in edges {
                if self.has_edge(edge.key()) {
                    self.merge_edge(&edge);
                } else {
                    self.add_edge(edge);
                }
            }
        }
    }

    // Check if both side of this node is either
    // 1. Branching
    // 2. Connected to branching node.
    // Reflexitive siblings, parents.
    fn zippable(&self, node: NodeIndex) -> Option<(NodeIndex, Vec<NodeIndex>)> {
        let (tail_side_par, tail_side_sibs) = self.get_reflex_nodes(node, Position::Tail, 6);
        let (head_side_par, head_side_sibs) = self.get_reflex_nodes(node, Position::Head, 6);
        if tail_side_sibs.len().max(head_side_sibs.len()) <= 1 {
            // This is not a "net-like structure"
            return None;
        }
        if head_side_par.is_empty() || tail_side_par.is_empty() {
            return None;
        }
        // Check if, on each side, there are exactly one type of unit.
        let unit_and_position = |(idx, pos)| (self.node(idx).unwrap().node.0, pos);
        let is_tail_side_par_unique = {
            let key = unit_and_position(tail_side_par[0]);
            tail_side_par.iter().all(|&e| unit_and_position(e) == key)
        };
        let is_head_side_par_unique = {
            let key = unit_and_position(head_side_par[0]);
            head_side_par.iter().all(|&e| unit_and_position(e) == key)
        };
        if !is_tail_side_par_unique || !is_head_side_par_unique {
            return None;
        }
        // Check if the siblings meets on both side.
        if tail_side_sibs.len() != head_side_sibs.len() {
            return None;
        }
        let tail_iter = tail_side_sibs.iter().map(|&e| unit_and_position(e).0);
        let head_iter = head_side_sibs.iter().map(|&e| unit_and_position(e).0);
        if tail_iter.zip(head_iter).any(|(x, y)| x != y) {
            return None;
        }
        // First determine the one to be retained.

        let mut sibs: Vec<_> = tail_side_sibs.iter().map(|(node, _)| *node).collect();
        sibs.sort_by_cached_key(|&node| self.node(node).unwrap().occ);
        let retain = sibs.pop().unwrap();
        assert!(!sibs.is_empty());
        Some((retain, sibs))
    }
    // Get (node,position), return the reflexitive parents and reflex siblings.
    // In other words, we define the k-th siblings S(k) and the k-th parents P(k) as
    // S(0) = {(node,position)}
    // P(k) = {(n,p) connected to S(k)}
    // S(k) = {(n,p) connected to P(k-1)}
    // And the reflex parent is define as P(\infty) and S(\infty).
    // In practice, for efficiency issue, we use P(cut) and S(cut).
    pub fn get_reflex_nodes(
        &self,
        node: NodeIndex,
        position: Position,
        cut: usize,
    ) -> (GraphBoundary, GraphBoundary) {
        let mut sibs = vec![(node, position)];
        let mut parents: Vec<_> = vec![];
        for _ in 0..cut {
            let par_len = parents.len();
            parents = sibs
                .iter()
                .flat_map(|&(n, p)| {
                    self.edges_from(n, p)
                        .iter()
                        .map(|e| (e.to, e.to_position))
                        .collect::<Vec<_>>()
                })
                .collect();
            parents.sort();
            parents.dedup();
            let sib_size = sibs.len();
            sibs = parents
                .iter()
                .flat_map(|&(n, p)| {
                    self.edges_from(n, p)
                        .iter()
                        .map(|e| (e.to, e.to_position))
                        .collect::<Vec<_>>()
                })
                .collect();
            sibs.sort();
            sibs.dedup();
            if sib_size == sibs.len() || par_len == parents.len() {
                break;
            }
        }
        parents.sort_unstable_by_key(|x| x.0);
        sibs.sort_unstable_by_key(|x| x.0);
        (parents, sibs)
    }
    /// Return the position where a simple path from given edge ends.
    pub fn destination(&self, edge: &DitchEdge) -> (NodeIndex, Position) {
        let mut current_node = edge.to;
        let mut current_pos = edge.to_position;
        loop {
            // Check if there's more than two indegree/outdegree.
            let (indeg, outdeg) =
                self.node(current_node)
                    .unwrap()
                    .edges
                    .iter()
                    .fold((0, 0), |(indeg, outdeg), e| {
                        match e.from_position == current_pos {
                            true => (indeg + 1, outdeg),
                            false => (indeg, outdeg + 1),
                        }
                    });
            assert!(0 < indeg);
            if 1 < indeg {
                break;
            }
            // Move position.
            current_pos = !current_pos;
            if outdeg != 1 {
                break;
            }
            // Move node. Unwrap never panics here.
            let mut edges = self.node(current_node).unwrap().edges.iter();
            let traced_edge = edges.find(|e| e.from_position == current_pos).unwrap();
            current_node = traced_edge.to;
            current_pos = traced_edge.to_position;
        }
        (current_node, current_pos)
    }
    /// Return simple path from the given edge to ends.
    /// Note that a node is added when the node is consumed in the simple path, i.e., we moved
    /// from a position to the opposite position.
    pub fn simple_path_from(&self, edge: &DitchEdge) -> Vec<NodeIndex> {
        let mut current_node = edge.to;
        let mut current_pos = edge.to_position;
        let mut nodes = vec![];
        loop {
            // Check if there's more than two indegree/outdegree.
            let (indeg, outdeg) =
                self.node(current_node)
                    .unwrap()
                    .edges
                    .iter()
                    .fold((0, 0), |(indeg, outdeg), e| {
                        match e.from_position == current_pos {
                            true => (indeg + 1, outdeg),
                            false => (indeg, outdeg + 1),
                        }
                    });
            assert!(0 < indeg);
            if 1 < indeg {
                break;
            }
            // Move position.
            current_pos = !current_pos;
            // Consumed!
            nodes.push(current_node);
            if outdeg != 1 {
                break;
            }
            // Move node. Unwrap never panics here.
            let mut edge = self.node(current_node).unwrap().edges.iter();
            let traced_edge = edge.find(|e| e.from_position == current_pos).unwrap();
            current_node = traced_edge.to;
            current_pos = traced_edge.to_position;
        }
        nodes
    }
    /// Return simple path start from the given node and position, and the destination nodes after this simple path.
    pub fn simple_path_and_dest(
        &self,
        mut node: NodeIndex,
        mut position: Position,
    ) -> (Vec<NodeIndex>, Vec<NodeIndex>) {
        let mut nodes = vec![];
        loop {
            // Move position.
            position = !position;
            nodes.push(node);
            let edges = self.node(node).unwrap().edges.iter();
            let mut edges = edges.filter(|e| e.from_position == position);
            let edge = match edges.next() {
                Some(res) => res,
                None => break,
            };
            // If the out-dgree is more than 1, return.
            // The destination nodes connected to (node,position).
            if edges.next().is_some() {
                break;
            }
            // Check the in-deg of the destination.
            // The destination nodes connected to (node,position).
            let indeg = self.count_edges(edge.to, edge.to_position);
            assert!(1 <= indeg);
            if 1 < indeg {
                break;
            }
            // Move node.
            node = edge.to;
            position = edge.to_position;
        }
        let mut dests: Vec<_> = self
            .node(node)
            .unwrap()
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .map(|e| e.to)
            .collect();
        dests.sort_unstable();
        (nodes, dests)
    }
    /// Remove small tips.
    /// In this function, for all node with the degree of zero,
    /// and its the copy number zero,
    /// we speculate the *local* coverage by traversing `diag` distance.
    /// If the coverage is less than `cov * thr`, the node would be removed.
    pub fn remove_tips(&mut self, thr: f64, diag: usize) {
        let to_remove: HashSet<_> = self
            .nodes()
            .filter(|(_, node)| matches!(node.copy_number, Some(0)))
            .filter_map(|(index, node)| {
                for &pos in &[Position::Head, Position::Tail] {
                    if node.edges.iter().filter(|e| e.from_position == pos).count() == 0 {
                        let coverage = self.local_coverage(index, pos, diag);
                        if (node.occ as f64) < coverage * thr {
                            return Some(index);
                        }
                    }
                }
                None
            })
            .collect();
        self.remove_nodes(&to_remove);
    }
    // compute the local coverage.
    // In other words, it is the average coverage with distance up to diag diameter.
    fn local_coverage(&self, node: NodeIndex, pos: Position, diag: usize) -> f64 {
        let (mut total_cov, mut total_node) = (0, 0);
        let mut current = vec![(node, pos)];
        for _ in 0..diag {
            let mut next_nodes = vec![];
            for &(node, pos) in current.iter() {
                if let Some(node) = self.node(node) {
                    total_cov += node.occ;
                    total_node += node.copy_number.unwrap_or(1);
                    for edge in node.edges.iter().filter(|e| e.from_position == !pos) {
                        next_nodes.push((edge.to, edge.to_position));
                    }
                }
            }
            next_nodes.sort();
            next_nodes.dedup();
            current = next_nodes;
        }
        total_cov as f64 / total_node as f64
    }
    fn remove_nodes(&mut self, to_remove: &HashSet<NodeIndex>) {
        for &idx in to_remove {
            self.delete(idx);
        }
    }
    fn remove_edge_between(
        &mut self,
        (from, from_pos): (NodeIndex, Position),
        (to, to_pos): (NodeIndex, Position),
    ) -> Option<DitchEdge> {
        let removed = self.node_mut(from).and_then(|node| {
            node.edges
                .iter()
                .find(|e| (e.from_position == from_pos && e.to == to && e.to_position == to_pos))
                .cloned()
        });
        if let Some(node) = self.node_mut(from) {
            node.edges.retain(|e| {
                !(e.from_position == from_pos && e.to == to && e.to_position == to_pos)
            });
        }
        // Remove in anti-direction.
        if let Some(node) = self.node_mut(to) {
            node.edges.retain(|e| {
                !(e.from_position == to_pos && e.to == from && e.to_position == from_pos)
            });
        }
        removed
    }
    // Removing `from` node if the copy number is zero, and the copy numbers
    // of the edges from this node are also zero.
    // Then, recursively call this function to all (previously) connected edges.
    fn remove_node_recursive(&'b mut self, from: NodeIndex) {
        if !self.is_deleted(from) {
            let node = self.node(from).unwrap();
            if matches!(node.copy_number, Some(0)) {
                let all_zero = node.edges.iter().all(|e| matches!(e.copy_number, Some(0)));
                if all_zero {
                    let mut affected: Vec<_> = node.edges.iter().map(|e| e.to).collect();
                    affected.sort_unstable();
                    affected.dedup();
                    let remove_edges: Vec<_> = node.edges.iter().map(|e| e.key()).collect();
                    for (f, t) in remove_edges {
                        self.remove_edge_between(f, t);
                    }
                    for node in affected {
                        self.remove_node_recursive(node);
                    }
                    self.delete(from);
                }
            }
        }
    }
    // Removing (from,from_pos)-(to,to_pos) edge.
    // Then, recursively removing the edges in the opposite position and `to` node itself
    // if the degree of `(to,to_pos)` is zero after removing.
    fn remove_edge_and_pruning(
        &'b mut self,
        (from, from_pos): (NodeIndex, Position),
        (to, to_pos): (NodeIndex, Position),
    ) {
        self.remove_edge_between((from, from_pos), (to, to_pos));
        // If the degree became zero and this node and the copy number is zero, remove recursively
        let removed_all_edges = self
            .node(to)
            .map(|node| node.edges.iter().all(|e| e.from_position == !to_pos))
            .unwrap_or(false);
        let to_copy_num = self.node(to).and_then(|n| n.copy_number);
        let is_zero_copy = matches!(to_copy_num, Some(0));
        if removed_all_edges && is_zero_copy {
            // Removing this node.
            // First, recursively call the removing function.
            let remove_edges: Vec<_> = self
                .node(to)
                .unwrap()
                .edges
                .iter()
                .map(|e| (e.to, e.to_position))
                .collect();
            for &(node, pos) in remove_edges.iter() {
                // The edge is from the *opposite* position of the `to` node!
                self.remove_edge_and_pruning((to, !to_pos), (node, pos));
            }
            // Then, removing this node itself.
            if !self.is_deleted(to) {
                self.delete(to);
            }
        }
    }
    /// Remove small connected components with the size less than `thr`.
    pub fn remove_small_component(&'b mut self, thr: usize) {
        let mut to_remove = HashSet::new();
        let mut arrived = HashSet::new();
        let mut cluster_num = 0;
        for (idx, _node) in self.nodes() {
            if arrived.contains(&idx) {
                continue;
            }
            let mut actives = HashSet::new();
            let mut stack = vec![idx];
            while !stack.is_empty() {
                let last = stack.pop().unwrap();
                actives.insert(last);
                arrived.insert(last);
                for to in self.node(last).unwrap().edges.iter().map(|x| x.to) {
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
    pub fn collapse_bubble(&'b mut self, c: &AssembleConfig) {
        let mut to_remove = HashSet::new();
        let mut queue = std::collections::VecDeque::new();
        for (index, _) in self.nodes() {
            queue.push_back((index, Position::Head));
            queue.push_back((index, Position::Tail));
        }
        while let Some((index, position)) = queue.pop_front() {
            let edges: Vec<_> = self.edges_from(index, position);
            if edges.len() <= 1 {
                continue;
            }
            let pos = edges[0].to_position;
            let unit = edges[0].to.0;
            let is_the_same_unit = edges.iter().all(|e| e.to_position == pos && e.to.0 == unit);
            let index_is_the_unique_parent = edges.iter().all(|e| {
                self.node(e.to)
                    .unwrap()
                    .edges
                    .iter()
                    .all(|f| (f.from_position != pos) | (f.to == index))
            });
            if is_the_same_unit && index_is_the_unique_parent {
                let (new_terminal, removed_nodes) = self.collapse_bubble_from(index, position, c);
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
    // Returned the index of the merged node, the opposite position, and the nodes merged.
    fn collapse_bubble_from(
        &'b mut self,
        root: NodeIndex,
        position: Position,
        _c: &AssembleConfig,
    ) -> ((NodeIndex, Position), Vec<NodeIndex>) {
        assert!(self.sanity_check());
        // Check the collapsing condition.
        let edges: Vec<_> = self.edges_from(root, position);
        assert!(edges.len() > 1);
        for edge in edges.iter() {
            assert_eq!(1, self.count_edges(edge.to, edge.to_position));
            let first = self
                .node(edge.to)
                .unwrap()
                .edges
                .iter()
                .find(|e| e.from_position == edge.to_position)
                .unwrap();
            assert_eq!(first.to, root);
            assert_eq!(first.to_position, position);
        }
        let mut edges: Vec<_> = edges
            .into_iter()
            .map(|e| (e.to, e.to_position, e.occ))
            .collect();
        edges.sort_by_key(|&(_, _, occ)| occ);
        // Merge all non-primary edges into the rist element. What is to merge is just occurance.
        let total_occ: usize = edges.iter().map(|&(_, _, occ)| occ).sum();
        let (primary, primary_pos, _) = edges.pop().unwrap();
        // Root -> primary edge.
        let root_node = self.node_mut(root).unwrap();
        root_node.edges.retain(|e| {
            e.from_position != position || (e.to, e.to_position) == (primary, primary_pos)
        });
        root_node
            .edges
            .iter_mut()
            .find(|e| (e.to, e.to_position) == (primary, primary_pos))
            .unwrap()
            .occ = total_occ;
        // Primary edge -> root.
        let primary_node = self.node_mut(primary).unwrap();
        // We do not need to remove any edges here. Just update.
        primary_node
            .edges
            .iter_mut()
            .find(|e| (e.to, e.to_position) == (root, position))
            .unwrap()
            .occ = total_occ;
        // Remove non-primary.
        let mut merged_edges: Vec<DitchEdge> = vec![];
        let mut removed_node_indices = vec![];
        for &(sec, sec_pos, _) in edges.iter() {
            let sec_node = self.node(sec).unwrap();
            let sec_edges = sec_node
                .edges
                .iter()
                .filter(|e| e.from_position == sec_pos)
                .filter(|e| (e.to, e.to_position) != (root, position))
                .cloned();
            merged_edges.extend(sec_edges);
            removed_node_indices.push(sec);
            self.delete(sec);
        }
        for edge in merged_edges {
            // Sib -> primary.
            let sib_node = self.node_mut(edge.to).unwrap();
            let removed_idx = sib_node
                .edges
                .iter()
                .position(|e| (e.to, e.to_position) == (edge.from, edge.from_position))
                .unwrap();
            let mut removed_edge = sib_node.edges.remove(removed_idx);
            removed_edge.to = primary;
            removed_edge.to_position = primary_pos;
            let mut sib_edges = sib_node.edges.iter_mut();
            match sib_edges.find(|e| (e.to, e.to_position) == (primary, !primary_pos)) {
                Some(e) => e.occ += edge.occ,
                None => sib_node.edges.push(removed_edge),
            }
            // Primary->Sib.
            let primary_node = self.node_mut(primary).unwrap();
            let mut primary_edges = primary_node.edges.iter_mut();
            match primary_edges.find(|e| (e.to, e.to_position) == (edge.to, edge.to_position)) {
                Some(e) => e.occ += edge.occ,
                None => {
                    let mut new_edge = edge.reverse();
                    new_edge.from = primary;
                    new_edge.from_node = primary_node.node;
                    primary_node.edges.push(new_edge)
                }
            }
        }
        ((primary, !primary_pos), removed_node_indices)
    }
    /// Squish short bubbles.
    /// In other words, if there is a branch,
    /// where the children have the same length,
    /// and the same content w.r.t their unit IDs,
    /// then it squish these children into one contig.
    /// Note that this proc is *not* bubble collapsing.
    /// For example, it squish the small contig B,C, as follows:
    ///          --|--B--|----- -|--D--|-- ...
    ///        /              /
    /// ---A---              /
    ///        \___|__C__|__/__|__E__|___ ...
    /// (Be careful. The B contig is connecting to D, not E, whereas
    /// C is connecting to both D and E.
    /// The function returns how to change the **clustering** on each unit.
    pub fn squish_bubbles(&self, len: usize) -> HashMap<Node, u64> {
        let mut squish_to: HashMap<Node, u64> = HashMap::new();
        for (_index, node) in self.nodes() {
            for pos in [Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                if edges.clone().count() <= 1 {
                    continue;
                }
                let path_and_dest: Vec<_> = edges
                    .map(|e| self.simple_path_and_dest(e.to, e.to_position))
                    .collect();
                let is_bubble = {
                    let paths: Vec<_> = path_and_dest
                        .iter()
                        .map(|x| {
                            let mut units: Vec<_> = x.0.iter().map(|x| x.0).collect();
                            units.sort_unstable();
                            units
                        })
                        .collect();
                    paths.iter().all(|us| us == &paths[0] && us.len() <= len)
                };
                if is_bubble {
                    let mut convert_table: HashMap<u64, u64> = HashMap::new();
                    for (path, _) in path_and_dest.iter() {
                        for &node in path.iter() {
                            let (unit, cluster) = self.node(node).unwrap().node;
                            convert_table
                                .entry(unit)
                                .and_modify(|x| *x = (*x).min(cluster))
                                .or_insert(cluster);
                        }
                    }
                    for (path, _) in path_and_dest.iter() {
                        for &node in path.iter() {
                            let (unit, cluster) = self.node(node).unwrap().node;
                            squish_to
                                .entry((unit, cluster))
                                .and_modify(|to| *to = (*to).min(convert_table[&unit]))
                                .or_insert_with(|| convert_table[&unit]);
                        }
                    }
                }
            }
        }
        squish_to
    }
    /// Z-selection of edges.
    /// In this function, we select edges based on the topology of the graph.
    /// Suppose a node with the degree more than 2 has a edge with no conflicting edge.
    /// Then, we remove edges from that node satisfying the conditions below:
    /// 1. It connect to a node with degree more than 1.
    /// 2. The connected node has at least one connected node with outdegree is equal to 1.
    /// In other words, we remove edges we could not select, otherwise the resulting graph would be split.
    pub fn z_edge_selection(&mut self) {
        let mut removed_edges = HashSet::new();
        let mut retain_edges = HashSet::new();
        let format_edge = |e: &DitchEdge| {
            let (f, t) = e.key();
            (f.min(t), f.max(t))
        };
        for (_index, node) in self.nodes() {
            for &pos in &[Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                let edges_selectable = edges.clone().map(|edge| self.can_select(edge));
                // If there's selectable edges, remove other edges.
                let num_of_selectable = edges_selectable.clone().filter(|&x| x).count();
                let num_of_removable = edges_selectable.clone().count() - num_of_selectable;
                if 0 < num_of_selectable && 0 < num_of_removable {
                    for (edge, selectable) in edges.zip(edges_selectable) {
                        match selectable {
                            true => retain_edges.insert(format_edge(edge)),
                            false => removed_edges.insert(format_edge(edge)),
                        };
                    }
                };
            }
        }
        for retain in retain_edges.iter() {
            removed_edges.remove(retain);
        }
        self.modify_by(|node| {
            node.edges
                .retain(|edge| !removed_edges.contains(&format_edge(edge)));
        });
    }
    // fn share_units(&self, boundary: &GraphBoundary) -> bool {
    //     let boundary: Vec<_> = boundary
    //         .iter()
    //         .map(|(n, pos)| (self.node(*n).unwrap().node.0, pos))
    //         .collect();
    //     if boundary.is_empty() {
    //         false
    //     } else {
    //         boundary.iter().all(|&elm| elm == boundary[0])
    //     }
    // }
    // fn z_edge_squish(&mut self) {
    //     let mut keys: Vec<_> = self.nodes().map(|n| n.0).collect();
    //     keys.sort_unstable();
    //     for node in keys {
    //         if self.is_deleted(node) || self.has_self_loop(node) {
    //             continue;
    //         }
    //         for position in [Position::Head, Position::Tail] {
    //             if self.count_edges(node, position) <= 1 {
    //                 continue;
    //             }
    //             let (parents, sibs) = self.get_reflex_nodes(node, position, 6);
    //             let par_copy = parents
    //                 .iter()
    //                 .filter_map(|&(n, _)| self.node(n).and_then(|n| n.copy_number))
    //                 .all(|cp| cp == 1);
    //             let sib_copy = sibs
    //                 .iter()
    //                 .filter_map(|&(n, _)| self.node(n).and_then(|n| n.copy_number))
    //                 .all(|cp| cp == 1);
    //             if !(par_copy && sib_copy) {
    //                 continue;
    //             }
    //             if !(self.share_units(&parents) && self.share_units(&sibs)) {
    //                 continue;
    //             }
    //             let sibs: Vec<_> = sibs.into_iter().filter(|&x| x.0 != node).collect();
    //             let retain = node;
    //             drop(node);
    //             // Make all the edges into sibs to retain.
    //             let (edges, increase_occ, increase_copy_num) = {
    //                 let (mut edges, mut occ, mut cp) = (vec![], 0, 0);
    //                 for &(node, _) in sibs.iter() {
    //                     let removed = self.node(node).unwrap();
    //                     edges.extend(removed.edges.clone());
    //                     occ += removed.occ;
    //                     cp += removed.copy_number.unwrap_or(0);
    //                     self.delete(node);
    //                 }
    //                 let retain_node = self.node(retain).unwrap().node;
    //                 for edge in edges.iter_mut() {
    //                     edge.from = retain;
    //                     edge.from_node = retain_node;
    //                 }
    //                 (edges, occ, cp)
    //             };
    //             {
    //                 let retain_node = self.node_mut(retain).unwrap();
    //                 retain_node.occ += increase_occ;
    //                 if let Some(cp) = retain_node.copy_number.as_mut() {
    //                     *cp += increase_copy_num;
    //                 }
    //             };
    //             for edge in edges {
    //                 if self.has_edge(edge.key()) {
    //                     self.merge_edge(&edge);
    //                 } else {
    //                     self.add_edge(edge);
    //                 }
    //             }
    //         }
    //     }
    // }
    /// Check if we can select this edge when we need to select an edge for each
    /// node while preserving the connectiviy of the graph.
    /// true if we could select this edge.
    pub fn can_select(&self, edge: &DitchEdge) -> bool {
        let to = self.node(edge.to).unwrap();
        // Filtering out the back edge.
        for to_edge in to
            .edges
            .iter()
            .filter(|to_e| to_e.from_position == edge.to_position && to_e.to != edge.from)
        {
            let sibling = self.node(to_edge.to).unwrap();
            // Check if `child.to` has only one parent, `to`.
            assert_eq!(to_edge.from, edge.to);
            assert_eq!(to_edge.from_position, edge.to_position);
            let parent = (to_edge.from, to_edge.from_position);
            let to_is_only_prent = sibling
                .edges
                .iter()
                .filter(|c_edge| c_edge.from_position == to_edge.to_position)
                .all(|c_edge| (c_edge.to, c_edge.to_position) == parent);
            if to_is_only_prent {
                return false;
            }
        }
        true
    }
    /// Remove lightweight edges with occurence less than or equal to `thr`.
    /// To retain that edge if the edge is the only edge from its terminal,
    /// set `retain_single_edge` to `true`.
    pub fn remove_lightweight_edges(&mut self, thr: usize, retain_single_edge: bool) {
        debug!("RM\t{thr}");
        let mut removed_edges = vec![];
        for (from, node) in self.nodes() {
            for &pos in &[Position::Head, Position::Tail] {
                if self.count_edges(from, pos) <= 1 {
                    continue;
                }
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                for e in edges.filter(|e| e.occ <= thr) {
                    // If retain mode, check this is not the only edge
                    let is_safe = self
                        .node(e.to)
                        .unwrap()
                        .edges
                        .iter()
                        .filter(|f| f.from_position == e.to_position)
                        .any(|f| thr < f.occ);
                    let removable = !retain_single_edge || is_safe;
                    if removable {
                        let (f, t) = e.key();
                        let norm_key = (f.min(t), f.max(t));
                        removed_edges.push(norm_key);
                    }
                }
            }
        }
        removed_edges.sort_unstable();
        removed_edges.dedup();
        for (f, t) in removed_edges {
            self.remove_edge_between(f, t);
        }
    }
}

// fn get_bridge(nodes: usize, edges: &[Vec<usize>]) -> Vec<Vec<usize>> {
//     warn!("Have you tested this code? I think not.");
//     let mut is_arrived = vec![false; nodes];
//     let mut bridges = vec![vec![]; nodes];
//     let mut order = vec![-1; nodes];
//     let mut low = vec![-1; nodes];
//     let mut parent = vec![0; nodes];
//     let mut count = 0;
//     for i in 0..nodes {
//         if is_arrived[i] {
//             continue;
//         }
//         let mut stack = vec![i];
//         'dfs: while !stack.is_empty() {
//             let last = *stack.last().unwrap();
//             if !is_arrived[last] {
//                 is_arrived[last] = true;
//                 order[last] = count as i64;
//                 count += 1;
//             }
//             for &to in edges[last].iter() {
//                 if !is_arrived[to] {
//                     parent[to] = last;
//                     stack.push(to);
//                     continue 'dfs;
//                 }
//             }
//             // Postorder, we have arrived all the nodes below this node.
//             let last = stack.pop().unwrap();
//             for &to in edges[last].iter() {
//                 if parent[to] == last {
//                     low[last] = low[last].min(low[to]);
//                     if order[last] < low[to] {
//                         bridges[last].push(to);
//                         bridges[to].push(last);
//                     }
//                 } else {
//                     // This is a back edge.
//                     low[last] = low[last].min(order[to]);
//                 }
//             }
//         }
//     }
//     bridges
// }

#[cfg(test)]
mod tests {
    use super::*;
    use definitions::EncodedRead;
    use definitions::ReadType;
    use rand::Rng;
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128StarStar;
    fn gen_read<R: Rng>(id: u64, rng: &mut R, hap: &[u64]) -> EncodedRead {
        let original_length = 5_000 + rng.gen_range(0..10_000);
        let start_pos = rng.gen_range(0..hap.len());
        let seq = vec![b'A'; 2_000];
        let cl = 2;
        let nodes: Vec<_> = hap
            .iter()
            .cycle()
            .skip(start_pos)
            .take(original_length / 2_000)
            .enumerate()
            .map(|(idx, &unit)| {
                let position = idx as usize * 2_000;
                let cigar = vec![definitions::Op::Match(2_000)];
                definitions::Node::new(unit, true, seq.clone(), cigar, position, cl)
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
                let label = Vec::new();
                definitions::Edge {
                    from: from.unit,
                    to: to.unit,
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
        let units: Vec<_> = node_cp
            .iter()
            .enumerate()
            .map(|(id, &cp)| Unit::new(id as u64, vec![], cp))
            .collect();
        let hap: Vec<_> = vec![0, 1, 2, 4, 3, 2, 6, 5, 2, 6, 5, 2, 7, 8];
        let read_num = 2 * 2_000 * 30 * hap.len() / 10_000;
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let reads: Vec<_> = (0..read_num)
            .map(|i| gen_read(i as u64, &mut rng, &hap))
            .collect();
        let total_units: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let cov = (total_units / hap.len() / 2) as f64;
        // let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let assemble_config = AssembleConfig::new(100, false, false, 6, 1f64);
        let graph = DitchGraph::new(&reads, &units, ReadType::CCS, &assemble_config);
        let (nodes, _) = graph.copy_number_estimation_gbs(cov);
        for (i, &cp) in node_cp.iter().enumerate() {
            let node = (i as u64, 0);
            let index = graph.nodes_index[&node];
            assert_eq!(cp, nodes[index]);
        }
    }
    #[test]
    fn from_reads_2() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let units: Vec<_> = node_cp
            .iter()
            .enumerate()
            .map(|(id, &cp)| Unit::new(id as u64, vec![], cp))
            .collect();
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
                    .filter(|w| w[0].unit == 1 && w[1].unit == 3)
                    .count()
            })
            .sum();
        println!("(1,3)\t{}", count);
        let count: usize = reads
            .iter()
            .map(|r| {
                r.nodes
                    .windows(2)
                    .filter(|w| w[0].unit == 1 && w[1].unit == 2)
                    .count()
            })
            .sum();
        println!("(1,2)\t{}", count);
        let total_units: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let cov = (total_units / (hap1.len() + hap2.len())) as f64;
        // let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let assemble_config = AssembleConfig::new(100, false, false, 6, 1f64);
        let graph = DitchGraph::new(&reads, &units, ReadType::CCS, &assemble_config);
        let (nodes, _) = graph.copy_number_estimation_gbs(cov);
        for (i, &cp) in node_cp.iter().enumerate() {
            let index = graph.nodes_index[&(i as u64, 0)];
            assert_eq!(cp, nodes[index]);
        }
    }
    #[test]
    fn from_reads_2_mst() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
        let units: Vec<_> = node_cp
            .iter()
            .enumerate()
            .map(|(id, &cp)| Unit::new(id as u64, vec![], cp))
            .collect();
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
                    .filter(|w| w[0].unit == 1 && w[1].unit == 3)
                    .count()
            })
            .sum();
        println!("(1,3)\t{}", count);
        let count: usize = reads
            .iter()
            .map(|r| {
                r.nodes
                    .windows(2)
                    .filter(|w| w[0].unit == 1 && w[1].unit == 2)
                    .count()
            })
            .sum();
        println!("(1,2)\t{}", count);
        let total_units: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let cov = (total_units / (hap1.len() + hap2.len())) as f64;
        let assemble_config = AssembleConfig::new(100, false, false, 6, 1f64);
        let mut graph = DitchGraph::new(&reads, &units, ReadType::CCS, &assemble_config);
        assert!(graph.sanity_check());
        println!("graph:{graph:?}");
        graph.assign_copy_number_mst(cov, &mut rng);
        println!("{node_cp:?}");
        let mut nodes: Vec<_> = graph
            .nodes()
            .map(|(_, n)| (n.node.0, n.copy_number))
            .collect();
        nodes.sort_unstable_by_key(|x| x.0);
        for (node, cp) in nodes {
            assert_eq!(cp.unwrap(), node_cp[node as usize]);
        }
    }
}
