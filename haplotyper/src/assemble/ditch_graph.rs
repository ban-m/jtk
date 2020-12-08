use super::AssembleConfig;
use crate::find_union;
use definitions::{Edge, EncodedRead};
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;
const DRAFT_COVERAGE: usize = 20;
const DRAFT_REP_NUM: usize = 7;
/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
#[derive(Clone)]
pub struct DitchGraph<'a, 'b, 'c> {
    nodes: Vec<DitchNode<'a, 'b, 'c>>,
    index: HashMap<NodeIndex, usize>,
}

impl<'a, 'b, 'c> std::fmt::Display for DitchGraph<'a, 'b, 'c> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        use histgram_viz::Histgram;
        let nodes = self.nodes.len();
        let edges = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        writeln!(f, "Node:{}, Edges:{}", nodes, edges)?;
        let occs: Vec<_> = self.nodes.iter().map(|n| n.nodes.len()).collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Node Occs:{}", hist.format(20, 20))?;
        let occs: Vec<_> = self
            .nodes
            .iter()
            .flat_map(|n| n.edges.iter().map(|e| e.edges.len()))
            .collect();
        let hist = Histgram::new(&occs);
        writeln!(f, "Edge Occs:{}", hist.format(20, 20))?;
        let degrees = {
            let mut degs: HashMap<usize, usize> = HashMap::new();
            for node in self.nodes.iter() {
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

impl<'a, 'b, 'c> std::fmt::Debug for DitchGraph<'a, 'b, 'c> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edge = self.nodes.iter().map(|e| e.edges.len()).sum::<usize>();
        writeln!(f, "Nodes:{}\tEdges:{}\n", self.nodes.len(), edge)?;
        let lines: Vec<_> = self
            .nodes
            .iter()
            .map(|node| {
                let idx = self
                    .index
                    .get(&NodeIndex::new(node.unit, node.cluster))
                    .unwrap();
                format!("{}-{}", idx, node)
            })
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
pub struct DitchNode<'a, 'b, 'c> {
    unit: u64,
    cluster: u64,
    // CAUTION!!!!!! The `unit` and `cluster` members should not be look-upped!
    nodes: Vec<&'a definitions::Node>,
    edges: Vec<DitchEdge<'b>>,
    tips: Vec<DitchTip<'c>>,
}

impl<'a, 'b, 'c> std::fmt::Display for DitchNode<'a, 'b, 'c> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "---------{}:{}---------", self.unit, self.cluster)?;
        writeln!(f, "Seq:")?;
        for (idx, n) in self.nodes.iter().enumerate() {
            writeln!(f, "{:3}:{}:{}", idx, n.seq, n.is_forward)?;
        }
        let lines: Vec<_> = self.edges.iter().map(|e| format!("{}", e)).collect();
        write!(f, "Edges:\n{}", lines.join("\n"))
    }
}

impl<'a, 'b, 'c> DitchNode<'a, 'b, 'c> {
    fn new(unit: u64, cluster: u64) -> Self {
        Self {
            unit,
            cluster,
            nodes: vec![],
            edges: vec![],
            tips: vec![],
        }
    }
    fn seq(&self) -> String {
        self.nodes
            .get(0)
            .map(|node| node.seq.clone())
            .unwrap_or_else(String::new)
    }
    fn consensus(&self) -> String {
        let result = self.consensus_with(100, DRAFT_REP_NUM);
        String::from_utf8(result).unwrap()
    }
    fn consensus_with(&self, len: usize, num: usize) -> Vec<u8> {
        let chunks: Vec<_> = self
            .nodes
            .iter()
            .take(DRAFT_COVERAGE)
            .map(|n| crate::local_clustering::node_to_subchunks(n, len))
            .collect();
        let max_pos = chunks
            .iter()
            .flat_map(|cs| cs.iter().map(|c| c.pos).max())
            .max()
            .unwrap();
        let mut subseqs = vec![vec![]; max_pos + 1];
        for cs in chunks.iter() {
            for c in cs.iter() {
                subseqs[c.pos].push(c.seq.as_slice());
            }
        }
        subseqs
            .par_iter()
            .flat_map(|subseq| consensus(subseq, num))
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct DitchEdge<'a> {
    from: usize,
    to: usize,
    from_position: Position,
    to_position: Position,
    // If true, the edge is `forward` one.
    // In other words, if true,
    // you can spell this edge by its label,
    // otherwise you should take rev-cmp of the label.
    edges: Vec<(&'a definitions::Edge, bool)>,
}

impl<'a> std::fmt::Display for DitchEdge<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}:{})->", self.from, self.from_position)?;
        write!(f, "({}:{})", self.to, self.to_position)?;
        let edgs: Vec<_> = self
            .edges
            .iter()
            .map(|x| format!("{}", x.0.label.len()))
            .collect();
        write!(f, " {} threads[{}].", edgs.len(), edgs.join(","))
    }
}

impl<'a> std::cmp::PartialEq for DitchEdge<'a> {
    fn eq<'b>(&self, other: &DitchEdge<'b>) -> bool {
        self.from == other.from
            && self.to == other.to
            && self.from_position == other.from_position
            && self.to_position == other.to_position
    }
}

impl<'a> std::cmp::Eq for DitchEdge<'a> {}

impl<'a> DitchEdge<'a> {
    fn new(from: Position, from_i: usize, to: Position, to_i: usize) -> Self {
        Self {
            from: from_i,
            to: to_i,
            from_position: from,
            to_position: to,
            edges: vec![],
        }
    }
    fn push(&mut self, x: (&'a Edge, bool)) {
        self.edges.push(x);
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

#[derive(Debug, Clone, Eq, PartialEq, Copy, Hash)]
pub struct NodeIndex {
    unit: u64,
    cluster: u64,
}
impl NodeIndex {
    pub fn new(unit: u64, cluster: u64) -> Self {
        Self { unit, cluster }
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
    None,
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
}

impl<'a> DitchGraph<'a, 'a, 'a> {
    pub fn new(reads: &[&'a EncodedRead], c: &AssembleConfig) -> Self {
        // Allocate all nodes.
        let mut index = HashMap::new();
        let mut nodes = vec![];
        for node in reads.iter().flat_map(|r| r.nodes.iter()) {
            let node_index = NodeIndex::new(node.unit, node.cluster);
            index.entry(node_index).or_insert_with(|| {
                nodes.push(DitchNode::new(node.unit, node.cluster));
                nodes.len() - 1
            });
            let node_index = *index.get(&node_index).unwrap();
            nodes[node_index].nodes.push(node);
        }
        let mut graph = Self { index, nodes };
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
            let index = *self.index.get(&NodeIndex::new(first.unit, first.cluster))?;
            let direction = true;
            let position = if first.is_forward { Head } else { Tail };
            self.nodes[index]
                .tips
                .push(DitchTip::new(&read.leading_gap, position, direction));
        }
        if let Some(last) = read.nodes.last() {
            let index = *self.index.get(&NodeIndex::new(last.unit, last.cluster))?;
            let direction = true;
            let position = if last.is_forward { Head } else { Tail };
            self.nodes[index]
                .tips
                .push(DitchTip::new(&read.trailing_gap, position, direction));
        }
        Some(())
    }
    fn append_read(&mut self, read: &'a EncodedRead, _c: &AssembleConfig) -> Option<()> {
        for (pairs, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
            let (from, to) = (&pairs[0], &pairs[1]);
            let from_index = *self.index.get(&NodeIndex::new(from.unit, from.cluster))?;
            let to_index = *self.index.get(&NodeIndex::new(to.unit, to.cluster))?;
            use Position::*;
            let from_pos = if from.is_forward { Tail } else { Head };
            let to_pos = if to.is_forward { Head } else { Tail };
            let mut dedg = DitchEdge::new(from_pos, from_index, to_pos, to_index);
            match self.nodes[from_index]
                .edges
                .iter_mut()
                .find(|x| x == &&dedg)
            {
                Some(res) => res.push((edge, true)),
                None => {
                    dedg.push((edge, true));
                    self.nodes[from_index].edges.push(dedg);
                }
            }
            let mut dedg = DitchEdge::new(to_pos, to_index, from_pos, from_index);
            match self.nodes[to_index].edges.iter_mut().find(|x| x == &&dedg) {
                Some(res) => res.push((edge, false)),
                None => {
                    dedg.push((edge, false));
                    self.nodes[to_index].edges.push(dedg);
                }
            }
        }
        Some(())
    }
}

impl<'a, 'b, 'c> DitchGraph<'a, 'b, 'c> {
    fn enumerate_candidates(&self) -> Vec<(usize, Position)> {
        let mut selected = vec![false; self.nodes.len()];
        let mut primary_candidates: Vec<_> = (0..self.nodes.len())
            .filter_map(|i| {
                for &position in &[Position::Head, Position::Tail] {
                    if self.nodes[i]
                        .edges
                        .iter()
                        .filter(|e| e.from_position == position)
                        .count()
                        != 1
                    {
                        selected[i] = true;
                        return Some((i, position));
                    }
                }
                None
            })
            .collect();
        let secondary_candidates =
            (0..self.nodes.len())
                .filter(|&i| !selected[i])
                .filter_map(|i| {
                    for &position in &[Position::Head, Position::Tail] {
                        // The unwrap() never panics.
                        let grand_child = self.nodes[i]
                            .edges
                            .iter()
                            .find(|e| e.from_position == position)
                            .map(|e| {
                                let count = self.nodes[e.to]
                                    .edges
                                    .iter()
                                    .filter(|f| f.from_position == e.to_position)
                                    .count();
                                if count == 0 {
                                    debug!("{},{}", i, position);
                                    debug!("Edge:{}", e);
                                    for edge in &self.nodes[e.to].edges {
                                        debug!("Child:{}", edge)
                                    }
                                    panic!()
                                }
                                count
                            })
                            .unwrap();
                        if grand_child > 1 {
                            return Some((i, position));
                        }
                    }
                    None
                });
        primary_candidates.extend(secondary_candidates);
        primary_candidates
    }
    // Assemble the ditch graph.
    // In other words, it reduce the simple path, currently.
    pub fn spell(
        &self,
        c: &AssembleConfig,
        cl: usize,
    ) -> (
        Vec<gfa::Segment>,
        Vec<gfa::Edge>,
        gfa::Group,
        Vec<ContigSummary>,
    ) {
        let mut arrived = vec![false; self.nodes.len()];
        let mut sids = vec![ContigTag::None; self.nodes.len()];
        let (mut g_segs, mut g_edges, mut summaries) = (vec![], vec![], vec![]);
        let candidates = self.enumerate_candidates();
        for (i, p) in candidates {
            if arrived[i] {
                continue;
            }
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges, summary) =
                self.traverse_from(&mut arrived, &mut sids, i, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
        }
        for i in 0..self.nodes.len() {
            if arrived[i] {
                continue;
            }
            let p = Position::Head;
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges, summary) =
                self.traverse_from(&mut arrived, &mut sids, i, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
        }
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.eid.clone()))
            .collect();
        let uid = Some(format!("group-{}", cl));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group, summaries)
    }
    pub fn remove_tips(&mut self) {
        let sum = self.nodes.iter().map(|e| e.nodes.len()).sum::<usize>();
        let mean = sum / self.nodes.len();
        let thr = mean / 8;
        debug!("Removing nodes less than {} occ.", thr);
        let to_remove: Vec<_> = self.nodes.iter().map(|e| e.nodes.len() < thr).collect();
        assert_eq!(*self.index.values().max().unwrap() + 1, self.nodes.len());
        self.remove_nodes(&to_remove);
        assert_eq!(*self.index.values().max().unwrap() + 1, self.nodes.len());
    }
    // Removing nodes and re-map all of the indices.
    fn remove_nodes(&mut self, to_remove: &[bool]) {
        assert_eq!(self.nodes.len(), to_remove.len());
        assert_eq!(*self.index.values().max().unwrap() + 1, self.nodes.len());
        let mapping = {
            let mut mapping = vec![];
            let mut index = 0;
            for &b in to_remove.iter() {
                mapping.push(index);
                index += !b as usize;
            }
            mapping
        };
        self.index.retain(|_, v| {
            if to_remove[*v] {
                false
            } else {
                *v = mapping[*v];
                true
            }
        });
        let mut idx = 0;
        self.nodes.retain(|_| {
            idx += 1;
            !to_remove[idx - 1]
        });
        assert_eq!(self.index.len(), self.nodes.len());
        self.nodes.iter_mut().for_each(|node| {
            node.edges.retain(|e| !to_remove[e.to]);
            node.edges.iter_mut().for_each(|e| {
                e.from = mapping[e.from];
                e.to = mapping[e.to];
            });
        });
        for node in self.nodes.iter() {
            for edge in node.edges.iter() {
                assert!(self.nodes[edge.to].edges.iter().any(|f| f.to == edge.from));
            }
        }
    }
    /// Resolve repetitive units. (TODO.)
    pub fn resolve_repeats(&self) {}
    pub fn remove_redundant_edges(&mut self, thr: usize) {
        let mut removed_edges = vec![vec![]; self.nodes.len()];
        for (from, node) in self.nodes.iter().enumerate() {
            for &pos in &[Position::Head, Position::Tail] {
                let to_cands: HashSet<_> = node
                    .edges
                    .iter()
                    .filter(|e| e.from_position == pos)
                    .map(|e| (e.to_position, e.to))
                    .collect();
                if to_cands.len() > 1 {
                    for e in node
                        .edges
                        .iter()
                        .filter(|e| e.from_position == pos && e.edges.len() <= thr)
                    {
                        let to = e.to;
                        let t_pos = e.to_position;
                        removed_edges[from].push((pos, to, t_pos));
                        removed_edges[to].push((t_pos, from, pos));
                    }
                }
            }
        }
        self.nodes.iter_mut().enumerate().for_each(|(idx, n)| {
            let to_be_removed = &removed_edges[idx];
            n.edges.retain(|x| {
                let probe = (x.from_position, x.to, x.to_position);
                !to_be_removed.contains(&probe)
            })
        });
    }
    pub fn remove_small_component(&mut self, thr: usize) {
        let mut to_remove = vec![false; self.nodes.len()];
        use find_union::FindUnion;
        let mut cluster = FindUnion::new(self.nodes.len());
        for node in self.nodes.iter() {
            for edge in node.edges.iter() {
                cluster.unite(edge.from, edge.to);
            }
        }
        for node in 0..self.nodes.len() {
            if cluster.find(node).unwrap() != node {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| cluster.find(x).unwrap() == node)
                .count();
            if count < thr {
                for (component, to_remove) in to_remove.iter_mut().enumerate() {
                    if cluster.find(component).unwrap() == node {
                        *to_remove = true;
                    }
                }
            }
        }
        self.remove_nodes(&to_remove);
    }
    pub fn collapse_buddle(&mut self, c: &AssembleConfig) {
        let mut to_remove = vec![false; self.nodes.len()];
        let mut queue = std::collections::VecDeque::new();
        for idx in 0..self.nodes.len() {
            queue.push_back((idx, Position::Head));
            queue.push_back((idx, Position::Tail));
        }
        while let Some((idx, position)) = queue.pop_front() {
            let edges: Vec<_> = self.nodes[idx]
                .edges
                .iter()
                .filter(|e| e.from_position == position && e.from == idx)
                .collect();
            if edges.len() <= 1 {
                continue;
            }
            let pos = edges[0].to_position;
            let unit = self.nodes[edges[0].to].unit;
            let is_the_same_unit = edges
                .iter()
                .all(|e| e.to_position == pos && unit == self.nodes[e.to].unit);
            let idx_is_the_unique_parent = edges.iter().all(|e| {
                self.nodes[e.to]
                    .edges
                    .iter()
                    .filter(|f| f.from_position == pos)
                    .all(|f| f.to == idx)
            });
            if is_the_same_unit && idx_is_the_unique_parent {
                let (new_terminal, removed_nodes) = self.collapse_bubble_from(idx, position, c);
                for i in removed_nodes {
                    to_remove[i] = true;
                }
                queue.push_back(new_terminal)
            }
        }
        self.remove_nodes(&to_remove);
    }
    // Collapse bubble from the given root.
    fn collapse_bubble_from(
        &mut self,
        root: usize,
        position: Position,
        _c: &AssembleConfig,
    ) -> ((usize, Position), Vec<usize>) {
        // Check the collapsing condition.
        let edges: Vec<_> = self.nodes[root]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .collect();
        // Check.
        assert!(edges.len() > 1);
        let (first_id, first_pos, first_unit, first_cluster) = {
            let first = edges.first().unwrap();
            let first_unit = self.nodes[first.to].unit;
            let first_cluster = self.nodes[first.to].cluster;
            (first.to, first.to_position, first_unit, first_cluster)
        };
        assert!(edges.iter().all(|e| e.to_position == first_pos));
        assert!(edges.iter().all(|e| self.nodes[e.to].unit == first_unit));
        let merged_nodes_ids: Vec<_> = edges.iter().map(|e| e.to).collect();
        let mut new_node = DitchNode::new(first_unit, first_cluster);
        let mut edge_to_root = DitchEdge::new(first_pos, first_id, position, root);
        // All of these edges should start from (first_id, !first_pos).
        let mut edges_otherside: HashMap<(usize, Position), Vec<_>> = HashMap::new();
        let node_otherside: HashSet<_> = merged_nodes_ids
            .iter()
            .flat_map(|&idx| self.nodes[idx].edges.iter())
            .filter(|e| !(e.to == root && e.to_position == position))
            .map(|e| (e.to, e.to_position))
            .collect();
        for &node in merged_nodes_ids.iter() {
            new_node.nodes.append(&mut self.nodes[node].nodes);
            while let Some(edge) = self.nodes[node].edges.pop() {
                if edge.to == root {
                    assert!(edge.to_position == position);
                    assert!(edge.from_position == first_pos);
                    edge_to_root.edges.extend(edge.edges);
                } else {
                    assert!(edge.from_position == !first_pos);
                    edges_otherside
                        .entry((edge.to, edge.to_position))
                        .or_default()
                        .extend(edge.edges);
                }
            }
        }
        new_node.edges.push(edge_to_root);
        for ((to, to_pos), edges) in edges_otherside {
            let mut edge = DitchEdge::new(!first_pos, first_id, to_pos, to);
            edge.edges.extend(edges);
            new_node.edges.push(edge);
        }
        self.nodes[first_id] = new_node;
        {
            let mut edge_to_removednodes = vec![];
            self.nodes[root]
                .edges
                .iter_mut()
                .filter(|e| e.from_position == position)
                .for_each(|edge| {
                    edge_to_removednodes.append(&mut edge.edges);
                });
            self.nodes[root].edges.retain(|e| !e.edges.is_empty());
            let mut new_edge = DitchEdge::new(position, root, first_pos, first_id);
            new_edge.edges = edge_to_removednodes;
            self.nodes[root].edges.push(new_edge);
        }
        for (other, other_position) in node_otherside {
            let mut edge_to_removednodes = vec![];
            self.nodes[other]
                .edges
                .iter_mut()
                .filter(|e| merged_nodes_ids.contains(&e.to))
                .for_each(|edge| {
                    edge_to_removednodes.append(&mut edge.edges);
                });
            self.nodes[other].edges.retain(|e| !e.edges.is_empty());
            let mut new_edge = DitchEdge::new(other_position, other, !first_pos, first_id);
            new_edge.edges = edge_to_removednodes;
            self.nodes[other].edges.push(new_edge);
        }
        let mut merged_nodes_ids = merged_nodes_ids;
        merged_nodes_ids.retain(|&x| x != first_id);
        ((first_id, !first_pos), merged_nodes_ids)
    }
    // Traverse from the given `start` node.
    fn traverse_from(
        &self,
        arrived: &mut [bool],
        sids: &mut [ContigTag],
        start: usize,
        start_position: Position,
        seqname: String,
        c: &AssembleConfig,
    ) -> (gfa::Segment, Vec<gfa::Edge>, ContigSummary) {
        // Find edges.
        let mut edges: Vec<_> = self.nodes[start]
            .edges
            .iter()
            .filter(|e| e.from_position == start_position)
            .filter_map(|e| {
                assert!(e.from == start);
                let (sid2, beg2) = match sids[e.to] {
                    ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(0, false);
                let a = None;
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            })
            .collect();
        let (mut node, mut position) = (start, start_position);
        let mut seq = self.initial_sequence(start, start_position);
        // Start traveresing.
        let mut unit_names = vec![];
        loop {
            arrived[node] = true;
            // Move forward.
            let cons = match position {
                Position::Head if c.to_polish => self.nodes[node].consensus(),
                Position::Head => self.nodes[node].seq(),
                Position::Tail if c.to_polish => revcmp_str(&self.nodes[node].consensus()),
                Position::Tail => revcmp_str(&self.nodes[node].seq()),
            };
            {
                let direction = match position {
                    Position::Head => true,
                    Position::Tail => false,
                };
                let elm = ContigElement {
                    unit: self.nodes[node].unit,
                    cluster: self.nodes[node].cluster,
                    strand: direction,
                };
                unit_names.push(elm);
            }
            seq += &cons;
            position = !position;
            // Check.
            let num_edges = self.nodes[node]
                .edges
                .iter()
                .filter(|e| e.from_position == position)
                .count();
            if num_edges == 0 || num_edges > 1 {
                break;
            }
            // There is only one child.
            let selected_edge = self.nodes[node]
                .edges
                .iter()
                .find(|e| e.from_position == position)
                .unwrap();
            assert_eq!(selected_edge.from, node);
            assert_eq!(selected_edge.from_position, position);
            // Succeed along with the edge.
            let (offset, label) = {
                let offset = selected_edge
                    .edges
                    .iter()
                    .map(|(e, _)| e.offset)
                    .sum::<i64>();
                let offset = offset / selected_edge.edges.len() as i64;
                let label = if offset <= 0 {
                    String::new()
                } else {
                    selected_edge
                        .edges
                        .iter()
                        .max_by_key(|x| x.0.label().len())
                        .and_then(|&(ref e, is_forward)| {
                            if is_forward {
                                String::from_utf8(e.label().to_vec()).ok()
                            } else {
                                String::from_utf8(bio_utils::revcmp(e.label())).ok()
                            }
                        })
                        .unwrap_or_else(String::new)
                    // let xs: Vec<_> = selected_edge
                    //     .edges
                    //     .iter()
                    //     .map(|&(ref e, is_forward)| {
                    //         if is_forward {
                    //             e.label().to_vec()
                    //         } else {
                    //             bio_utils::revcmp(e.label())
                    //         }
                    //     })
                    //     .collect();
                    // let xs: Vec<_> = xs
                    //     .iter()
                    //     .map(|x| x.as_slice())
                    //     .filter(|x| !x.is_empty())
                    //     .collect();
                    // consensus(&xs, 10).into_iter().map(|x| x as char).collect()
                };
                (offset, label)
            };
            if offset >= 0 {
                seq += &label;
            } else {
                (0..(-offset) as usize).all(|_| seq.pop().is_some());
            }
            let (next, next_position) = (selected_edge.to, selected_edge.to_position);
            // Check the number of child.
            let num_children = self.nodes[next]
                .edges
                .iter()
                .filter(|e| e.from_position == next_position)
                .count();
            if num_children >= 2 || arrived[next] {
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
        if start == node {
            sids[node] = ContigTag::Both(seqname.clone(), start_position, position, seq.len());
        } else {
            sids[start] = ContigTag::Start(seqname.clone(), start_position, seq.len());
            sids[node] = ContigTag::End(seqname.clone(), position, seq.len());
        }
        let seg = gfa::Segment::from(seqname.clone(), seq.len(), Some(seq.clone()));
        // Add gfa edges.
        let tail_edges = self.nodes[node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .filter_map(|e| {
                assert!(e.from == node);
                let (sid2, beg2) = match sids[e.to] {
                    ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::End(ref name, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(0, false))
                    }
                    ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
                        (gfa::RefID::from(name, true), gfa::Position::from(len, true))
                    }
                    _ => return None,
                };
                let eid = None;
                let sid1 = gfa::RefID::from(&seqname, true);
                let beg1 = gfa::Position::from(seq.len(), true);
                let a = None;
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            });
        edges.extend(tail_edges);
        (seg, edges, summary)
    }
    fn initial_sequence(&self, node: usize, position: Position) -> String {
        let num_edges = self.nodes[node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count();
        if num_edges > 0 {
            String::new()
        } else {
            // let mut seq: Vec<_> = self.nodes[node]
            //     .tips
            //     .iter()
            //     .filter(|tip| tip.position == position && !tip.seq.is_empty())
            //     .map(|tip| match tip.in_direction {
            //         true => tip.seq.to_vec(),
            //         false => bio_utils::revcmp(tip.seq),
            //     })
            //     .collect();
            // seq.sort_by_key(|x| x.len());
            // seq.reverse();
            // let result = super::naive_consensus::consensus(&seq);
            self.nodes[node]
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
            //String::from_utf8(result).unwrap()
        }
    }
    fn trailing_sequence(&self, node: usize, position: Position) -> String {
        let num_edges = self.nodes[node]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .count();
        if num_edges > 0 {
            String::new()
        } else {
            // let mut seq: Vec<_> = self.nodes[node]
            //     .tips
            //     .iter()
            //     .filter(|tip| tip.position == position && !tip.seq.is_empty())
            //     .map(|tip| match tip.in_direction {
            //         true => bio_utils::revcmp(tip.seq),
            //         false => tip.seq.to_vec(),
            //     })
            //     .collect();
            // seq.sort_by_key(|x| x.len());
            // seq.reverse();
            // String::from_utf8(super::naive_consensus::consensus(&seq)).unwrap()
            self.nodes[node]
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

pub fn consensus(xs: &[&[u8]], num: usize) -> Vec<u8> {
    if xs.is_empty() {
        return vec![];
    } else if xs.iter().map(|xs| xs.len()).max().unwrap() < 10 {
        return xs.iter().max_by_key(|x| x.len()).unwrap().to_vec();
    }
    let param = (-2, -2, &|x, y| if x == y { 2 } else { -4 });
    use rand_xoshiro::Xoroshiro128StarStar;
    let cs: Vec<_> = (0..num as u64)
        .into_par_iter()
        .map(|s| {
            use rand::seq::SliceRandom;
            let mut rng: Xoroshiro128StarStar = rand::SeedableRng::seed_from_u64(s);
            let mut cs: Vec<_> = xs.to_vec();
            cs.shuffle(&mut rng);
            let max_len = cs.iter().map(|s| s.len()).max().unwrap_or(0);
            let node_num_thr = (max_len as f64 * 1.5).floor() as usize;
            cs.iter()
                .fold(poa_hmm::POA::default(), |x, y| {
                    let res = if x.nodes().len() > node_num_thr {
                        x.add(y, 1., param).remove_node(0.4)
                    } else {
                        x.add(y, 1., param)
                    };
                    res
                })
                .remove_node(0.4)
                .finalize()
                .consensus()
        })
        .collect();
    let cs: Vec<_> = cs.iter().map(|cs| cs.as_slice()).collect();
    poa_hmm::POA::default()
        .update_thr(&cs, &vec![1.; cs.len()], param, 0.8, 1.5)
        .consensus()
}

#[cfg(test)]
mod tests {
    #![allow(dead_code)]
    use definitions::Node;
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
        let _graph = DitchGraph::new(&reads, &c);
    }
    #[test]
    fn generate() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
        let _ = graph.spell(&c, 0);
    }
    #[test]
    fn validate() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        let (segment, _edge, _group) = graph.spell(&c, 0);
        eprintln!("{:?}", segment);
        eprintln!("{:?}", answer);
        assert_eq!(segment.len(), 1);
        assert_eq!(segment[0].sequence.as_ref().unwrap(), answer[0].as_str());
    }
    #[test]
    fn validate_large() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_large();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let mut graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        graph.collapse_buddle(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let mut graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        graph.collapse_buddle(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let mut graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        graph.collapse_buddle(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let mut graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        graph.collapse_buddle(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
        let mut graph = DitchGraph::new(&reads, &c);
        eprintln!("{:?}", graph);
        graph.collapse_buddle(&c);
        eprintln!("{:?}", graph);
        let (segments, edges, group) = graph.spell(&c, 0);
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
                .map(gfa::Content::Edge)
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
}
