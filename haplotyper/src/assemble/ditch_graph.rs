use super::AssembleConfig;
use definitions::{Edge, EncodedRead};
use std::collections::HashMap;
/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
#[derive(Clone)]
pub struct DitchGraph<'a, 'b> {
    nodes: Vec<DitchNode<'a, 'b>>,
    index: HashMap<NodeIndex, usize>,
}

impl<'a, 'b> std::fmt::Display for DitchGraph<'a, 'b> {
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

impl<'a, 'b> std::fmt::Debug for DitchGraph<'a, 'b> {
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
pub struct DitchNode<'a, 'b> {
    unit: u64,
    cluster: u64,
    // CAUTION!!!!!! The `unit` and `cluster` members should not be look-upped!
    nodes: Vec<&'a definitions::Node>,
    edges: Vec<DitchEdge<'b>>,
}

impl<'a, 'b> std::fmt::Display for DitchNode<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "---------{}:{}---------", self.unit, self.cluster)?;
        writeln!(f, "Seq:")?;
        for (idx, n) in self.nodes.iter().enumerate() {
            writeln!(f, "{:3}:{}", idx, n.seq)?;
        }
        let lines: Vec<_> = self.edges.iter().map(|e| format!("{}", e)).collect();
        write!(f, "Edges:\n{}", lines.join("\n"))
    }
}

impl<'a, 'b> DitchNode<'a, 'b> {
    fn new(unit: u64, cluster: u64) -> Self {
        Self {
            unit,
            cluster,
            nodes: vec![],
            edges: vec![],
        }
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

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
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

fn revcmp(seq: &str) -> String {
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

impl<'a> DitchGraph<'a, 'a> {
    pub fn new(reads: &[&'a EncodedRead], _c: &AssembleConfig) -> Self {
        // Allocate all nodes.
        let (index, mut nodes) = {
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
            (index, nodes)
        };
        // Append edges.
        for read in reads.iter() {
            if read.nodes.len() != read.edges.len() + 1 {
                debug!("{}\t{}\t{}", read.id, read.nodes.len(), read.edges.len());
            }
            for (pairs, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
                let (from, to) = (&pairs[0], &pairs[1]);
                let from_index = *index.get(&NodeIndex::new(from.unit, from.cluster)).unwrap();
                let to_index = *index.get(&NodeIndex::new(to.unit, to.cluster)).unwrap();
                let from_pos = if from.is_forward {
                    Position::Tail
                } else {
                    Position::Head
                };
                let to_pos = if to.is_forward {
                    Position::Head
                } else {
                    Position::Tail
                };
                let mut dedg = DitchEdge::new(from_pos, from_index, to_pos, to_index);
                match nodes[from_index].edges.iter_mut().find(|x| x == &&dedg) {
                    Some(res) => res.push((edge, true)),
                    None => {
                        dedg.push((edge, true));
                        nodes[from_index].edges.push(dedg);
                    }
                }
                let mut dedg = DitchEdge::new(to_pos, to_index, from_pos, from_index);
                match nodes[to_index].edges.iter_mut().find(|x| x == &&dedg) {
                    Some(res) => res.push((edge, false)),
                    None => {
                        dedg.push((edge, false));
                        nodes[to_index].edges.push(dedg);
                    }
                }
            }
        }
        Self { index, nodes }
    }
    fn enumerate_candidates(&self) -> Vec<(usize, Position)> {
        let mut selected = vec![false; self.nodes.len()];
        let mut primary_candidates: Vec<_> = (0..self.nodes.len())
            .filter_map(|i| {
                for position in vec![Position::Head, Position::Tail] {
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
                    for position in vec![Position::Head, Position::Tail] {
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
    ) -> (Vec<gfa::Segment>, Vec<gfa::Edge>, gfa::Group) {
        let mut arrived = vec![false; self.nodes.len()];
        let mut sids = vec![ContigTag::None; self.nodes.len()];
        let (mut g_segs, mut g_edges) = (vec![], vec![]);
        let candidates = self.enumerate_candidates();
        for (i, p) in candidates {
            if arrived[i] {
                continue;
            }
            debug!("First loop:{}", i);
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges) = self.traverse_from(&mut arrived, &mut sids, i, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
        }
        for i in 0..self.nodes.len() {
            if arrived[i] {
                continue;
            }
            debug!("Second loop:{}", i);
            let p = Position::Head;
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges) = self.traverse_from(&mut arrived, &mut sids, i, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
        }
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.eid.clone()))
            .collect();
        let uid = Some(format!("group-{}", cl));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group)
    }
    // Removing nodes and re-map all of the indices.
    fn remove_nodes(&mut self, to_remove: &[bool]) {
        let mapping = {
            let mut mapping = vec![];
            let mut index = 0;
            for &b in to_remove {
                mapping.push(index);
                index += !b as usize;
            }
            mapping
        };
        self.index.iter_mut().for_each(|(_, v)| *v = mapping[*v]);
        {
            let mut idx = 0;
            self.nodes.retain(|_| {
                idx += 1;
                !to_remove[idx - 1]
            });
        }
        self.nodes.iter_mut().for_each(|node| {
            node.edges.iter_mut().for_each(|e| {
                e.from = mapping[e.from];
                e.to = mapping[e.to];
            })
        });
    }
    pub fn collapse_buddle(&mut self, c: &AssembleConfig) {
        let mut to_remove = vec![false; self.nodes.len()];
        for i in 0..self.nodes.len() {
            for position in vec![Position::Head, Position::Tail] {
                let edges: Vec<_> = self.nodes[i]
                    .edges
                    .iter()
                    .filter(|e| e.from_position == position && e.from == i)
                    .collect();
                if edges.len() > 1 {
                    // Never panic.
                    let pos = edges[0].to_position;
                    let unit = self.nodes[edges[0].to].unit;
                    if edges
                        .iter()
                        .all(|n| n.to_position == pos && unit == self.nodes[n.to].unit)
                    {
                        // The head side of the i-th node has two or mode cluster,
                        // which has the same distination.
                        // Check the collaption criteria, and collapse if possible.
                        eprintln!(
                            "Removing a bubble starting from {}-{}-{:?}",
                            self.nodes[i].unit, self.nodes[i].cluster, position
                        );
                        for i in self.collapse_bubble_from(i, position, c) {
                            // eprintln!("Removing {}", i);
                            to_remove[i] = true;
                        }
                    }
                }
            }
        }
        self.remove_nodes(&to_remove);
    }
    fn collapse_bubble_from(
        &mut self,
        i: usize,
        position: Position,
        _c: &AssembleConfig,
    ) -> Vec<usize> {
        // Check the collapsing condition.
        let edges: Vec<_> = self.nodes[i]
            .edges
            .iter()
            .filter(|e| e.from_position == position)
            .inspect(|e| assert!(e.from == i))
            .collect();
        // Check.
        assert!(edges.len() > 1);
        let (first_id, first_pos, first_unit) = {
            let first = edges.first().unwrap();
            let first_unit = self.nodes[first.to].unit;
            (first.to, first.to_position, first_unit)
        };
        let merged_nodes_ids: Vec<_> = edges.iter().map(|e| e.to).collect();
        assert!(edges
            .iter()
            .all(|e| e.to_position == first_pos && first_unit == self.nodes[e.to].unit));
        // Get the candidate child.
        let pos = !first_pos;
        let (child_id, child_pos) = match self.nodes[first_id]
            .edges
            .iter()
            .find(|e| e.from_position == pos && e.from == first_id)
        {
            Some(child) => (child.to, child.to_position),
            None => return vec![],
        };
        // Check all distinations of `i` has only one parent, `i`.
        // Check all otherside of the distinations of `i` has only one child.
        let is_collapsable_bubble = edges.iter().all(|e| {
            let res = self.nodes[e.to].edges.iter().all(|f| {
                (f.to_position == position && f.to == i)
                    || (f.to_position == child_pos && f.to == child_id)
            });
            res && self.nodes[e.to].edges.len() == 2
        });
        if !is_collapsable_bubble {
            vec![]
        } else {
            // Add to the `first_id` from other edges.
            let remove_nodes: Vec<_> = edges.into_iter().skip(1).map(|e| e.to).collect();
            let mut node_result = vec![];
            let mut edge_result: Vec<(&definitions::Edge, bool)> = vec![];
            for &remove_node in remove_nodes.iter() {
                node_result.append(&mut self.nodes[remove_node].nodes);
                self.nodes[remove_node]
                    .edges
                    .iter_mut()
                    .for_each(|ditch_edge| {
                        edge_result.append(&mut ditch_edge.edges);
                    });
            }
            self.nodes[first_id].nodes.extend(node_result);
            self.nodes[first_id]
                .edges
                .iter_mut()
                .find(|ditch_edge| ditch_edge.to == i && ditch_edge.to_position == position)
                .unwrap()
                .edges
                .extend(edge_result);
            // Change edges from nodes[i]
            let edge_result = self.nodes[i]
                .edges
                .iter_mut()
                .filter(|e| e.from_position == position)
                .skip(1)
                .fold(vec![], |mut x, ditch_edge| {
                    x.append(&mut ditch_edge.edges);
                    x
                });
            self.nodes[i]
                .edges
                .iter_mut()
                .find(|e| e.from_position == position)
                .unwrap()
                .edges
                .extend(edge_result);
            self.nodes[i]
                .edges
                .retain(|ditch_edge| !ditch_edge.edges.is_empty());
            let edge_result = self.nodes[child_id]
                .edges
                .iter_mut()
                .filter(|e| merged_nodes_ids.contains(&e.to))
                .skip(1)
                .fold(vec![], |mut x, d_edge| {
                    x.append(&mut d_edge.edges);
                    x
                });
            // Never panic.
            self.nodes[child_id]
                .edges
                .iter_mut()
                .find(|e| merged_nodes_ids.contains(&e.to))
                .unwrap()
                .edges
                .extend(edge_result);
            self.nodes[child_id]
                .edges
                .retain(|ditch_edge| !ditch_edge.edges.is_empty());
            remove_nodes
        }
    }
    // Traverse from the given `start` node.
    fn traverse_from(
        &self,
        arrived: &mut [bool],
        sids: &mut [ContigTag],
        start: usize,
        start_position: Position,
        seqname: String,
        _c: &AssembleConfig,
    ) -> (gfa::Segment, Vec<gfa::Edge>) {
        // Find edges.
        let mut edges: Vec<_> = self.nodes[start]
            .edges
            .iter()
            .filter(|e| e.from_position == start_position)
            .filter_map(|e| {
                assert!(e.from == start);
                let (sid2, beg2) = match &sids[e.to] {
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
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            })
            .collect();
        let (mut node, mut position) = (start, start_position);
        let mut seq = String::new();
        // Start traveresing.
        let mut unit_names = vec![];
        loop {
            unit_names.push((self.nodes[node].unit, self.nodes[node].cluster));
            arrived[node] = true;
            // Move forward.
            match position {
                Position::Head => seq += &self.nodes[node].nodes[0].seq,
                Position::Tail => {
                    let s = &self.nodes[node].nodes[0].seq;
                    seq += &revcmp(s);
                }
            };
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
            // if let Some(&(ref edge, is_forward)) =
            let &(ref edge, is_forward) = selected_edge.edges.first().unwrap();
            if edge.offset >= 0 {
                if is_forward {
                    seq += &edge.label;
                } else {
                    seq += &revcmp(&edge.label);
                }
            } else {
                assert!(-edge.offset > 0);
                (0..(-edge.offset) as usize).all(|_| seq.pop().is_some());
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
        let unit_names: Vec<_> = unit_names
            .iter()
            .map(|(u, c)| format!("{}:{}", u, c))
            .collect();
        debug!("{}\t{}", seqname, unit_names.join("\t"));
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
                let (sid2, beg2) = match &sids[e.to] {
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
                Some(gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a))
            });
        edges.extend(tail_edges);
        (seg, edges)
    }
}

#[cfg(test)]
mod tests {
    #![allow(dead_code)]
    use rand::SeedableRng;
    use definitions::Node;
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
    fn gen_remove_test_2() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(129180);
        let units: Vec<String> = (0..4)
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
        let (mut read4, mut read5) = (read.clone(), read.clone());
        for unit in 1..4 {
            let cl = if unit == 1 { 1 } else { 2 };
            let seq = units[unit].clone();
            let unit = unit as u64;
            read4.nodes.push(gen_node(unit, cl, true, seq));
            if unit < 3 {
                read5.edges.push(gen_edge(unit, unit + 1, 0, edge.clone()));
            }
        }
        for unit in (1..4).rev() {
            let cl = if unit == 1 { 1 } else { 2 };
            let seq: Vec<_> = units[unit].chars().collect();
            let seq: String = revcmp(&seq).iter().collect();
            let unit = unit as u64;
            read5.nodes.push(gen_node(unit, cl, false, seq));
            if 0 < unit {
                read5.edges.push(gen_edge(unit, unit - 1, 0, edge.clone()));
            }
        }
        let reads = vec![read0, read1, read2, read3, read4, read5];
        let mut answers = units[0..3].to_vec();
        answers.push(units[2].clone() + &units[3]);
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
        let units: Vec<String> = (0..3)
            .map(|_| (0..100).filter_map(|_| BASES.choose(&mut rng)).collect())
            .collect();
        let edge = String::new();
        let read = EncodedRead::default();
        let (mut read0, mut read1, mut read2) = (read.clone(), read.clone(), read.clone());
        for unit in 0..3 {
            let seq = units[unit].clone();
            let unit = unit as u64;
            let (cl0, cl1) = if unit == 1 { (0, 1) } else { (0, 0) };
            let cl2 = if unit == 0 { 0 } else { 2 };
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
            let (cl3, cl4) = if unit == 1 { (0, 1) } else { (0, 0) };
            let cl5 = if unit == 0 { 0 } else { 2 };
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
        let mut answers = vec![];
        answers.push(units[0].clone());
        answers.push(units[1].clone() + &units[2]);
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
        assert_eq!(edges.len(), 5);
        assert_eq!(segments.len(), 5);
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
