use definitions::*;
use gfa::GFA;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct AssembleConfig {}
impl std::default::Default for AssembleConfig {
    fn default() -> Self {
        Self {}
    }
}

pub trait Assemble {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA;
}

impl Assemble for DataSet {
    fn assemble_as_gfa(&self, c: &AssembleConfig) -> GFA {
        let mut clusters: HashMap<_, Vec<_>> = HashMap::new();
        {
            // ID=>Cluster
            let id2cluster: HashMap<u64, usize> = self
                .assignments
                .iter()
                .map(|asn| (asn.id, asn.cluster))
                .collect();
            for read in self.encoded_reads.iter() {
                if let Some(&cluster) = id2cluster.get(&read.id) {
                    clusters.entry(cluster).or_default().push(read);
                }
            }
        }
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let records: Vec<_> = clusters
            .into_iter()
            .flat_map(|(cl, reads)| cluster_to_gfa(cl, reads, c))
            .chain(std::iter::once(header))
            .collect();
        GFA::from_records(records)
    }
}

fn cluster_to_gfa(cl: usize, reads: Vec<&EncodedRead>, c: &AssembleConfig) -> Vec<gfa::Record> {
    let graph = DitchGraph::new(&reads, c);
    let mut records = vec![];
    let (nodes, edges, group) = graph.spell(c, cl);
    let nodes = nodes
        .into_iter()
        .map(gfa::Content::Seg)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(nodes);
    let edges = edges
        .into_iter()
        .map(gfa::Content::Edge)
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(edges);
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    records.push(group);
    records
}

/// Ditch Graph
/// Each unit of the dataset consists of a pair of node, named
/// tail-node and head-node, and it is represented as a 'node'in a
/// Ditch graph.
/// Each node has several edges, induced by the connection inside reads.
#[derive(Debug, Clone)]
pub struct DitchGraph<'a, 'b> {
    nodes: Vec<DitchNode<'a, 'b>>,
    index: HashMap<NodeIndex, usize>,
}

impl<'a, 'b> std::fmt::Display for DitchGraph<'a, 'b> {
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
    fn new(reads: &[&'a EncodedRead], _c: &AssembleConfig) -> Self {
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
        // Append edge.
        for read in reads.iter() {
            if read.nodes.len() != read.edges.len() {
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
                let h_edge_num = self.nodes[i]
                    .edges
                    .iter()
                    .filter(|e| e.from_position == Position::Head)
                    .count();
                let t_edge_num = self.nodes[i]
                    .edges
                    .iter()
                    .filter(|e| e.from_position == Position::Tail)
                    .count();
                if h_edge_num != 1 {
                    selected[i] = true;
                    Some((i, Position::Head))
                } else if t_edge_num != 1 {
                    selected[i] = true;
                    Some((i, Position::Tail))
                } else {
                    None
                }
            })
            .collect();
        let secondary_candidates =
            (0..self.nodes.len())
                .filter(|&i| !selected[i])
                .filter_map(|i| {
                    // There is only one edge from (i, Position::Head).
                    // Thus, the unwrap() never panics.
                    let head_grand_child = self.nodes[i]
                        .edges
                        .iter()
                        .find(|e| e.from_position == Position::Head)
                        .map(|e| {
                            let (next, next_position) = (e.to, e.to_position);
                            self.nodes[next]
                                .edges
                                .iter()
                                .filter(|e| e.from_position == next_position)
                                .count()
                        })
                        .unwrap();
                    assert!(head_grand_child != 0);
                    if head_grand_child > 1 {
                        return Some((i, Position::Head));
                    }
                    let tail_grand_child = self.nodes[i]
                        .edges
                        .iter()
                        .find(|e| e.from_position == Position::Tail)
                        .map(|e| {
                            let (next, next_position) = (e.to, e.to_position);
                            self.nodes[next]
                                .edges
                                .iter()
                                .filter(|e| e.from_position == next_position)
                                .count()
                        })
                        .unwrap();
                    assert!(tail_grand_child != 0);
                    if tail_grand_child > 1 {
                        return Some((i, Position::Tail));
                    }
                    None
                });
        primary_candidates.extend(secondary_candidates);
        primary_candidates
    }
    // Assemble the ditch graph.
    // In other words, it reduce the simple path, currently.
    fn spell(
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
            eprintln!("First loop:{}", i);
            let name = format!("tig_{:03}_{:03}", cl, g_segs.len());
            let (contig, edges) = self.traverse_from(&mut arrived, &mut sids, i, p, name, c);
            g_segs.push(contig);
            g_edges.extend(edges);
        }
        for i in 0..self.nodes.len() {
            if arrived[i] {
                continue;
            }
            eprintln!("Second loop:{}", i);
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
        loop {
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
        eprintln!("{}", graph);
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
        eprintln!("{}", graph);
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
        eprintln!("{}", graph);
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
        eprintln!("{}", graph);
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
}
