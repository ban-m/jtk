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
        .map(|n| gfa::Content::Seg(n))
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(nodes);
    let edges = edges
        .into_iter()
        .map(|n| gfa::Content::Edge(n))
        .map(|n| gfa::Record::from_contents(n, vec![]));
    records.extend(edges);
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    records.push(group);
    records
}

#[derive(Debug, Clone)]
pub struct DitchGraph<'a, 'b> {
    nodes: Vec<DitchNode<'a, 'b>>,
    index: HashMap<NodeIndex, usize>,
}

impl<'a, 'b> std::fmt::Display for DitchGraph<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edge = self.nodes.iter().map(|e| e.edges.len()).sum::<usize>();
        writeln!(f, "Nodes:{}\tEdges:{}", self.nodes.len(), edge)?;
        for node in self.nodes.iter() {
            let idx = self
                .index
                .get(&NodeIndex::new(node.unit, node.cluster))
                .unwrap();
            write!(f, "{}\t{}", idx, node)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct DitchNode<'a, 'b> {
    unit: u64,
    cluster: u64,
    nodes: Vec<&'a definitions::Node>,
    edges: Vec<DitchEdge<'b>>,
}

impl<'a, 'b> std::fmt::Display for DitchNode<'a, 'b> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{}:{}", self.unit, self.cluster)?;
        writeln!(f, "Seq:")?;
        for (idx, n) in self.nodes.iter().enumerate() {
            writeln!(f, "{:3}:{}", idx, n.seq)?;
        }
        for (idx, e) in self.edges.iter().enumerate() {
            writeln!(f, "{:3}:->{}", idx, e)?;
        }
        Ok(())
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
    edges: Vec<&'a definitions::Edge>,
}

impl<'a> std::fmt::Display for DitchEdge<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}:{})->", self.from, self.from_position)?;
        write!(f, "({}:{}):", self.to, self.to_position)?;
        let edgs: Vec<_> = self
            .edges
            .iter()
            .map(|x| format!("{}", x.label.len()))
            .collect();
        write!(f, "{}", edgs.join(","))
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
    fn new(from: &Node, from_i: usize, to: &Node, to_i: usize) -> Self {
        let from_position = if from.is_forward {
            Position::Tail
        } else {
            Position::Head
        };
        let to_position = if to.is_forward {
            Position::Head
        } else {
            Position::Tail
        };
        Self {
            from: from_i,
            to: to_i,
            from_position,
            to_position,
            edges: vec![],
        }
    }
    fn push(&mut self, x: &'a Edge) {
        self.edges.push(x);
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
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

impl<'a> DitchGraph<'a, 'a> {
    fn new(reads: &[&'a EncodedRead], _c: &AssembleConfig) -> Self {
        // Allocate all nodes.
        let (index, mut nodes) = {
            let mut index = HashMap::new();
            let mut nodes = vec![];
            for node in reads.iter().flat_map(|r| r.nodes.iter()) {
                let node_index = NodeIndex::new(node.unit, node.cluster);
                if !index.contains_key(&node_index) {
                    index.insert(node_index, nodes.len());
                    nodes.push(DitchNode::new(node.unit, node.cluster));
                }
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
                let mut dedg = DitchEdge::new(from, from_index, to, to_index);
                match nodes[from_index].edges.iter_mut().find(|x| x == &&dedg) {
                    Some(res) => res.push(edge),
                    None => {
                        dedg.push(edge);
                        nodes[from_index].edges.push(dedg);
                    }
                }
            }
        }
        Self { index, nodes }
    }
    fn spell(
        &self,
        c: &AssembleConfig,
        cl: usize,
    ) -> (Vec<gfa::Segment>, Vec<gfa::Edge>, gfa::Group) {
        let mut starts = self.enumerate_candidates(c);
        // Arrived. Head/Tail.
        let mut tail_arrived = vec![false; self.nodes.len()];
        let mut head_arrived = vec![false; self.nodes.len()];
        let mut sids = vec![None; self.nodes.len()];
        let (mut g_segs, mut g_edges) = (vec![], vec![]);
        'outer: while let Some(start) = starts.pop() {
            if tail_arrived[start] || head_arrived[start] {
                continue;
            }
            // Start traversing.
            
        }
        let uid = Some(format!("group-{}", cl));
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.eid.clone()))
            .collect();
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group)
    }
    fn enumerate_candidates(&self, c: &AssembleConfig) -> Vec<(usize, Position)> {
        let head_degree = vec![0; self.nodes.len()];
        let tail_degree = vec![0; self.nodes.len()];
        for node in self.nodes.iter() {
            for edge in node.edges.iter() {
                match edge.to_position {
                    Position::Head => head_degree[edge.to] += 1,
                    Position::Tail => tail_degree[edge.to] += 1,
                }
            }
        }
        let head_degree = head_degree
            .iter()
            .enumerate()
            .filter(|&(_, deg)| deg > 1)
            .map(|&(idx, _)| (idx, Position::Head));
        let tail_degree = tail_degree
            .iter()
            .enumerate()
            .filter(|&(_, deg)| deg > 1)
            .map(|&(idx, _)| (idx, Position::Tail));
        let mut candidates: Vec<_> = head_degree.chain(tail_degree).collect();
        debug!("{} starts candidates.", starts.len());
        // These are "phony" targets. This is added to treat simple cycle.
        // In other words, if there is a simple cycle, all nodes in that cycle
        // would have degree of 1, and thus never be a seed of a path.
        // Therefore, I put all the element into the candidates.
        candidates.extend((0..self.nodes.len()).map(|i| (i, Position::Head)));
        candidates.extend((0..self.nodes.len()).map(|i| (i, Position::Tail)));
        candidates.reverse();
        candidates
    }
}

#[cfg(test)]
mod tests {
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
    fn gen_edge(offset: i64, label: String) -> Edge {
        let mut edge = Edge::default();
        edge.offset = offset;
        edge.label = label;
        edge
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
        read.edges.push(gen_edge(1, seq[10..11].iter().collect()));
        read.edges.push(gen_edge(-6, String::new()));
        read.edges.push(gen_edge(2, seq[27..29].iter().collect()));
        let reads = vec![read];
        let seg = vec![gfa::Segment {
            sid: "test".to_string(),
            slen: 30,
            sequence: Some(seq.iter().collect()),
        }];
        (reads, vec![seq.iter().collect()])
    }
    fn gen_mock_large() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        gen_mock1()
    }
    fn gen_mock_complex() -> DataPack {
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(10);
        gen_mock1()
    }
    use super::*;
    #[test]
    fn create() {
        let c = AssembleConfig::default();
        let (reads, _) = gen_mock1();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
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
        let (segment, edge, group) = graph.spell(&c, 0);
        assert_eq!(segment.len(), 1);
        assert_eq!(segment[0].sequence.as_ref().unwrap(), answer[0].as_str());
    }
    #[test]
    fn validate_large() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_large();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
        let (segment, edge, group) = graph.spell(&c, 0);
        // Assertion.
    }
    #[test]
    fn validate_complex() {
        let c = AssembleConfig::default();
        let (reads, answer) = gen_mock_complex();
        let reads: Vec<_> = reads.iter().collect();
        let graph = DitchGraph::new(&reads, &c);
        let (segment, edge, group) = graph.spell(&c, 0);
        // Assertion.
    }
}
