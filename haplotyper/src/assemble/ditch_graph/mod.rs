use super::copy_number::CoverageCalibrator;
use super::AssembleConfig;
use definitions::ReadType;
use definitions::{EncodedRead, Unit};
use rayon::prelude::*;
mod sequence_generation;
use std::collections::HashMap;
use std::collections::HashSet;
pub mod dg_test;
const ERROR_PROB: f64 = 0.05;
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
type EdgeBetweenSimplePath = (usize, bool, usize, bool, f64);

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

impl std::fmt::Display for Focus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let &Self {
            from,
            from_position,
            to,
            to_position,
            dist,
            log_likelihood_ratio,
        } = self;
        let (fu, fc) = from;
        let (tu, tc) = to;
        write!(
            f,
            "{},{},{},{},{},{},{},{:.2}",
            fu, fc, from_position, tu, tc, to_position, dist, log_likelihood_ratio
        )
    }
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
    // If this is not empty, then this label contains these node, from the front to the end,
    // as indicated direction(true=>forward, head->tail, false=>reverse, tail->head).
    // It also contains the occurence of the node, estimated.
    proxying: Vec<(Node, bool, usize)>,
}

impl std::fmt::Display for DitchEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({:?}:{})->", self.from, self.from_position)?;
        write!(f, "({:?}:{})", self.to, self.to_position)?;
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
            proxying: Vec::new(),
        }
    }
    fn with_proxy(
        from: Node,
        from_position: Position,
        to: Node,
        to_position: Position,
        seq: EdgeLabel,
        proxying: Vec<(Node, bool, usize)>,
    ) -> Self {
        Self {
            from,
            to,
            from_position,
            to_position,
            occ: 0,
            seq,
            copy_number: None,
            proxying,
        }
    }

    // Reverse this edge.
    fn reverse(&self) -> Self {
        let proxying: Vec<_> = self
            .proxying
            .iter()
            .rev()
            .map(|&(n, d, o)| (n, !d, o))
            .collect();
        Self {
            from: self.to,
            from_position: self.to_position,
            to: self.from,
            to_position: self.from_position,
            seq: self.seq.reverse(),
            occ: 0,
            copy_number: self.copy_number,
            proxying,
        }
    }
    // reduce the copy number by one,
    // and decrease the .occ parameter appropriately.
    fn reduce_one_copy_number(&mut self) {
        match self.copy_number.as_mut() {
            Some(cp) if *cp == 1 => {
                let occ_reduce = self.occ / *cp;
                self.occ -= occ_reduce;
                *cp -= 1;
            }
            _ => {}
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
        write!(f, "{}\t{}", self.id, line.join("\t"))
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

fn take_consensus<R: std::borrow::Borrow<EncodedRead>>(
    reads: &[R],
    units: Option<&[Unit]>,
    read_type: &definitions::ReadType,
    c: &AssembleConfig,
) -> (NodeConsensus, EdgeConsensus) {
    let mut nodes: HashMap<_, Vec<_>> = HashMap::new();
    let mut edges: HashMap<_, Vec<_>> = HashMap::new();
    for read in reads.iter().map(|r| r.borrow()) {
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
            let (tuple, value) = if (from, from_pos) <= (to, to_pos) {
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
    // use kiley::gphmm::*;
    // let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    // let config = kiley::PolishConfig::with_model(100, 0, 0, 0, 0, hmm);
    let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let node_consensus: HashMap<_, _> = nodes
        .into_par_iter()
        .map(|(key, value)| {
            let draft = match units.get(&(key.0)) {
                Some(res) => res.to_vec(),
                None => value.get(0).map(|x| x.seq().to_vec()).unwrap(),
            };
            let cons = if c.to_polish {
                let seqs: Vec<_> = value.iter().map(|x| x.seq()).collect();
                use definitions::Op;
                let mut ops: Vec<Vec<_>> = value
                    .iter()
                    .map(|x| {
                        x.cigar
                            .iter()
                            .flat_map(|op| {
                                let (op, len) = match op {
                                    Op::Match(l) => (kiley::Op::Match, *l),
                                    Op::Del(l) => (kiley::Op::Del, *l),
                                    Op::Ins(l) => (kiley::Op::Ins, *l),
                                };
                                std::iter::repeat(op).take(len)
                            })
                            .collect()
                    })
                    .collect();
                let band = read_type.band_width(draft.len());
                hmm.polish_until_converge_with(&draft, &seqs, &mut ops, band)
                // kiley::polish_chunk_by_parts(&draft, &seqs, &config)
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
                    let mut seqs: Vec<_> = value
                        .iter()
                        .map(|(ed, is_forward)| match is_forward {
                            true => ed.label().to_vec(),
                            false => bio_utils::revcmp(ed.label()),
                        })
                        .filter(|e| !e.is_empty())
                        .collect();
                    let consensus = if seqs.len() > 3 {
                        let draft = seqs.pop().unwrap();
                        let band = read_type.band_width(draft.len());
                        hmm.polish_until_converge(&draft, &seqs, band)
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
    // selecting a branche from branches by using phaing sets.
    // In other words, whenever faced with branches,
    // the algorithm searches further node in the same phaseset.
    // If the search successed within `len` time BFS,
    // it select an arbitrary path from the current branching point to the found node,
    // removing other unnessesary edges.
    pub fn resolve_by_phaseset(&mut self, phaseset: &[HashSet<Node>], len: usize) -> Option<()> {
        let mut cand_nodes: Vec<_> = vec![];
        for node in self.nodes().filter(|n| n.copy_number == Some(1)) {
            cand_nodes.push((node.node, Position::Tail));
            cand_nodes.push((node.node, Position::Head));
        }
        // Make phase-set into node->phaseblock.
        let phaseset: HashMap<(u64, u64), usize> = phaseset
            .iter()
            .enumerate()
            .flat_map(|(block, nodes)| nodes.iter().map(|&node| (node, block)).collect::<Vec<_>>())
            .collect();
        for (from, from_pos) in cand_nodes {
            // Check if this is branching.
            let have_sibs = self
                .get_edges(from, from_pos)
                .all(|edge| 1 < self.get_edges(edge.to, edge.to_position).count());
            let num_edges = self.get_edges(from, from_pos).count();
            if num_edges < 2 || !have_sibs {
                continue;
            }
            // resolving branches.
            let from = (from, from_pos);
            let path_phasing = self.bfs_to_the_same_phase(from, len, &phaseset)?;
            let to = *path_phasing.last().unwrap();
            let (label, proxying) = self.spell_along_path(&path_phasing, from);
            debug!("PHASEPATH\tSpan\t{:?}\t{:?}\t{}", from.0, to.0, label.len());
            self.span_region(from, to, &path_phasing, label, proxying);
        }
        Some(())
    }
    // POINTER: ASSEMBLEIMPL
    // TODO: Tune this.
    // TODO: Maybe we should turn the optimization part into gibbs sampling...?
    pub fn clean_up_graph_for_assemble(
        &mut self,
        cov: f64,
        lens: &[usize],
        reads: &[&EncodedRead],
        c: &super::AssembleConfig,
        read_type: definitions::ReadType,
    ) {
        self.assign_copy_number(cov, lens);
        self.remove_zero_copy_elements(lens, 0.2);
        self.assign_copy_number(cov, lens);
        self.remove_zero_copy_elements(lens, 0.5);
        if read_type == ReadType::CLR {
            self.zip_up_overclustering();
        }
        // From good Likelihood ratio focus, to weaker ones.
        let min_llr = c.span_likelihood_ratio;
        let llr_stream = (0..10).rev().map(|x| x as f64).take_while(|&x| min_llr < x);
        for (i, llr) in llr_stream.enumerate() {
            debug!("REPEATRESOLVE\t{}", i);
            self.assign_copy_number_mcmc(cov, lens);
            self.resolve_repeats(reads, c, llr as f64);
            if log_enabled!(log::Level::Trace) {
                dump(self, i, c);
            }
        }
        self.remove_lightweight_edges(0, true);
        // if read_type == ReadType::CLR {
        // self.zip_up_overclustering();
        // }
        self.assign_copy_number_mcmc(cov, lens);
        self.remove_zero_copy_elements(lens, 0.8);
        self.remove_zero_copy_path(0.3);
    }
    pub fn new<R: std::borrow::Borrow<EncodedRead>>(
        reads: &'a [R],
        units: Option<&[Unit]>,
        read_type: definitions::ReadType,
        c: &AssembleConfig,
    ) -> Self {
        // Take a consensus of nodes and edges.
        let (nodes_seq, edge_seq) = take_consensus(reads, units, &read_type, c);
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
        for read in reads.iter().map(|r| r.borrow()) {
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

fn dump(graph: &DitchGraph, i: usize, c: &AssembleConfig) {
    let (segments, edge, group, summaries) = graph.spell(c);
    let nodes = segments.into_iter().map(|node| {
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
                log::debug!(
                    "CONUNIT\t{}\t{}\t{}",
                    contigsummary.id,
                    total / contigsummary.summary.len(),
                    ids.join("\t")
                );
                vec![coverage]
            }
            None => Vec::new(),
        };
        gfa::Record::from_contents(gfa::Content::Seg(node), tags)
    });
    let edges = edge
        .into_iter()
        .map(|(edge, tags)| gfa::Record::from_contents(gfa::Content::Edge(edge), tags));
    let group = gfa::Record::from_contents(gfa::Content::Group(group), vec![]);
    let header = gfa::Content::Header(gfa::Header::default());
    let header = gfa::Record::from_contents(header, vec![]);
    let records = std::iter::once(header)
        .chain(std::iter::once(group))
        .chain(nodes)
        .chain(edges)
        .collect();
    let gfa = gfa::GFA::from_records(records);
    if let Ok(mut wtr) = std::fs::File::create(format!("{}.gfa", i)).map(std::io::BufWriter::new) {
        use std::io::Write;
        if let Err(why) = writeln!(&mut wtr, "{}", gfa) {
            trace!("{:?}", why);
        }
    }
}

impl<'a> DitchGraph<'a> {
    pub fn dump(&self) {
        if log::log_enabled!(log::Level::Trace) {
            for edge in self.edges() {
                trace!("{}", edge.occ);
            }
        }
    }
    pub fn nodes(&self) -> impl std::iter::Iterator<Item = &DitchNode> {
        self.nodes.values()
    }
    pub fn edges(&self) -> impl std::iter::Iterator<Item = &DitchEdge> {
        self.nodes.values().flat_map(|n| n.edges.iter())
    }
    // Return iterator yeilding the edge from a specified node and position.
    fn get_edges(&self, node: Node, pos: Position) -> impl std::iter::Iterator<Item = &DitchEdge> {
        match self.nodes.get(&node) {
            Some(node) => node.edges.iter().filter(move |e| e.from_position == pos),
            None => panic!(
                "`get_edge({:?})` is called, but there's no such node,",
                node
            ),
        }
    }
    // Decrease the copy number of the given edge by one,
    // tune the occurence of these edges appropriately.
    fn decrement_edge_copy_number(&mut self, from: (Node, Position), to: (Node, Position)) {
        if let Some(node) = self.nodes.get_mut(&from.0) {
            node.edges
                .iter_mut()
                .filter(|edge| edge.from_position == from.1 && (edge.to, edge.to_position) == to)
                .for_each(|edge| edge.reduce_one_copy_number());
        }
        if let Some(node) = self.nodes.get_mut(&to.0) {
            node.edges
                .iter_mut()
                .filter(|edge| edge.from_position == to.1 && (edge.to, edge.to_position) == from)
                .for_each(|edge| edge.reduce_one_copy_number());
        }
    }
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
            for position in [Position::Head, Position::Tail] {
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
            for position in [Position::Head, Position::Tail] {
                let grand_child = node
                    .edges
                    .iter()
                    .find(|e| e.from_position == position)
                    .map(|e| {
                        let count = self.get_edges(e.to, e.to_position).count();
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
    // Partition edges whether or not it is in simple path/not.
    // Only (from, _, to, _) with from <= to edges would be in these vectors.
    fn partition_edges_by_simple_path(&self) -> (Vec<&DitchEdge>, Vec<&DitchEdge>) {
        let degrees: HashMap<_, (u32, u32)> =
            self.nodes
                .iter()
                .map(|(&key, node)| {
                    let (head, tail) = node.edges.iter().fold((0, 0), |(head, tail), edge| {
                        match edge.from_position {
                            Position::Head => (head + 1, tail),
                            Position::Tail => (head, tail + 1),
                        }
                    });
                    (key, (head, tail))
                })
                .collect();
        self.nodes
            .values()
            .flat_map(|node| node.edges.iter())
            .filter(|edge| edge.from <= edge.to)
            .partition(|edge| {
                let from_deg = match edge.from_position {
                    Position::Head => degrees[&edge.from].0,
                    Position::Tail => degrees[&edge.from].1,
                };
                let to_deg = match edge.to_position {
                    Position::Head => degrees[&edge.to].0,
                    Position::Tail => degrees[&edge.to].1,
                };
                from_deg == 1 && to_deg == 1
            })
    }
    // Return Node->serialized simple-path id hashmapping and
    // edges between simple paths.
    fn reduce_simple_path(&self) -> (HashMap<Node, usize>, Vec<&DitchEdge>) {
        let node_index: HashMap<Node, usize> = self
            .nodes
            .keys()
            .enumerate()
            .map(|(i, k)| (*k, i))
            .collect();
        let (edges_in_simple_path, edges_between_simple_path) =
            self.partition_edges_by_simple_path();
        use crate::find_union::FindUnion;
        let mut fu = FindUnion::new(node_index.len());
        for edge in edges_in_simple_path.iter() {
            let mut nodes = std::iter::once(edge.from)
                .chain(std::iter::once(edge.to))
                .map(|node| node_index[&node]);
            // Never panic
            let mut current = nodes.next().unwrap();
            for next in nodes {
                fu.unite(current, next);
                current = next;
            }
        }
        let cluster_index: HashMap<_, _> = (0..node_index.len())
            .filter(|&n| fu.find(n) == Some(n))
            .enumerate()
            .map(|(idx, key)| (key, idx))
            .collect();
        {
            let mut plug_num = vec![HashSet::new(); cluster_index.len()];
            for edge in edges_between_simple_path.iter() {
                let from_index = cluster_index[&fu.find(node_index[&edge.from]).unwrap()];
                plug_num[from_index].insert((edge.from, edge.from_position));
                let to_index = cluster_index[&fu.find(node_index[&edge.to]).unwrap()];
                plug_num[to_index].insert((edge.to, edge.to_position));
            }
            for plugs in plug_num {
                assert!(plugs.len() <= 2, "{:?}", plugs);
            }
        }
        let node_to_cluster: HashMap<_, _> = node_index
            .iter()
            .map(|(&node, &index)| {
                let cluster = fu.find(index).unwrap();
                (node, cluster_index[&cluster])
            })
            .collect();
        (node_to_cluster, edges_between_simple_path)
    }
    fn convert_connecting_edges(
        &self,
        node_to_pathid: &HashMap<Node, usize>,
        edges: &[&DitchEdge],
        calibrator: Option<&CoverageCalibrator>,
    ) -> (Vec<Vec<(Node, Position)>>, Vec<EdgeBetweenSimplePath>) {
        let path_num: usize = *node_to_pathid.values().max().unwrap_or(&0) + 1;
        let mut terminals: Vec<_> = vec![vec![]; path_num];
        let edges: Vec<_> = edges
            .iter()
            .map(|edge| {
                let from_index = node_to_pathid[&edge.from];
                let to_index = node_to_pathid[&edge.to];
                // Decide the "position" of these nodes.
                // There should not be more than two plugs on each simple path.
                let from_node = (edge.from, edge.from_position);
                let fp = match terminals[from_index].iter().position(|n| n == &from_node) {
                    Some(idx) => idx == 1,
                    None => {
                        terminals[from_index].push(from_node);
                        terminals[from_index].len() == 2
                    }
                };
                assert!(
                    terminals[from_index].len() <= 2,
                    "{:?}",
                    terminals[from_index]
                );
                let to_node = (edge.to, edge.to_position);
                let tp = match terminals[to_index].iter().position(|n| n == &to_node) {
                    Some(idx) => idx == 1,
                    None => {
                        terminals[to_index].push(to_node);
                        terminals[to_index].len() == 2
                    }
                };
                assert!(terminals[to_index].len() <= 2);
                let unit_len =
                    self.nodes[&edge.from].seq().len() + self.nodes[&edge.to].seq().len();
                let gap_len = (unit_len as i64 + edge.seq.len()) as usize;
                let weight = match calibrator {
                    Some(calibrator) => calibrator.calib(edge.occ, gap_len),
                    None => edge.occ as f64,
                };
                (from_index, fp, to_index, tp, weight)
            })
            .collect();
        (terminals, edges)
    }
    fn convert_path_weight(
        &self,
        node_to_pathid: &HashMap<Node, usize>,
        calibrator: Option<&CoverageCalibrator>,
    ) -> Vec<(f64, usize)> {
        let path_num: usize = *node_to_pathid.values().max().unwrap_or(&0) + 1;
        let mut node_weights = vec![(0f64, 0usize); path_num];
        for node in self.nodes.values() {
            let simple_path_index = node_to_pathid[&node.node];
            node_weights[simple_path_index].1 += 1;
            let unit_len = node.seq().len();
            node_weights[simple_path_index].0 += match calibrator {
                Some(calibrator) => calibrator.calib(node.occ, unit_len),
                None => node.occ as f64,
            };
        }
        node_weights
            .iter()
            .map(|&(sum, len)| (sum / len as f64, len))
            .collect()
    }
    fn gather_answer(
        &self,
        edges: &[(usize, bool, usize, bool, f64)],
        node_cp: &[usize],
        edge_cp: &[usize],
        node_to_pathid: &HashMap<Node, usize>,
        terminals: &[Vec<(Node, Position)>],
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        //let mut node_copy_number: HashMap<_, _> = node_to_pathid
        let node_copy_number: HashMap<_, _> = node_to_pathid
            .iter()
            .map(|(&node, &pathid)| (node, node_cp[pathid]))
            .collect();
        let mut edge_copy_number = HashMap::new();
        for (&edge_cp, &(from, fplus, to, tplus, _)) in edge_cp.iter().zip(edges.iter()) {
            let from = terminals[from][fplus as usize];
            let to = terminals[to][tplus as usize];
            edge_copy_number.insert((from, to), edge_cp);
            edge_copy_number.insert((to, from), edge_cp);
        }
        for edge in self.nodes.values().flat_map(|node| node.edges.iter()) {
            let from = node_to_pathid[&edge.from];
            let to = node_to_pathid[&edge.to];
            if from == to {
                let cp = node_cp[from];
                let from = (edge.from, edge.from_position);
                let to = (edge.to, edge.to_position);
                edge_copy_number.insert((from, to), cp);
                // Do we really need this? If a node is in a proxying edge,
                // it is already exhausted, and its copy number should not affect the
                // copy number of the original node.
                // Increment node in proxying edge.
                // Note that we should avoid double-count of the same edge.
                // if (edge.from, edge.from_position) <= (edge.to, edge.to_position) {
                //     for node in edge.proxying.iter().map(|x| x.0) {
                //         *node_copy_number.entry(node).or_default() += cp;
                //     }
                // }
            }
        }
        (node_copy_number, edge_copy_number)
    }
    fn generate_coverage_calib(&self, naive_cov: f64, lens: &[usize]) -> (CoverageCalibrator, f64) {
        let calibrator = CoverageCalibrator::new(lens);
        let unit_len_sum: usize = self.nodes.values().map(|n| n.seq().len()).sum();
        let cov = calibrator.calib_f64(naive_cov, unit_len_sum / self.nodes.len());
        debug!("COPYNUM\tCOVERAGE\t{:.3}\t{:.3}", naive_cov, cov,);
        (calibrator, cov)
    }
    /// Estimoate copy number of nodes and edges.
    /// *This function does not modify the graph content*. If you want
    /// to assign copy number to each node, call `assign_copy_number` instead.
    pub fn copy_number_estimation(
        &self,
        naive_cov: f64,
        lens: &[usize],
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (calibrator, cov) = self.generate_coverage_calib(naive_cov, lens);
        let calibrator = Some(&calibrator);
        let (terminals, mut edges) =
            self.convert_connecting_edges(&node_to_pathid, &connecting_edges, calibrator);
        let mut nodes = self.convert_path_weight(&node_to_pathid, calibrator);
        edges.iter_mut().for_each(|x| x.4 /= cov);
        nodes.iter_mut().for_each(|x| x.0 /= cov);
        let (node_cp, edge_cp) = super::copy_number::estimate_copy_number(&nodes, &edges);
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
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
    /// Estimoate copy number of nodes and edges by a gibbs sampler.
    /// *This function does not modify the graph content*.
    /// If you want to assign copy number to each node, call `assign_copy_number_gbs` instead.
    pub fn copy_number_estimation_gbs(
        &self,
        naive_cov: f64,
        lens: &[usize],
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (calibrator, cov) = self.generate_coverage_calib(naive_cov, lens);
        let calibrator = Some(&calibrator);
        let (terminals, edges) =
            self.convert_connecting_edges(&node_to_pathid, &connecting_edges, calibrator);
        let nodes = self.convert_path_weight(&node_to_pathid, calibrator);
        let nodes: Vec<_> = nodes.iter().map(|x| x.0).collect();
        let (node_cp, edge_cp) = super::copy_number::estimate_copy_number_gbs(&nodes, &edges, cov);
        if log_enabled!(log::Level::Trace) {
            trace!("COVCP\tType\tCov\tCp");
            for (n, cp) in nodes.iter().zip(node_cp.iter()) {
                trace!("COVCP\tNODE\t{:.2}\t{}", n, cp,);
            }
            for (e, cp) in edges.iter().zip(edge_cp.iter()) {
                trace!("COVCP\tEDGE\t{:.2}\t{}", e.4, cp);
            }
        }
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number_gbs(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) = self.copy_number_estimation_gbs(naive_cov, lens);
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
    /// Estimoate copy number of nodes and edges by MCMC.
    /// *This function does not modify the graph content*.
    /// If you want to assign copy number to each node, call `assign_copy_number_gbs` instead.
    pub fn copy_number_estimation_mcmc(
        &self,
        naive_cov: f64,
        _lens: &[usize],
    ) -> (HashMap<Node, usize>, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (terminals, edges) =
            self.convert_connecting_edges(&node_to_pathid, &connecting_edges, None);
        let nodes = self.convert_path_weight(&node_to_pathid, None);
        let (node_cp, edge_cp) =
            super::copy_number::estimate_copy_number_mcmc(&nodes, &edges, naive_cov);
        if log_enabled!(log::Level::Trace) {
            trace!("COVCP\tType\tCov\tCp");
            for ((n, len), cp) in nodes.iter().zip(node_cp.iter()) {
                trace!("COVCP\tNODE\t{}\t{}\t{}", n, len, cp,);
            }
            for (e, cp) in edges.iter().zip(edge_cp.iter()) {
                trace!("COVCP\tEDGE\t{}\t{}", e.4, cp);
            }
        }
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }
    // TODO:Fasten this function.
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number_mcmc(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) =
            self.copy_number_estimation_mcmc(naive_cov, lens);
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
    pub fn remove_zero_copy_elements(&mut self, lens: &[usize], thr: f64) {
        // (from,from_position,to,to_position) and from.0 <= to.0.
        fn format_edge(e: &DitchEdge) -> (Node, Position, Node, Position) {
            if (e.from, e.from_position) <= (e.to, e.to_position) {
                (e.from, e.from_position, e.to, e.to_position)
            } else {
                (e.to, e.to_position, e.from, e.from_position)
            }
        }
        // Check the node violating for right-left condition.
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
                match (plus, minus) {
                    (0, 0) => Some(key),
                    (0, _) | (_, 0) => None,
                    (x, y) if x != y => Some(key),
                    _ => None,
                }
            })
            .collect();
        debug!("UNSOUND\t{}\t{}", unsound_nodes.len(), self.nodes.len());
        let calibrator = CoverageCalibrator::new(lens);
        let ave_unit_len = {
            let sum_unit_len: usize = self.nodes.values().map(|node| node.seq().len()).sum();
            sum_unit_len / self.nodes.len()
        };
        let calib_occ = |edge: &DitchEdge| {
            let len = (edge.seq.len() + 2 * ave_unit_len as i64) as usize;
            calibrator.calib(edge.occ, len)
        };
        let mut is_ok_to_remove = HashSet::new();
        let mut retain_edges = HashSet::new();
        for (key, node) in self.nodes.iter() {
            // If the estimation of the copy number is is poor, do not remove edges.
            if unsound_nodes.contains(key) {
                retain_edges.extend(node.edges.iter().map(format_edge));
                continue;
            }
            for position in [Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == position);
                // It is ok to take max by this way because if the max is 0, there's no edge.
                let max = edges.clone().map(calib_occ).fold(0f64, |x, y| x.max(y));
                // It automatically drop the heviest edge from the selection.
                for edge in edges {
                    if matches!(edge.copy_number, Some(0)) && calib_occ(edge) / max < thr {
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
        debug!("REMOVED\t{}\t{}", removed, max_removed_weight);
        // If the condition 2. and 3. hold, then the degree should be 0.
        self.nodes.retain(|_, node| {
            node.copy_number.map(|cp| cp != 0).unwrap_or(true) || !node.edges.is_empty()
        });
    }
    /// Remove zero-copy path.
    /// 1. Find a zero-copy path.
    /// 2. Check whether the destination of the zero-copy path are reachable from other paths.
    /// 3. If so, remove that path.
    pub fn remove_zero_copy_path(&mut self, thr: f64) {
        // enumerate parent nodes of the single-copy-path.
        let mut parents_of_zc = vec![];
        for (&key, node) in self.nodes.iter() {
            // Parents should have non-zero-copy.
            if !matches!(node.copy_number,Some(x) if 0 < x ) {
                continue;
            }
            for pos in [Position::Head, Position::Tail] {
                let eds = node.edges.iter().filter(|e| e.from_position == pos);
                if 2 <= eds.clone().count() {
                    let has_zc = eds
                        .filter_map(|e| self.nodes.get(&e.to))
                        .any(|n| matches!(n.copy_number, Some(0)));
                    if has_zc {
                        parents_of_zc.push((key, pos));
                    }
                }
            }
        }
        // Check reachability.
        for (parent, par_pos) in parents_of_zc {
            if !self.nodes.contains_key(&parent) {
                continue;
            }
            if self.get_edges(parent, par_pos).count() <= 1 {
                continue;
            }
            let has_zc = self
                .get_edges(parent, par_pos)
                .filter_map(|e| self.nodes.get(&e.to).and_then(|n| n.copy_number))
                .any(|cp| cp == 0);
            if has_zc {
                // Check the destination of the zero-copy path and non-zero-copy paths.
                let (zc_edges, nzc_edges): (Vec<_>, Vec<_>) =
                    self.get_edges(parent, par_pos).partition(|e| {
                        matches!(self.nodes.get(&e.to).and_then(|n| n.copy_number), Some(0))
                    });
                let zc_dests: HashSet<_> = zc_edges
                    .iter()
                    .flat_map(|e| self.simple_path_and_dest(e.to, e.to_position).1)
                    .collect();
                let zc_max = zc_edges.iter().map(|e| self.nodes[&e.to].occ).max();
                let nzc_dests: HashSet<_> = nzc_edges
                    .iter()
                    .flat_map(|e| self.simple_path_and_dest(e.to, e.to_position).1)
                    .collect();
                let nzc_max = nzc_edges.iter().map(|e| self.nodes[&e.to].occ).max();
                // zc_max is Some(_), because there's at least one zc_edge.
                let zc_nzc_ratio = match (zc_max, nzc_max) {
                    (Some(x), Some(y)) => x as f64 / y as f64,
                    (Some(_), _) => 1f64,
                    _ => panic!(),
                };
                if zc_dests.is_subset(&nzc_dests) && zc_nzc_ratio < thr {
                    debug!(
                        "RZCP\t{:?}\t{}\t{}\t{:.3}",
                        parent,
                        par_pos,
                        zc_dests.is_subset(&nzc_dests),
                        zc_nzc_ratio
                    );
                    let targets: Vec<_> = zc_edges.iter().map(|e| (e.to, e.to_position)).collect();
                    for (node, pos) in targets {
                        self.remove_edge_and_pruning(parent, par_pos, node, pos);
                    }
                }
            }
        }
    }
    /// Remove transitive edge if the copy_number of that edge is zero.
    #[allow(dead_code)]
    pub fn transitive_edge_reduction(&mut self) {
        let mut removed_edges: HashMap<_, Vec<_>> =
            self.nodes.keys().map(|&k| (k, vec![])).collect();
        for (from, node) in self.nodes.iter() {
            for &pos in &[Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                let to_cands = edges.clone().count();
                if to_cands > 1 {
                    for e in edges.clone() {
                        // Check whether or not this edge is transitive edge.
                        let is_transitive =
                            edges.clone().any(|edge| match self.nodes.get(&edge.to) {
                                Some(child) => child
                                    .edges
                                    .iter()
                                    .filter(|gc| gc.from_position == !edge.to_position)
                                    .any(|gc| gc.to == e.to && gc.to_position == e.to_position),
                                None => false,
                            });
                        if is_transitive && matches!(e.copy_number, Some(0)) {
                            let (to, t_pos) = (e.to, e.to_position);
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
        }
        self.nodes.iter_mut().for_each(|(key, node)| {
            let to_be_removed = &removed_edges[key];
            node.edges.retain(|x| {
                let probe = (x.from_position, x.to, x.to_position);
                !to_be_removed.contains(&probe)
            })
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
                (*node, (is_artic, are_bridge))
            })
            .collect()
    }
    /// Zip up overclustered regions.
    /// The graph should have estimeted repeat numbers on each node.
    pub fn zip_up_overclustering(&mut self) {
        let mut to_remove = HashSet::new();
        for node in self.nodes.values() {
            if node.copy_number.map(|cp| cp != 1).unwrap_or(true) {
                continue;
            }
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
                if first_dest == second_dest {
                    debug!("ZIPPINGUP\t{:?}\t{}", node.node, pos);
                    let simple_path = self.simple_path_from(edges.clone().next().unwrap());
                    to_remove.extend(simple_path);
                }
            }
        }
        self.remove_nodes(&to_remove);
    }
    /// Return the position where a simple path from given edge ends.
    pub fn destination(&self, edge: &DitchEdge) -> (Node, Position) {
        let mut current_node = edge.to;
        let mut current_pos = edge.to_position;
        loop {
            // Check if there's more than two indegree/outdegree.
            let (indeg, outdeg) =
                self.nodes[&current_node]
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
            let traced_edge = self.get_edges(current_node, current_pos).next().unwrap();
            current_node = traced_edge.to;
            current_pos = traced_edge.to_position;
        }
        (current_node, current_pos)
    }
    /// Return simple path from the given edge to ends.
    /// Note that a node is added when the node is consumed in the simple path, i.e., we moved
    /// from a position to the opposite position.
    pub fn simple_path_from(&self, edge: &DitchEdge) -> Vec<Node> {
        let mut current_node = edge.to;
        let mut current_pos = edge.to_position;
        let mut nodes = vec![];
        loop {
            // Check if there's more than two indegree/outdegree.
            let (indeg, outdeg) =
                self.nodes[&current_node]
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
            let traced_edge = self.get_edges(current_node, current_pos).next().unwrap();
            current_node = traced_edge.to;
            current_pos = traced_edge.to_position;
        }
        nodes
    }
    /// Return simple path start from the given node and position, and the destination nodes after this simple path.
    pub fn simple_path_and_dest(
        &self,
        mut node: Node,
        mut position: Position,
    ) -> (Vec<Node>, Vec<Node>) {
        let mut nodes = vec![];
        loop {
            // Move position.
            position = !position;
            nodes.push(node);
            let mut edges = self.get_edges(node, position);
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
            let indeg = self.get_edges(edge.to, edge.to_position).count();
            assert!(1 <= indeg);
            if 1 < indeg {
                break;
            }
            // Move node.
            node = edge.to;
            position = edge.to_position;
            nodes.extend(edge.proxying.iter().map(|(n, _, _)| n));
        }
        let mut dests: Vec<_> = self.get_edges(node, position).map(|e| e.to).collect();
        dests.sort_unstable();
        (nodes, dests)
    }

    /// Remove small tips.
    /// In this function, for all node with the degree of zero,
    /// we speculate the *local* coverage by traversing DIAG distance.
    /// If the coverage is less than `cov * thr`, the node would be removed.
    pub fn remove_tips(&mut self, thr: f64, diag: usize) {
        let to_remove: HashSet<_> = self
            .nodes
            .iter()
            .filter_map(|(key, node)| {
                for &pos in &[Position::Head, Position::Tail] {
                    if node.edges.iter().filter(|e| e.from_position == pos).count() == 0 {
                        let coverage = self.local_coverage(key, pos, diag);
                        if (node.occ as f64) < coverage * thr {
                            return Some(*key);
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
    fn local_coverage(&self, node: &Node, pos: Position, diag: usize) -> f64 {
        let (mut total_cov, mut total_node) = (0, 0);
        let mut current = vec![(node, pos)];
        for _ in 0..diag {
            let mut next_nodes = vec![];
            for &(node, pos) in current.iter() {
                if let Some(node) = self.nodes.get(node) {
                    total_cov += node.occ;
                    total_node += 1;
                    for edge in node.edges.iter().filter(|e| e.from_position == !pos) {
                        next_nodes.push((&edge.to, edge.to_position));
                    }
                }
            }
            next_nodes.sort();
            next_nodes.dedup();
            current = next_nodes;
        }
        total_cov as f64 / total_node as f64
    }
    // Removing nodes and re-map all of the indices.
    fn remove_nodes(&mut self, to_remove: &HashSet<Node>) {
        self.nodes.retain(|key, _| !to_remove.contains(key));
        for node in self.nodes.values_mut() {
            node.edges.retain(|edge| !to_remove.contains(&edge.to));
        }
    }
    /// Check the dynamic of foci.
    pub fn survey_foci(&mut self, foci: &[Focus]) -> usize {
        foci.iter()
            .filter(|focus| {
                let (key, pos) = (focus.from, focus.from_position);
                // Check if the targets are present in the graph.
                let from_copy_num = self.nodes.get(&focus.from).and_then(|n| n.copy_number);
                let to_copy_num = self.nodes.get(&focus.to).and_then(|node| node.copy_number);
                if !matches!(from_copy_num, Some(1)) || !matches!(to_copy_num, Some(1)) {
                    return false;
                }
                // Check if this is branching.
                let num_edges = self.get_edges(key, pos).count();
                if num_edges != 1 {
                    return false;
                }
                let edge = self.get_edges(key, pos).next().unwrap();
                let sibs = self.get_edges(edge.to, edge.to_position).count();
                if sibs <= 1 {
                    return false;
                }
                debug!("FOCUS\tTRY\t{}", focus);
                self.survey_focus(focus).is_some()
            })
            .count()
    }
    // If resolveing successed, return Some(())
    // othewise, it encountered a stronger focus, and return None.
    // fn survey_focus(&mut self, focus: &Focus, info: &HashMap<Node, f64>) -> Option<()> {
    fn survey_focus(&mut self, focus: &Focus) -> Option<()> {
        let path_to_focus = self.bfs_to_the_target(focus)?;
        let from = (focus.from, focus.from_position);
        let to = (focus.to, focus.to_position);
        let (label, proxying) = self.spell_along_path(&path_to_focus, from);
        debug!("FOCUS\tSPAN\t{}\t{}bp", focus, label.len());
        self.span_region(from, to, &path_to_focus, label, proxying);
        Some(())
    }
    fn span_region(
        &mut self,
        (from, from_pos): (Node, Position),
        (to, to_pos): (Node, Position),
        path_spaning: &[(Node, Position)],
        label: EdgeLabel,
        proxying: Vec<(Node, bool, usize)>,
    ) {
        assert_eq!(*path_spaning.last().unwrap(), (to, to_pos));
        let mut edge = DitchEdge::with_proxy(from, from_pos, to, to_pos, label, proxying);
        edge.occ = self.nodes[&from].edges.iter().map(|e| e.occ).sum();
        // First, allocate edge to be removed.
        // Note that sometimes the newly added edge is already in this node(i.e. resolving some short branch),
        // So we should be care that the edge is not in the removing edges.
        let remove_from: Vec<_> = self
            .get_edges(from, from_pos)
            .filter(|e| !(e.to == to && e.to_position == to_pos))
            .map(|e| (e.to, e.to_position))
            .collect();
        let remove_to: Vec<_> = self
            .get_edges(to, to_pos)
            .filter(|e| !(e.to == from && e.to_position == from_pos))
            .map(|e| (e.to, e.to_position))
            .collect();
        // Next, add the new edge. If there's already the same egdge, just update it.
        if let Some(node) = self.nodes.get_mut(&to) {
            let edge = edge.reverse();
            match node.edges.iter_mut().find(|e| e == &&edge) {
                Some(e) => *e = edge,
                None => node.edges.push(edge),
            }
        }
        if let Some(node) = self.nodes.get_mut(&from) {
            match node.edges.iter_mut().find(|e| e == &&edge) {
                Some(e) => *e = edge,
                None => node.edges.push(edge),
            }
        }
        // Next, removing the rest of the edges.
        for (node, pos) in remove_from {
            self.remove_edge_and_pruning(from, from_pos, node, pos);
        }
        for (node, pos) in remove_to {
            self.remove_edge_and_pruning(to, to_pos, node, pos);
        }
        // Removing nodes along path, if nessesarry.
        for &(node, pos) in path_spaning {
            // Sometimes, this node is already removed from the graph.
            if !self.nodes.contains_key(&node) {
                continue;
            }
            // If it is not zero-copy node, do not anything.
            if !matches!(self.nodes[&node].copy_number, Some(0)) {
                continue;
            }
            let parents: Vec<_> = self
                .get_edges(node, pos)
                .map(|e| (e.to, e.to_position))
                .collect();
            for (parent, par_pos) in parents {
                if self
                    .get_edges(parent, par_pos)
                    .any(|e| !(e.to == node && e.to_position == pos))
                {
                    self.remove_edge_and_pruning(parent, par_pos, node, pos);
                }
            }
        }
    }
    // Fron (start,pos) position to the node with the same phase block.
    // Note that, the path would contain the (target,target pos) at the end of the path and
    // not contain the (start,pos) position itself.
    fn bfs_to_the_same_phase(
        &self,
        (start, pos): (Node, Position),
        len: usize,
        phaseset: &HashMap<Node, usize>,
    ) -> Option<Vec<(Node, Position)>> {
        // Breadth first search until get to the target.
        // Just !pos to make code tidy.
        let target_phase = Some(phaseset.get(&start)?);
        let mut nodes_at = vec![vec![(start, !pos)]];
        let mut parents = vec![vec![]];
        'outer: for dist in 0..len {
            let mut next_nodes = vec![];
            let mut parent = vec![];
            assert_eq!(dist + 1, nodes_at.len());
            for (idx, &(node, pos)) in nodes_at[dist].iter().enumerate() {
                if phaseset.get(&node) == target_phase && node != start {
                    break 'outer;
                }
                // Move to the end of the node.
                let pos = !pos;
                for edge in self.get_edges(node, pos) {
                    next_nodes.push((edge.to, edge.to_position));
                    parent.push(idx);
                }
            }
            assert_eq!(next_nodes.len(), parent.len());
            nodes_at.push(next_nodes);
            parents.push(parent);
        }
        // Back-track.
        assert_eq!(nodes_at.len(), parents.len());
        let mut dist = nodes_at.len() - 1;
        let mut idx = match nodes_at[dist]
            .iter()
            .position(|(n, _)| phaseset.get(n) == target_phase)
        {
            Some(idx) => idx,
            None => {
                debug!("PHASEPATH\tFailedToFindTarget\t{:?}", start);
                return None;
            }
        };
        let mut back_track = vec![];
        while 0 < dist {
            back_track.push(nodes_at[dist][idx]);
            idx = parents[dist][idx];
            dist -= 1;
        }
        back_track.reverse();
        for (i, (node, pos)) in back_track.iter().enumerate() {
            debug!("PHASEPATH\tTrace\t{}\t{:?}\t{}", i, node, pos);
        }
        Some(back_track)
    }
    // Fron (start,pos) position to the focus target(start=focus.node, position=focus.position).
    // Note that, the path would contain the (target,target pos) at the end of the path and
    // not contain the (start,pos) position itself.
    // TODO: Maybe this function failed to find the target focus,
    // because of the removed nodes.
    // If it happened, eigther this focus or previous
    // branching is false. Maybe just dropping this focus would be better?
    fn bfs_to_the_target(&self, focus: &Focus) -> Option<Vec<(Node, Position)>> {
        // Breadth first search until get to the target.
        // Just !pos to make code tidy.
        let (start, pos) = (focus.from, focus.from_position);
        let mut nodes_at = vec![vec![(start, !pos)]];
        let mut parents = vec![vec![]];
        'outer: for dist in 0..focus.dist + 1 {
            let mut next_nodes = vec![];
            let mut parent = vec![];
            assert_eq!(dist + 1, nodes_at.len());
            for (idx, &(node, pos)) in nodes_at[dist].iter().enumerate() {
                if node == focus.to && pos == focus.to_position && focus.dist == dist {
                    break 'outer;
                }
                // Move to the end of the node.
                let pos = !pos;
                for edge in self.get_edges(node, pos) {
                    next_nodes.push((edge.to, edge.to_position));
                    parent.push(idx);
                }
            }
            assert_eq!(next_nodes.len(), parent.len());
            nodes_at.push(next_nodes);
            parents.push(parent);
        }
        // Back-track.
        assert_eq!(nodes_at.len(), parents.len());
        let mut dist = nodes_at.len() - 1;
        let mut idx = match nodes_at[dist]
            .iter()
            .position(|&(n, p)| n == focus.to && p == focus.to_position)
        {
            Some(idx) => idx,
            None => {
                debug!("FOCUS\tFailedToFindTarget{}", focus);
                return None;
            }
        };
        let mut back_track = vec![];
        while 0 < dist {
            back_track.push(nodes_at[dist][idx]);
            idx = parents[dist][idx];
            dist -= 1;
        }
        back_track.reverse();
        for (i, (node, pos)) in back_track.iter().enumerate() {
            debug!("FOCUS\tTrace\t{}\t{:?}\t{}", i, node, pos);
        }
        Some(back_track)
    }
    // Spell along the given path,
    // reducing the copy number of the nodes, the occs of the nodes, and the
    // occs of the edges.
    fn spell_along_path(
        &mut self,
        path: &[(Node, Position)],
        (start, pos): (Node, Position),
    ) -> (EdgeLabel, Vec<(Node, bool, usize)>) {
        let first_elm = path[0];
        let first_label = self
            .get_edges(start, pos)
            .find(|edge| (edge.to, edge.to_position) == first_elm)
            .map(|edge| edge.seq.clone())
            .unwrap();
        let mut path_label = vec![];
        let mut processed_nodes = vec![];
        for window in path.windows(2) {
            // !!!Here, each position is the start position of each node!!!!
            let (from, from_pos) = window[0];
            let (to, to_pos) = window[1];
            let mut node_seq = match from_pos {
                Position::Head => self.nodes[&from].seq().to_vec(),
                Position::Tail => bio_utils::revcmp(self.nodes[&to].seq()),
            };
            // Reduce the occ and copy number of the node, and occ of the edge.
            let occ = {
                // Reduce te occ of the edge.
                // The reverse direction should be cared.
                let reduced: Vec<_> = self
                    .nodes
                    .get(&from)
                    .unwrap()
                    .edges
                    .iter()
                    .filter(|e| e.to == to && e.to_position == to_pos)
                    .filter(|e| matches!(e.copy_number, Some(cp) if cp > 0))
                    .map(|edge| ((edge.from, edge.from_position), (edge.to, edge.to_position)))
                    .collect();
                for (from, to) in reduced {
                    self.decrement_edge_copy_number(from, to)
                }
                // This unwrap never panics.
                let node = self.nodes.get_mut(&from).unwrap();
                match node.copy_number {
                    Some(cp) if 0 < cp => {
                        if [536, 786].contains(&node.node.0) {
                            debug!("DUMP\t{}", node);
                        }
                        // TODO: Is this OK? Maybe the copy number is not reliable.
                        // The problem is not the reduction of the copy number -- it should be zero --
                        // but the reduction of the `occ`. Maybe we should substract average coverage
                        // along the path from the `node.occ`.
                        let occ = node.occ / cp;
                        node.occ -= occ;
                        node.copy_number = Some(cp - 1);
                        occ
                    }
                    _ => 0,
                }
            };
            processed_nodes.push((from, from_pos == Position::Head, occ));
            let label = self
                .get_edges(from, !from_pos)
                .find(|edge| edge.to == to && edge.to_position == to_pos)
                .map(|edge| &edge.seq)
                .unwrap();
            match label {
                EdgeLabel::Ovlp(l) => {
                    node_seq.reverse();
                    let truncate = (node_seq.len() as i64 + *l).max(0);
                    node_seq.truncate(truncate as usize);
                    node_seq.reverse();
                }
                EdgeLabel::Seq(seq) => node_seq.extend(seq),
            }
            path_label.extend(node_seq);
        }
        // Finalize the label of the edge.
        let label = match first_label {
            EdgeLabel::Ovlp(l) if l + (path_label.len() as i64) < 0 => {
                EdgeLabel::Ovlp(l + path_label.len() as i64)
            }
            EdgeLabel::Ovlp(l) => {
                path_label.reverse();
                for _ in 0..(-l) as usize {
                    path_label.pop();
                }
                path_label.reverse();
                EdgeLabel::Seq(path_label)
            }
            EdgeLabel::Seq(seq) => EdgeLabel::Seq(seq.iter().copied().chain(path_label).collect()),
        };
        (label, processed_nodes)
    }
    // Removing (from,from_pos)-(to,to_pos) edge.
    // Then, recursively removing the edges in the opposite position and `to` node itself
    // if the degree of `(to,to_pos)` is zero after removing.
    fn remove_edge_and_pruning(
        &mut self,
        from: Node,
        from_pos: Position,
        to: Node,
        to_pos: Position,
    ) {
        debug!(
            "Removing ({},{},{})-({},{},{})",
            from.0, from.1, from_pos, to.0, to.1, to_pos
        );
        if let Some(node) = self.nodes.get_mut(&from) {
            node.edges.retain(|e| {
                !(e.from_position == from_pos && e.to == to && e.to_position == to_pos)
            });
        }
        // Remove in anti-direction.
        if let Some(node) = self.nodes.get_mut(&to) {
            node.edges.retain(|e| {
                !(e.from_position == to_pos && e.to == from && e.to_position == from_pos)
            });
        }
        // If the degree became zero and this node and the copy number is zero, remove recursively
        let removed_all_edges = self
            .nodes
            .get(&to)
            .map(|node| node.edges.iter().all(|e| e.from_position == !to_pos))
            .unwrap_or(false);
        let to_copy_num = self.nodes.get(&to).and_then(|n| n.copy_number);
        let is_zero_copy = matches!(to_copy_num, Some(0));
        if removed_all_edges && is_zero_copy {
            // Removing this node.
            // First, recursively call the removing function.
            let remove_edges: Vec<_> = self.nodes[&to]
                .edges
                .iter()
                .map(|e| (e.to, e.to_position))
                .collect();
            for &(node, pos) in remove_edges.iter() {
                // The edge is from the *opposite* position of the `to` node!
                self.remove_edge_and_pruning(to, !to_pos, node, pos);
            }
            debug!("Removing {:?}", to);
            // Then, removing this node itself.
            if self.nodes.contains_key(&to) {
                assert!(self.nodes[&to].edges.is_empty());
                self.nodes.remove(&to);
            }
        }
    }
    /// Resolve repeats by "foci" algorithm.
    /// The algorithm consists three step.
    /// 1. For each branching node u, find another node v where reads outgoing from u are
    /// "gathering". In other words, we check the likelihood ratio between "gathering"-model
    /// and null-model and check whether or not the likelihood ration is larger than `thr`.
    /// 2. Select an arbitrary path u -> v, make an new edge connecting u and v,
    /// and remove all other edges from u and to v. By doing this, u and v would be in the same simple path(repeat resolved).
    /// 3. Recursively remove nodes and their edges, whose in-degree are zero.
    /// TODO: This code has bug when the graph contains loops?
    pub fn resolve_repeats(&mut self, reads: &[&EncodedRead], config: &AssembleConfig, thr: f64) {
        debug!("FOCI\tRESOLVE\t{:.3}\t{}", thr, config.min_span_reads);
        let mut count = 1;
        while count != 0 {
            let mut foci = self.get_foci(reads, config);
            let prev = foci.len();
            foci.retain(|val| thr < val.llr());
            foci.retain(|val| val.llr() != std::f64::INFINITY);
            if foci.is_empty() {
                break;
            }
            debug!("FOCI\tNUM\t{}\t{}", prev, foci.len());
            foci.sort_by(|x, y| match y.llr().partial_cmp(&x.llr()).unwrap() {
                std::cmp::Ordering::Equal => (y.dist, y.from).cmp(&(x.dist, x.from)),
                x => x,
            });
            count = self.survey_foci(&foci);
            debug!("FOCI\tTryAndSuccess\t{}\t{}", foci.len(), count);
        }
    }
    /// Return a hash map containing all foci, thresholded by `config` parameter.
    /// For each (node, position), keep the strongest focus.
    pub fn get_foci(&self, reads: &[&EncodedRead], config: &AssembleConfig) -> Vec<Focus> {
        let mut foci = vec![];
        for node in self
            .nodes
            .values()
            .filter(|n| matches!(n.copy_number, Some(1)))
        {
            for pos in [Position::Head, Position::Tail] {
                let edge_num = node
                    .edges
                    .iter()
                    .filter(|edge| edge.from_position == pos)
                    .count();
                if edge_num != 1 {
                    continue;
                }
                let edge = node.edges.iter().find(|e| e.from_position == pos).unwrap();
                let siblings = self.get_edges(edge.to, edge.to_position).count();
                assert!(1 <= siblings);
                if siblings <= 1 {
                    continue;
                }
                if let Some(focus) = self.examine_focus(node, pos, reads, config) {
                    foci.push(focus);
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
        config: &AssembleConfig,
    ) -> Option<Focus> {
        // TODO:Maybe we can faster this process by using reverse indexing.
        let reads: Vec<_> = reads
            .iter()
            .filter(|r| r.contains(node.node))
            .copied()
            .collect();
        let dist_nodes = self.enumerate_candidate_nodes(&reads, config.min_span_reads, node, pos);
        let mut focus: Option<Focus> = None;
        for (dist, nodes) in dist_nodes.iter().enumerate().filter(|(_, ns)| 1 < ns.len()) {
            // TODO: Maybe we need to clarify condition so that to remove this if-statement...
            let total_occs: usize = nodes.iter().map(|n| self.nodes[&n.0].occ).sum();
            if total_occs == 0 {
                continue;
            }
            let dist = dist + 1;
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
            let ith_ln = {
                let mut probs: Vec<_> = nodes.iter().map(|n| self.nodes[&n.0].occ as f64).collect();
                let sum: f64 = probs.iter().sum();
                probs.iter_mut().for_each(|x| *x = (*x / sum).ln());
                assert!(probs.iter().all(|x| !x.is_nan()), "{:?}", probs);
                probs
            };
            let null_prob: f64 = occs
                .iter()
                .zip(ith_ln.iter())
                .map(|(&occ, ln)| match occ {
                    0 => 0f64,
                    _ => occ as f64 * ln,
                })
                .sum();
            assert!(!null_prob.is_nan(), "{:?}\t{:?}", ith_ln, occs);
            assert!(null_prob < std::f64::INFINITY, "{:?}\t{:?}", ith_ln, occs);
            let choice_num = nodes.len() as f64;
            let correct_lk = ((1f64 - ERROR_PROB).powi(2) + ERROR_PROB / choice_num).ln();
            let error_lk = {
                let correct_to_error =
                    (1.0 - ERROR_PROB) * ERROR_PROB * (choice_num - 1.0) / choice_num;
                let error_to_error =
                    ERROR_PROB / choice_num * ERROR_PROB * (choice_num - 1f64) / choice_num;
                (correct_to_error + error_to_error).ln()
            };
            assert!(!correct_lk.is_nan());
            assert!(!error_lk.is_nan());
            if let Some((lk_ratio, (to, to_pos))) = nodes
                .iter()
                .enumerate()
                .filter(|(_, (node, _))| matches!(self.nodes[node].copy_number, Some(1)))
                .map(|(k, &n)| {
                    let lk: f64 = occs
                        .iter()
                        .map(|&x| x as f64)
                        .enumerate()
                        .map(|(i, occ)| occ * (if i == k { correct_lk } else { error_lk }))
                        .sum();
                    assert!(!lk.is_nan());
                    (lk - null_prob, n)
                })
                .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
            {
                if focus.as_ref().map(|f| f.llr()).unwrap_or(0f64) < lk_ratio {
                    focus = Some(Focus::new(node.node, pos, to, to_pos, dist, lk_ratio));
                }
            }
        }
        focus
    }
    pub fn enumerate_candidate_nodes(
        &self,
        reads: &[&EncodedRead],
        min_span: usize,
        node: &DitchNode,
        pos: Position,
    ) -> Vec<Vec<(Node, Position)>> {
        let max_reach = Self::max_reach(reads, min_span, node, pos);
        let mut nodes_at = vec![vec![]; max_reach];
        let mut active_nodes: Vec<_> = node
            .edges
            .iter()
            .filter(|edge| edge.from_position == pos)
            .map(|edge| (edge.proxying.len() + 1, edge.to, edge.to_position))
            .filter(|&(d, _, _)| d <= max_reach)
            .collect();
        while !active_nodes.is_empty() {
            let (dist, node, p) = active_nodes.pop().unwrap();
            let to_nodes = self
                .get_edges(node, !p)
                .map(|e| (dist + e.proxying.len() + 1, e.to, e.to_position))
                .filter(|&(d, _, _)| d <= max_reach);
            active_nodes.extend(to_nodes);
            nodes_at[dist - 1].push((node, p));
        }
        nodes_at.iter_mut().for_each(|nodes| {
            nodes.sort();
            nodes.dedup();
        });
        nodes_at
    }
    fn max_reach(
        reads: &[&EncodedRead],
        min_span: usize,
        node: &DitchNode,
        pos: Position,
    ) -> usize {
        let mut dist_node: HashMap<_, _> = HashMap::new();
        for read in reads.iter() {
            let start = read
                .nodes
                .iter()
                .position(|n| (n.unit, n.cluster) == node.node)
                .unwrap();
            match (read.nodes[start].is_forward, pos) {
                (true, Position::Tail) | (false, Position::Head) => {
                    // ---> (Here to start) -- <---> -- <---> .... Or
                    // <--- (Here to start) -- <---> -- <---> ....
                    for (dist, node) in read.nodes.iter().skip(start).enumerate() {
                        let node = (node.unit, node.cluster);
                        *dist_node.entry((node, dist)).or_default() += 1;
                    }
                }
                (true, Position::Head) | (false, Position::Tail) => {
                    let read = read.nodes.iter().take(start + 1).rev();
                    for (dist, node) in read.enumerate() {
                        let node = (node.unit, node.cluster);
                        *dist_node.entry((node, dist)).or_default() += 1;
                    }
                }
            }
        }
        dist_node
            .iter()
            .filter(|&(_, &count)| min_span <= count)
            .fold(0, |dist, (&(_, d), _)| d.max(dist))
    }
    // fn get_covered_nodes(reads: &[&EncodedRead], thr: usize) -> HashSet<Node> {
    //     let mut counts: HashMap<_, usize> = HashMap::new();
    //     for node in reads.iter().flat_map(|r| r.nodes.iter()) {
    //         *counts.entry((node.unit, node.cluster)).or_default() += 1;
    //     }
    //     counts.retain(|_, val| thr <= *val);
    //     counts.keys().copied().collect()
    // }
    // /// Return nodes achievable from the given node.
    // /// The i-th vector is the nodes reachable from `node` in i+1 hop.
    // /// Fot the i-th vector, one of the node shoule be in the
    // /// covered region.
    // pub fn enumerate_node_upto(
    //     &self,
    //     node: &DitchNode,
    //     pos: Position,
    //     covered: &HashSet<Node>,
    // ) -> Vec<Vec<(Node, Position)>> {
    //     let mut has_arrived = HashSet::new();
    //     let initial_nodes: Vec<_> = node
    //         .edges
    //         .iter()
    //         .filter(|edge| edge.from_position == pos)
    //         .map(|edge| (edge.to, edge.to_position))
    //         .collect();
    //     let mut nodes_at: Vec<Vec<_>> = vec![initial_nodes];
    //     for i in 0.. {
    //         let mut next_nodes = vec![];
    //         for &(node, p) in nodes_at[i].iter() {
    //             has_arrived.insert(node);
    //             let to_nodes = self
    //                 .get_edges(node, !p)
    //                 .filter(|e| !has_arrived.contains(&e.to))
    //                 .map(|e| (e.to, e.to_position));
    //             next_nodes.extend(to_nodes);
    //         }
    //         next_nodes.sort();
    //         next_nodes.dedup();
    //         if next_nodes.iter().all(|n| !covered.contains(&(n.0))) {
    //             break;
    //         } else {
    //             nodes_at.push(next_nodes);
    //         }
    //     }
    //     nodes_at
    // }
    // /// Remove edge from given reads.
    // pub fn remove_edges_from_short_reads(&mut self, reads: &[&EncodedRead]) {
    //     for read in reads.iter() {
    //         for w in read.nodes.windows(2) {
    //             let from = (w[0].unit, w[0].cluster);
    //             use Position::*;
    //             let from_pos = if w[0].is_forward { Tail } else { Head };
    //             let to = (w[1].unit, w[1].cluster);
    //             let to_pos = if w[1].is_forward { Tail } else { Head };
    //             let mock_edge = DitchEdge::new(from, from_pos, to, to_pos, EdgeLabel::Ovlp(0));
    //             let forward_edge = self
    //                 .nodes
    //                 .get_mut(&from)
    //                 .and_then(|n| n.edges.iter_mut().find(|e| e == &&mock_edge));
    //             if let Some(forward_edge) = forward_edge {
    //                 forward_edge.occ -= 1;
    //             }
    //             let mock_edge = mock_edge.reverse();
    //             let reverse_edge = self
    //                 .nodes
    //                 .get_mut(&to)
    //                 .and_then(|n| n.edges.iter_mut().find(|e| e == &&mock_edge));
    //             if let Some(reverse_edge) = reverse_edge {
    //                 reverse_edge.occ -= 1;
    //             }
    //         }
    //     }
    // }
    /// Remove lightweight edges with occurence less than or equal to `thr`.
    /// To retain that edge if the edge is the only edge from its terminal,
    /// set `retain_single_edge` to `true`.
    pub fn remove_lightweight_edges(&mut self, thr: usize, retain_single_edge: bool) {
        let mut removed_edges: HashMap<_, Vec<_>> =
            self.nodes.keys().map(|&k| (k, vec![])).collect();
        for (from, node) in self.nodes.iter() {
            for &pos in &[Position::Head, Position::Tail] {
                let edges = node.edges.iter().filter(|e| e.from_position == pos);
                if edges.clone().count() <= 1 {
                    continue;
                };
                let thr = {
                    let max_occ = edges.clone().map(|x| x.occ).max().unwrap();
                    if retain_single_edge {
                        thr.min(max_occ - 1)
                    } else {
                        thr
                    }
                };
                for e in edges.clone().filter(|e| e.occ <= thr) {
                    // If retain mode, check this is not the only edge
                    let removable = !retain_single_edge
                        || self.get_edges(e.to, e.to_position).any(|f| e.occ < f.occ);
                    if removable {
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
            let edges: Vec<_> = self.get_edges(key, position).collect();
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
        let edges: Vec<_> = self.get_edges(root, position).collect();
        assert!(edges.len() > 1);
        // TODO: Maybe we need to consensus these bubbles.
        let (first_node, first_pos, seq, edgelabel) = {
            let first = edges.iter().max_by_key(|x| x.occ).unwrap();
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
            while let Some(mut edge) = self.nodes.get_mut(node).unwrap().edges.pop() {
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
        let root_to_edge = new_node
            .edges
            .iter()
            .find(|e| e.to == root && e.to_position == position)
            .unwrap()
            .reverse();
        self.nodes.get_mut(&root).unwrap().edges.push(root_to_edge);
        self.nodes.insert(first_node, new_node);
        let num_edge = self.get_edges(root, position).count();
        assert_eq!(num_edge, 1);
        for (other, other_position) in node_otherside {
            let (removed_occ, label) = self
                .get_edges(other, other_position)
                .filter(|e| merged_nodes.contains(&e.to))
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
    pub fn squish_bubbles(&self, len: usize) -> HashMap<Node, u64> {
        let mut squish_to: HashMap<Node, u64> = HashMap::new();
        for node in self.nodes.values() {
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
                // debug!("BUBBLE\t{:?}\t{}\t{}", node.node, pos, is_bubble);
                if is_bubble {
                    let mut convert_table: HashMap<u64, u64> = HashMap::new();
                    for (path, _) in path_and_dest.iter() {
                        for &(unit, cluster) in path.iter() {
                            convert_table
                                .entry(unit)
                                .and_modify(|x| *x = (*x).min(cluster))
                                .or_insert(cluster);
                        }
                    }
                    // debug!("CONVERT\t{:?}\t{}\t{:?}", node.node, pos, convert_table);
                    for (path, _) in path_and_dest.iter() {
                        for &(unit, cluster) in path.iter() {
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
            if e.from <= e.to {
                (e.from, e.from_position, e.to, e.to_position)
            } else {
                (e.to, e.to_position, e.from, e.from_position)
            }
        };
        for node in self.nodes.values() {
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
        for node in self.nodes.values_mut() {
            node.edges
                .retain(|edge| !removed_edges.contains(&format_edge(edge)));
        }
    }
    /// Check if we can select this edge when we need to select an edge for each
    /// node while preserving the connectiviy of the graph.
    /// true if we could select this edge.
    pub fn can_select(&self, edge: &DitchEdge) -> bool {
        let to = &self.nodes[&edge.to];
        // Filtering out the back edge.
        for to_edge in to
            .edges
            .iter()
            .filter(|to_e| to_e.from_position == edge.to_position && to_e.to != edge.from)
        {
            let sibling = &self.nodes[&to_edge.to];
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

#[cfg(test)]
mod tests {
    use super::*;
    use definitions::EncodedRead;
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
                definitions::Node::new(unit, true, &seq, cigar, position, cl)
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
                let label = String::new();
                definitions::Edge {
                    from: from.unit,
                    to: to.unit,
                    offset: start as i64 - end as i64,
                    label,
                }
            })
            .collect();
        let rem = original_length - original_length / 2_000 * 2_000;
        EncodedRead {
            id,
            original_length,
            leading_gap: vec![],
            trailing_gap: vec![b'A'; rem],
            edges,
            nodes,
        }
    }
    #[test]
    fn from_reads_1() {
        // Generating reads from looping_case
        let node_cp: Vec<_> = vec![2, 2, 8, 2, 2, 4, 4, 2, 2];
        let hap: Vec<_> = vec![0, 1, 2, 4, 3, 2, 6, 5, 2, 6, 5, 2, 7, 8];
        let read_num = 2 * 2_000 * 30 * hap.len() / 10_000;
        let mut rng: Xoroshiro128StarStar = SeedableRng::seed_from_u64(482304);
        let reads: Vec<_> = (0..read_num)
            .map(|i| gen_read(i as u64, &mut rng, &hap))
            .collect();
        let total_units: usize = reads.iter().map(|r| r.nodes.len()).sum();
        let cov = (total_units / hap.len() / 2) as f64;
        let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let assemble_config = AssembleConfig::new(1, 100, false, false, 6, 1f64);
        let graph = DitchGraph::new(&reads, None, ReadType::CCS, &assemble_config);
        let (nodes, _) = graph.copy_number_estimation_gbs(cov, &lens);
        for (i, &cp) in node_cp.iter().enumerate() {
            assert_eq!(cp, nodes[&(i as u64, 0)]);
        }
    }
    #[test]
    fn from_reads_2() {
        let node_cp: Vec<_> = vec![2, 3, 2, 3, 2, 3, 2, 2];
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
        let lens: Vec<_> = reads.iter().map(|r| r.original_length).collect();
        let assemble_config = AssembleConfig::new(1, 100, false, false, 6, 1f64);
        let graph = DitchGraph::new(&reads, None, ReadType::CCS, &assemble_config);
        let (nodes, _) = graph.copy_number_estimation_gbs(cov, &lens);
        for (i, &cp) in node_cp.iter().enumerate() {
            assert_eq!(cp, nodes[&(i as u64, 0)]);
        }
    }
}
