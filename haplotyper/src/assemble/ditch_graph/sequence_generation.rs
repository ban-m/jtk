//! module defines how a ditch graph would generate sequence, or reduce simple paths.
//! 2022/07/07: I changed the strategy.
use super::*;

/// A summary of a contig. It tells us
/// the unit, the cluster, and the direction
/// it spelled.
#[derive(Debug, Clone)]
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

#[derive(Debug, Clone)]
pub struct ContigElement {
    pub unit: u64,
    pub cluster: u64,
    pub strand: bool,
    /// The "coverage."
    pub occ: usize,
    pub copy_number: Option<usize>,
}

impl ContigElement {
    fn new(node: &DitchNode, position: Position) -> Self {
        let (unit, cluster) = node.node;
        let strand = position == Position::Head;
        let occ = node.occ;
        let copy_number = node.copy_number;
        Self {
            unit,
            cluster,
            strand,
            occ,
            copy_number,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ContigEncoding {
    pub id: String,
    tiles: Vec<UnitAlignmentInfo>,
    // Is this OK?
    cache: HashMap<(Node, u8), usize>,
}

impl ContigEncoding {
    pub fn tiles(&self) -> &[UnitAlignmentInfo] {
        self.tiles.as_slice()
    }
    pub fn matches<'a>(&'a self, node: Node, direction: bool) -> MatchIter<'a> {
        MatchIter {
            inner: self,
            node,
            direction,
            hit: 0,
        }
    }
    fn new(id: &str) -> Self {
        Self {
            id: id.to_string(),
            tiles: Vec::new(),
            cache: HashMap::new(),
        }
    }
    fn push(&mut self, tile: UnitAlignmentInfo) {
        let k = tile.unit_info();
        let hit = (0..).find(|&i| !self.cache.contains_key(&(k, i))).unwrap();
        self.cache.insert((k, hit), self.tiles.len());
        self.tiles.push(tile);
    }
}

#[derive(Debug, Clone)]
pub struct MatchIter<'a> {
    inner: &'a ContigEncoding,
    node: Node,
    direction: bool,
    hit: u8,
}

impl<'a> std::iter::Iterator for MatchIter<'a> {
    type Item = (usize, &'a UnitAlignmentInfo);
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.inner.cache.get(&(self.node, self.hit)) {
                Some(&idx) if self.inner.tiles[idx].unit_direction == self.direction => break,
                None => return None,
                Some(_) => self.hit += 1,
            }
        }
        self.hit += 1;
        let idx = *self.inner.cache.get(&(self.node, self.hit - 1)).unwrap();
        self.inner.tiles.get(idx).map(|t| (idx, t))
    }
}

#[derive(Debug, Clone, Default)]
pub struct UnitAlignmentInfo {
    unit: u64,
    cluster: u64,
    contig_start: usize,
    contig_end: usize,
    // If false, it should be rev-comped.
    unit_direction: bool,
    unit_start: usize,
    unit_end: usize,
    unit_len: usize,
}

impl UnitAlignmentInfo {
    pub fn unit_len(&self) -> usize {
        self.unit_len
    }
    pub fn contig_range(&self) -> (usize, usize) {
        (self.contig_start, self.contig_end)
    }
    pub fn unit_range(&self) -> (bool, usize, usize) {
        (self.unit_direction, self.unit_start, self.unit_end)
    }
    pub fn unit_info(&self) -> Node {
        (self.unit, self.cluster)
    }
    pub fn unit_and_dir_info(&self) -> (Node, bool) {
        ((self.unit, self.cluster), self.unit_direction)
    }
    pub fn new(
        (unit, cluster): (u64, u64),
        (unit_direction, unit_start, unit_end): (bool, usize, usize),
        (contig_start, contig_end): (usize, usize),
        unit_len: usize,
    ) -> Self {
        Self {
            unit,
            cluster,
            contig_start,
            contig_end,
            unit_direction,
            unit_start,
            unit_end,
            unit_len,
        }
    }
}

#[derive(Debug, Copy, Clone, Default)]
struct UnitAlnInfoBuilder {
    pub unit: Option<u64>,
    pub cluster: Option<u64>,
    pub contig_start: Option<usize>,
    pub contig_end: Option<usize>,
    /// If false, it should be rev-comped.
    pub unit_direction: Option<bool>,
    pub unit_start: Option<usize>,
    pub unit_end: Option<usize>,
    pub unit_len: Option<usize>,
}

impl UnitAlnInfoBuilder {
    fn set_node(mut self, (unit, cluster): Node) -> Self {
        self.unit = Some(unit);
        self.cluster = Some(cluster);
        self
    }
    fn set_contig_start(mut self, start: usize) -> Self {
        self.contig_start = Some(start);
        self
    }
    fn set_contig_end(mut self, end: usize) -> Self {
        self.contig_end = Some(end);
        self
    }
    fn set_unit_direction(mut self, position: Position) -> Self {
        self.unit_direction = Some(position == Position::Head);
        self
    }
    fn set_unit_start(mut self, start: usize) -> Self {
        self.unit_start = Some(start);
        self
    }
    fn set_unit_end(mut self, end: usize) -> Self {
        self.unit_end = Some(end);
        self
    }
    fn set_unit_len(mut self, len: usize) -> Self {
        self.unit_len = Some(len);
        self
    }
    fn build(&self) -> UnitAlignmentInfo {
        let Self {
            unit,
            cluster,
            contig_start,
            contig_end,
            unit_direction,
            unit_start,
            unit_end,
            unit_len,
        } = *self;
        UnitAlignmentInfo {
            unit: unit.unwrap(),
            cluster: cluster.unwrap(),
            contig_start: contig_start.unwrap(),
            contig_end: contig_end.unwrap(),
            unit_direction: unit_direction.unwrap(),
            unit_start: unit_start.unwrap(),
            unit_end: unit_end.unwrap(),
            unit_len: unit_len.unwrap(),
        }
    }
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

impl<'a> super::DitchGraph<'a> {
    /// Reduce simple path of this graph and returns the edges and nodes of the reduced graph..
    /// The contig summary contains the units used to construct a contig in the traversing order.
    #[allow(clippy::type_complexity)]
    pub fn spell(
        &self,
        c: &AssembleConfig,
    ) -> (
        Vec<gfa::Segment>,
        Vec<(gfa::Edge, Vec<gfa::SamTag>)>,
        gfa::Group,
        Vec<ContigSummary>,
        Vec<ContigEncoding>,
    ) {
        let mut arrived = HashSet::new();
        let mut sids: HashMap<_, _> = HashMap::new();
        let (mut g_segs, mut g_edges, mut summaries) = (vec![], vec![], vec![]);
        let mut unit_position_on_contigs = vec![];
        let mut candidates = self.enumerate_candidates();
        candidates.sort();
        for (node, p) in candidates {
            if arrived.contains(&node) {
                continue;
            }
            let name = format!("tig_{:04}", g_segs.len());
            let contig_info = self.traverse_from(&mut arrived, &mut sids, node, p, name, c);
            let (contig, edges, summary, unit_positions) = contig_info;
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
            unit_position_on_contigs.push(unit_positions);
        }
        let mut nodes: Vec<_> = self.nodes().map(|x| x.0).collect();
        nodes.sort();
        for key in nodes {
            if arrived.contains(&key) {
                continue;
            }
            let p = Position::Head;
            let name = format!("tig_{:04}", g_segs.len());
            let contig_info = self.traverse_from(&mut arrived, &mut sids, key, p, name, c);
            let (contig, edges, summary, unit_positions) = contig_info;
            g_segs.push(contig);
            g_edges.extend(edges);
            summaries.push(summary);
            unit_position_on_contigs.push(unit_positions);
        }
        for (idx, _) in self.nodes() {
            assert!(arrived.contains(&idx));
        }
        let ids: Vec<_> = g_segs
            .iter()
            .map(|seg| seg.sid.clone())
            .chain(g_edges.iter().filter_map(|e| e.0.eid.clone()))
            .collect();
        let uid = Some(format!("group-{}", 0));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group, summaries, unit_position_on_contigs)
    }
    fn enumerate_adjacent_tag(
        &self,
        seqname: &str,
        start: NodeIndex,
        start_position: Position,
        sids: &HashMap<NodeIndex, ContigTag>,
        gfa_edge_start: gfa::Position,
    ) -> Vec<(gfa::Edge, Vec<gfa::SamTag>)> {
        let edges = self.edges_from(start, start_position);
        edges
            .into_iter()
            .filter_map(|e| {
                let (sid2, beg2) = match *sids.get(&e.to)? {
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
                let sid1 = gfa::RefID::from(seqname, true);
                let beg1 = gfa_edge_start;
                let a = None;
                let edge = gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a);
                let mut samtag = vec![
                    gfa::SamTag::new(format!("cv:i:{}", e.occ)),
                    gfa::SamTag::new(format!("ln:i:{}", e.seq.len())),
                ];
                if let Some(cp) = e.copy_number {
                    samtag.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
                }
                Some((edge, samtag))
            })
            .collect()
    }

    // Traverse from the given `start` node of `start_position` Position.
    // The retuened ContigSummary contains which (unit,cluster) elements is used to
    // construct this contig, in the order of apprearance.
    fn traverse_from(
        &self,
        arrived: &mut HashSet<NodeIndex>,
        sids: &mut HashMap<NodeIndex, ContigTag>,
        start: NodeIndex,
        start_position: Position,
        seqname: String,
        _c: &AssembleConfig,
    ) -> (
        gfa::Segment,
        Vec<(gfa::Edge, Vec<gfa::SamTag>)>,
        ContigSummary,
        ContigEncoding,
    ) {
        // Find edges.
        // I impled here!
        let gfa_pos = gfa::Position::from(0, false);
        let edges = self.enumerate_adjacent_tag(&seqname, start, start_position, sids, gfa_pos);
        let mut position_of_units = ContigEncoding::new(&seqname);
        let (mut node_index, mut position) = (start, start_position);
        let mut seq = self.initial_sequence(start, start_position);
        // Start traveresing.
        let mut unit_names = vec![];
        loop {
            let mut builder = UnitAlnInfoBuilder::default();
            let node = self.node(node_index).unwrap();
            builder = builder
                .set_node(node.node)
                .set_contig_start(seq.len())
                .set_unit_start(0)
                .set_unit_direction(position)
                .set_unit_len(node.seq().len());
            arrived.insert(node_index);
            // Move forward.
            let cons = match position {
                Position::Head => node.seq_as_string(),
                Position::Tail => revcmp_str(&node.seq_as_string()),
            };
            unit_names.push(ContigElement::new(node, position));
            seq += &cons;
            position = !position;
            // Check.
            if self.count_edges(node_index, position) != 1 {
                builder = builder
                    .set_contig_end(seq.len())
                    .set_unit_end(node.seq().len());
                position_of_units.push(builder.build());
                break;
            }
            // There is only one child.
            let selected_edge = &node
                .edges
                .iter()
                .find(|e| e.from_position == position)
                .unwrap();
            assert_eq!(selected_edge.from, node_index);
            assert_eq!(selected_edge.from_position, position);
            // Succeed along with the edge.
            match &selected_edge.seq {
                EdgeLabel::Ovlp(l) => {
                    let _ = (0..(-l)).filter_map(|_| seq.pop()).count();
                    builder = builder
                        .set_contig_end(seq.len())
                        .set_unit_end(node.seq().len() - (-l) as usize);
                }
                EdgeLabel::Seq(label) => {
                    builder = builder
                        .set_contig_end(seq.len())
                        .set_unit_end(node.seq().len());
                    seq.extend(label.iter().map(|&x| x as char));
                }
            };
            position_of_units.push(builder.build());
            let (next, next_position) = (selected_edge.to, selected_edge.to_position);
            // Check the number of child.
            let num_children = self.count_edges(next, next_position);
            // Or looping...
            if num_children >= 2 || arrived.contains(&next) {
                break;
            }
            // Jump to the next node.
            assert!(num_children == 1);
            node_index = next;
            position = next_position;
        }
        seq += &self.trailing_sequence(node_index, position);
        let summary = ContigSummary::new(&seqname, &unit_names);
        // Register start and tail node.
        // This if statement is no needed?
        if start == node_index {
            let tag = ContigTag::Both(seqname.clone(), start_position, !start_position, seq.len());
            sids.insert(start, tag);
        } else {
            let tag = ContigTag::Start(seqname.clone(), start_position, seq.len());
            sids.insert(start, tag);
            let tag = ContigTag::End(seqname.clone(), position, seq.len());
            sids.insert(node_index, tag);
        }
        let gfa_pos = gfa::Position::from(seq.len(), true);
        let seg = gfa::Segment::from(seqname.clone(), seq.len(), Some(seq));
        let tail_edges = self.enumerate_adjacent_tag(&seqname, node_index, position, sids, gfa_pos);
        let mut edges = edges;
        edges.extend(tail_edges);
        (seg, edges, summary, position_of_units)
    }
    fn initial_sequence(&self, node: NodeIndex, position: Position) -> String {
        if self.count_edges(node, position) > 0 {
            String::new()
        } else {
            self.node(node)
                .unwrap()
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
    fn trailing_sequence(&self, node: NodeIndex, position: Position) -> String {
        if self.count_edges(node, position) > 0 {
            String::new()
        } else {
            self.node(node)
                .unwrap()
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
