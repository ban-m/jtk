//! module defines how a ditch graph would generate sequence, or reduce simple paths.
use super::*;
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
    ) {
        let mut arrived = HashSet::new();
        let mut sids: HashMap<_, _> = HashMap::new();
        let (mut g_segs, mut g_edges, mut summaries) = (vec![], vec![], vec![]);
        let candidates = self.enumerate_candidates();
        for (node, p) in candidates {
            if arrived.contains(&node) {
                continue;
            }
            let name = format!("tig_{:04}", g_segs.len());
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
            let name = format!("tig_{:04}", g_segs.len());
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
        let uid = Some(format!("group-{}", 0));
        let group = gfa::Group::Set(gfa::UnorderedGroup { uid, ids });
        (g_segs, g_edges, group, summaries)
    }
    fn enumerate_adjacent_tag(
        &self,
        seqname: &str,
        start: Node,
        start_position: Position,
        sids: &HashMap<Node, ContigTag>,
        gfa_edge_start: gfa::Position,
    ) -> Vec<(gfa::Edge, Vec<gfa::SamTag>)> {
        self.get_edges(start, start_position)
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
                // let beg1 = gfa::Position::from(0, false);
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
    // let tail_edges = self.nodes[&node]
    // .edges
    // .iter()
    // .filter(|e| e.from_position == position)
    // .filter_map(|e| {
    //     let (sid2, beg2) = match sids.get(&e.to)? {
    //         &ContigTag::Start(ref name, pos, _) if pos == e.to_position => {
    //             (gfa::RefID::from(name, true), gfa::Position::from(0, false))
    //         }
    //         &ContigTag::End(ref name, pos, len) if pos == e.to_position => {
    //             (gfa::RefID::from(name, true), gfa::Position::from(len, true))
    //         }
    //         &ContigTag::Both(ref name, pos, _, _) if pos == e.to_position => {
    //             (gfa::RefID::from(name, true), gfa::Position::from(0, false))
    //         }
    //         &ContigTag::Both(ref name, _, pos, len) if pos == e.to_position => {
    //             (gfa::RefID::from(name, true), gfa::Position::from(len, true))
    //         }
    //         _ => return None,
    //     };
    //     let eid = None;
    //     let sid1 = gfa::RefID::from(&seqname, true);
    //     let beg1 = gfa::Position::from(seq.len(), true);
    //     let a = None;
    //     let edge = gfa::Edge::from(eid, sid1, sid2, beg1, beg1, beg2, beg2, a);
    //     let mut samtag = vec![
    //         gfa::SamTag::new(format!("cv:i:{}", e.occ)),
    //         gfa::SamTag::new(format!("ln:i:{}", e.seq.len())),
    //     ];
    //     if let Some(cp) = e.copy_number {
    //         samtag.push(gfa::SamTag::new(format!("cp:i:{}", cp)));
    //     }
    //     Some((edge, samtag))
    // });

    // Traverse from the given `start` node of `start_position` Position.
    // The retuened ContigSummary contains which (unit,cluster) elements is used to
    // construct this contig, in the order of apprearance.
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
        let gfa_pos = gfa::Position::from(0, false);
        let edges = self.enumerate_adjacent_tag(&seqname, start, start_position, sids, gfa_pos);
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
                    copy_number: self.nodes[&node].copy_number,
                };
                unit_names.push(elm);
            }
            seq += &cons;
            position = !position;
            // Check.
            let num_edges = self.get_edges(node, position).count();
            if num_edges == 0 || num_edges > 1 {
                break;
            }
            // There is only one child.
            let selected_edge = self.get_edges(node, position).next().unwrap();
            assert_eq!(selected_edge.from, node);
            assert_eq!(selected_edge.from_position, position);
            // Succeed along with the edge.
            match &selected_edge.seq {
                EdgeLabel::Ovlp(l) => assert!((0..(-l)).all(|_| seq.pop().is_some())),
                EdgeLabel::Seq(label) => seq.extend(label.iter().map(|&x| x as char)),
            };
            // Append nodes proxied by this edge to the last node. The order is correct(I created so.).
            for &((unit, cluster), strand, occ) in selected_edge.proxying.iter() {
                let elm = ContigElement {
                    unit,
                    cluster,
                    strand,
                    occ,
                    copy_number: None,
                };
                unit_names.push(elm);
            }
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
        let gfa_pos = gfa::Position::from(seq.len(), true);
        let tail_edges = self.enumerate_adjacent_tag(&seqname, node, position, sids, gfa_pos);
        let mut edges = edges;
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
