use definitions::*;
use rayon::prelude::*;
const IDEN: f64 = 0.3;
const OVLP: i32 = 3;
const MISMATCH: f64 = 0.15;
const MATCH: f64 = 1f64 - MISMATCH;
const MISM_SCORE: f64 = -1.897;
const GAP_SCORE: f64 = -4.605;
const MISM_UNIT: f64 = -10000000f64;
// TODO: Convert them into Vector type.
use crate::find_union::FindUnion;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct AssembleConfig {
    identity: f64,
    ovlp: i32,
}

impl AssembleConfig {
    pub fn new(ovlp: i32, identity: f64) -> Self {
        Self { identity, ovlp }
    }
}

pub const DEFAULT_CONFIG: AssembleConfig = AssembleConfig {
    identity: IDEN,
    ovlp: OVLP,
};

#[derive(Clone)]
pub struct StringGraph<'a> {
    // TODO: nodes should be (&Node, bool) to indicate the
    // direction of the reads. Maybe it fasten the alignment.
    pub edges: HashMap<u64, Vec<StrEdge>>,
    pub nodes: HashMap<u64, Vec<(&'a Node, bool)>>,
    pub reads: Vec<&'a EncodedRead>,
    pub parent: FindUnion,
}

impl<'a> std::fmt::Debug for StringGraph<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let reads = self.reads.len();
        let nodes = self.nodes.len();
        let edges = self.edges.values().map(|x| x.len()).sum::<usize>();
        write!(f, "R:{}\tN:{}\tE:{}", reads, nodes, edges)
    }
}

impl<'a> StringGraph<'a> {
    pub fn sanity_check(&self) -> bool {
        self.edges.values().all(|edges| {
            edges.iter().all(|edge| {
                self.edges[&edge.to].iter().any(|f| {
                    f.to == edge.from
                        && f.to_position == edge.from_position
                        && f.from == edge.to
                        && f.from_position == edge.to_position
                })
            })
        })
    }
    pub fn from_dataset(ds: &'a DataSet, c: &AssembleConfig) -> Self {
        // use crate::global_clustering::{
        //     error_correction::local_correction, GlobalClusteringConfig,
        // };
        // let config = GlobalClusteringConfig::new(3, 2, 1, -1, -1);
        // let reads = local_correction(ds, &config);
        // let assignments: HashMap<_, _> = ds.assignments.iter().map(|r| (r.id, r.cluster)).collect();
        let scoring_scheme = get_match_score(ds);
        let edges: HashMap<_, Result<Vec<_>, _>> = ds
            .encoded_reads
            .par_iter()
            .map(|r| {
                let mut edges = vec![];
                for q in ds.encoded_reads.iter().filter(|q| q.id != r.id)
                //  && assignments.get(&q.id) == assignments.get(&r.id))
                {
                    // match StrEdge::from_corrected_reads(r, q, c, &scoring_scheme) {
                    match StrEdge::from_reads(r, q, c, &scoring_scheme) {
                        Ok(res) => edges.push(res),
                        Err(is_contained) if is_contained => return (r.id, Err(q.id)),
                        _ => {}
                    }
                }
                (r.id, Ok(edges))
            })
            .collect();
        let reads: Vec<_> = ds.encoded_reads.iter().collect();
        Self::from_edges(edges, &reads)
    }
    pub fn new(reads: &'a [EncodedRead], c: &AssembleConfig) -> Self {
        let scoring: HashMap<_, _> = reads
            .iter()
            .flat_map(|read| read.nodes.iter())
            .map(|n| ((n.unit, n.cluster), 4f64.recip().ln()))
            .collect();
        let edges: HashMap<_, Result<Vec<_>, _>> = reads
            .par_iter()
            .map(|r| {
                let mut edges = vec![];
                for q in reads.iter().filter(|q| q.id != r.id) {
                    match StrEdge::from_reads(r, q, c, &scoring) {
                        Ok(res) => edges.push(res),
                        Err(is_contained) if is_contained => return (r.id, Err(q.id)),
                        _ => {}
                    }
                }
                (r.id, Ok(edges))
            })
            .collect();
        let reads: Vec<_> = reads.iter().collect();
        Self::from_edges(edges, &reads)
    }
    fn from_edges(
        edges: HashMap<u64, Result<Vec<StrEdge>, u64>>,
        reads: &[&'a EncodedRead],
    ) -> Self {
        let max_id = reads.iter().map(|r| r.id).max().unwrap();
        let mut fu = FindUnion::new((max_id + 1) as usize);
        {
            // DUMPS
            let total = reads.len();
            let available = edges
                .values()
                .filter(|x| match x {
                    Ok(res) => !res.is_empty(),
                    _ => false,
                })
                .count();
            let orphans = edges
                .values()
                .filter(|x| match x {
                    Ok(res) => res.is_empty(),
                    _ => false,
                })
                .count();
            let contained = edges.values().filter(|x| matches!(x, Err(_))).count();
            debug!("NODESUM\tTotal\tValid\tOrphan\tContained");
            debug!(
                "NODESUM\t{}\t{}\t{}\t{}",
                total, available, orphans, contained
            );
        }
        // Contained read removal.
        for (&from, edge) in edges.iter() {
            if let Err(to) = edge {
                fu.unite(from as usize, *to as usize).unwrap();
            }
        }
        let is_contained: HashSet<u64> = edges
            .iter()
            .filter_map(|(&k, v)| if v.is_err() { Some(k) } else { None })
            .collect();
        let isolated: HashSet<u64> = edges
            .values()
            .filter_map(|edge| edge.as_ref().err())
            .filter(|node| match edges[&node] {
                Ok(ref res) => res.is_empty(),
                _ => false,
            })
            .copied()
            .collect();
        debug!("{} Alignments..", edges.len());
        let edges: HashMap<u64, Vec<_>> = edges
            .into_iter()
            .filter_map(|(k, v)| {
                v.ok().and_then(|mut es| {
                    if es.iter().all(|e| is_contained.contains(&e.to)) {
                        for to in es.iter().map(|e| e.to) {
                            fu.unite(k as usize, to as usize).unwrap();
                        }
                    }
                    es.retain(|e| !is_contained.contains(&e.to));
                    if !es.is_empty() || isolated.contains(&k) {
                        assert!(es.iter().all(|e| e.from == k));
                        es.sort_by_key(|e| e.to);
                        Some((k, es))
                    } else {
                        None
                    }
                })
            })
            .collect();
        debug!(
            "Into {} alignments.",
            edges.values().map(|x| x.len()).sum::<usize>()
        );
        let is_used: HashSet<_> = edges.keys().copied().collect();
        let nodes: HashMap<u64, Vec<_>> = reads
            .iter()
            .filter(|r| is_used.contains(&r.id))
            .map(|r| (r.id, r.nodes.iter().map(|n| (n, true)).collect()))
            .collect();
        Self {
            reads: reads.to_vec(),
            edges,
            nodes,
            parent: fu,
        }
    }
    pub fn transitive_edge_reduction(&mut self) {
        assert!(self
            .edges
            .values()
            .all(|edges| edges.iter().map(|e| e.to).is_sorted()));
        let is_transitive_edge: HashMap<u64, Vec<bool>> = self
            .edges
            .iter()
            .map(|(&k, edges)| {
                let is_transitive: Vec<_> = edges
                    .iter()
                    .map(|edge| {
                        // TODO: More elegant code.
                        let mut f_ptr = edges.iter().peekable();
                        // At least there's one edge, the reverse edges of `edge` itself.
                        let mut t_ptr = self.edges[&edge.to].iter().peekable();
                        while let (Some(f_edge), Some(t_edge)) = (f_ptr.peek(), t_ptr.peek()) {
                            if f_edge.to < t_edge.to {
                                f_ptr.next();
                            } else if t_edge.to < f_edge.to {
                                t_ptr.next();
                            // this edge should not be the reverse edge of `edge`.
                            } else if f_edge.to == edge.to || t_edge.to == edge.from {
                                f_ptr.next();
                                t_ptr.next();
                            } else if edge.is_transitive(f_edge, t_edge) {
                                return true;
                            } else {
                                f_ptr.next();
                                t_ptr.next();
                            }
                        }
                        false
                    })
                    .collect();
                (k, is_transitive)
            })
            .collect();
        // Remove transitive edges.
        let num_edges = self.edges.values().map(|x| x.len()).sum::<usize>();
        self.edges.iter_mut().for_each(|(k, v)| {
            let transitive_flags = &is_transitive_edge[&k];
            let mut idx = 0;
            v.retain(|_| {
                idx += 1;
                !transitive_flags[idx - 1]
            });
        });
        let after_edges = self.edges.values().map(|x| x.len()).sum::<usize>();
        debug!("{}=>{}", num_edges, after_edges);
        for key in self.nodes.keys() {
            assert!(self.edges.contains_key(key));
        }
        // Sanity check.
        for edges in self.edges.values() {
            for edge in edges.iter() {
                assert!(self.edges[&edge.to].iter().any(|f| {
                    f.from == edge.to
                        && f.from_position == edge.to_position
                        && f.to == edge.from
                        && f.to_position == edge.from_position
                }));
            }
        }
    }
    pub fn simple_path_reduction(&mut self) {
        let ids: Vec<_> = self.nodes.keys().copied().collect();
        for id in ids {
            if self.edges.contains_key(&id) {
                self.reduce_edges_from(id).unwrap();
            }
        }
    }
    //Remove dead end.
    pub fn tip_removal(&mut self) {
        let mut to_be_removed: HashSet<_> = HashSet::new();
        for (_, edges) in self.edges.iter() {
            for position in vec![Position::Head, Position::Tail] {
                let num_edges = edges
                    .iter()
                    .filter(|ed| ed.from_position == position)
                    .count();
                if num_edges < 2 {
                    continue;
                }
                for ed in edges.iter().filter(|ed| ed.from_position == position) {
                    if self.edges[&ed.to].len() == 1 {
                        let tip_edge = self.edges[&ed.to][0];
                        assert_eq!(ed.from, tip_edge.to);
                        assert_eq!(ed.to, tip_edge.from);
                        assert_eq!(ed.from_position, tip_edge.to_position);
                        assert_eq!(ed.to_position, tip_edge.from_position);
                        // This node is to be removed.
                        to_be_removed.insert(ed.to);
                    }
                }
            }
        }
        for remove_node in to_be_removed.iter() {
            self.nodes.remove(remove_node);
            self.edges.remove(remove_node);
        }
        for edges in self.edges.values_mut() {
            edges.retain(|edge| !to_be_removed.contains(&edge.to));
        }
    }
    // TODO:We need DFS approach.
    fn reduce_edges_from(&mut self, id: u64) -> Option<()> {
        // True if some edge was marged.
        let mut marged = true;
        while marged {
            marged = false;
            for &pos in &[Position::Head, Position::Tail] {
                if self.edges[&id]
                    .iter()
                    .filter(|e| e.from_position == pos)
                    .count()
                    != 1
                {
                    continue;
                }
                // Never panic.
                let &edge = self.edges[&id].iter().find(|e| e.from_position == pos)?;
                let num_edges = self.edges[&edge.to]
                    .iter()
                    .filter(|e| e.from_position == edge.to_position)
                    .count();
                if num_edges != 1 || edge.to == edge.from {
                    continue;
                }
                marged = true;
                // Never panic.
                self.edges.get_mut(&id)?.retain(|e| e != &edge);
                self.reduce_edge(edge)?;
            }
        }
        Some(())
    }
    fn reduce_edge(&mut self, edge: StrEdge) -> Option<()> {
        let mut poped_seq = self.nodes.remove(&edge.to)?;
        let edge_seq = self.nodes.get_mut(&edge.from)?;
        use Position::*;
        // Elong segment.
        match (edge.from_position, edge.to_position) {
            (Tail, Head) => {
                for _ in 0..edge.offset {
                    edge_seq.pop();
                }
                edge_seq.extend(poped_seq);
            }
            (Tail, Tail) => {
                for _ in 0..edge.offset {
                    edge_seq.pop();
                }
                let push_seq = poped_seq.into_iter().rev().map(|(n, b)| (n, !b));
                edge_seq.extend(push_seq);
            }
            (Head, Tail) => {
                poped_seq.extend(edge_seq.iter().copied().skip(edge.offset));
                edge_seq.clear();
                edge_seq.extend(poped_seq);
            }
            (Head, Head) => {
                let mut ps = poped_seq;
                // Revcmp
                ps.reverse();
                ps.iter_mut().for_each(|x| x.1 = !x.1);
                ps.extend(edge_seq.iter().copied().skip(edge.offset));
                edge_seq.clear();
                edge_seq.extend(ps);
            }
        }
        // Remove reverse edge, append them to `edge.from` node..
        let mut new_edges = self.edges.remove(&edge.to)?;
        new_edges.retain(|f| {
            !(f.to == edge.from
                && f.to_position == edge.from_position
                && f.from == edge.to
                && f.from_position == edge.to_position)
        });
        new_edges.iter_mut().for_each(|e| {
            self.edges
                .get_mut(&e.to)
                .unwrap()
                .iter_mut()
                .filter(|f| f.to == e.from && f.to_position == e.from_position)
                .for_each(|f| {
                    f.to = edge.from;
                    f.to_position = edge.from_position;
                });
            e.from = edge.from;
            e.from_position = edge.from_position;
        });
        self.edges.get_mut(&edge.from)?.extend(new_edges);
        // Sort the edges.
        self.edges.get_mut(&edge.from)?.sort_by_key(|x| x.to);
        assert!(self.edges.get(&edge.from)?.is_sorted_by_key(|x| x.to));
        if log::log_enabled!(log::Level::Debug) {
            for edges in self.edges.values() {
                for edge in edges {
                    match self.edges.get(&edge.to) {
                        Some(edges) => {
                            let is_ok = edges.iter().any(|f| {
                                f.from == edge.to
                                    && f.from_position == edge.to_position
                                    && f.to == edge.from
                                    && f.to_position == edge.from_position
                            });
                            assert!(is_ok, "{:?}->{:?}", edge, edges);
                        }
                        None => panic!("{:?}", edge),
                    }
                }
            }
        }
        self.parent.unite(edge.from as usize, edge.to as usize)
    }
    pub fn light_node_trim(&mut self, thr: usize) {
        let nodes: Vec<_> = self.nodes.keys().copied().collect();
        let ids: Vec<_> = self.reads.iter().map(|r| r.id).collect();
        let is_lightweight: HashSet<u64> = nodes
            .into_iter()
            .filter(|&k| {
                let parent = self.parent.find(k as usize).unwrap();
                let coverage = ids
                    .iter()
                    .filter(|&&id| self.parent.find(id as usize).unwrap() == parent)
                    .count();
                coverage < thr
            })
            .collect();
        for edges in self.edges.values_mut() {
            edges.retain(|edge| {
                !is_lightweight.contains(&edge.from) && !is_lightweight.contains(&edge.to)
            });
        }
    }
    pub fn node_count(&mut self) {
        let mut cluster: HashMap<_, _> = HashMap::new();
        let mut cl = 0;
        let ids: Vec<_> = self.reads.iter().map(|r| r.id as usize).collect();
        for &id in ids.iter() {
            if self.parent.find(id).unwrap() == id {
                let members: Vec<_> = ids
                    .iter()
                    .filter(|&q| self.parent.find(*q).unwrap() == id)
                    .copied()
                    .collect();
                cluster.insert(cl, members);
                cl += 1;
            }
        }
        let nodes: Vec<_> = self.nodes.keys().copied().collect();
        let tot_reads = nodes
            .iter()
            .map(|&id| {
                cluster
                    .values()
                    .find(|mem| mem.contains(&(id as usize)))
                    .unwrap()
                    .len()
            })
            .sum::<usize>();
        let sum_cl = cluster.values().map(|xs| xs.len()).sum::<usize>();
        debug!("NODES\t{}\t{}\t{}", self.reads.len(), tot_reads, sum_cl);
    }
    pub fn assemble_as_gfa(&mut self) -> gfa::GFA {
        use gfa::*;
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let (nodes, edges, group, summaries) = self.assemble();
        let nodes_length: HashMap<_, _> = nodes.iter().map(|n| (n.sid.clone(), n.slen)).collect();
        let mut removed_nodes = HashSet::new();
        let group: Vec<_> = group
            .into_iter()
            .filter(|g| {
                let short_contig = g
                    .set()
                    .unwrap()
                    .iter()
                    .all(|id| nodes_length.get(id).unwrap() < &100_000);
                let singleton = g.len() == 1;
                if short_contig && singleton {
                    for id in g.set().unwrap().iter() {
                        removed_nodes.insert(id.clone());
                    }
                    false
                } else {
                    true
                }
            })
            .collect();
        let nodes: Vec<_> = nodes
            .into_iter()
            .filter(|node| !removed_nodes.contains(&node.sid))
            .collect();
        let total = removed_nodes.iter().map(|k| nodes_length[k]).sum::<u64>();
        debug!("Removed {} contigs({}bp)", removed_nodes.len(), total);
        // Edge would be the same, unless it is circular loop.
        let edges: Vec<_> = edges
            .into_iter()
            .filter(|edge| {
                !removed_nodes.contains(&edge.sid1.id) && !removed_nodes.contains(&edge.sid2.id)
            })
            .collect();
        for summary in summaries.iter().filter(|n| !removed_nodes.contains(&n.id)) {
            debug!("CONTIG\t{}\t{}", summary.id, summary.reads.len());
        }
        let mut records = vec![];
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
        let group = group
            .into_iter()
            .map(gfa::Content::Group)
            .map(|g| gfa::Record::from_contents(g, vec![]));
        records.extend(group);
        let mut header = vec![header];
        header.extend(records);
        GFA::from_records(header)
    }
    fn node_consensus(&self) -> HashMap<(u64, u64), Vec<u8>> {
        let mut slots: HashMap<_, Vec<_>> = HashMap::new();
        for read in self.reads.iter() {
            for node in read.nodes.iter() {
                slots
                    .entry((node.unit, node.cluster))
                    .or_default()
                    .push(node);
            }
        }
        slots
            .into_par_iter()
            .map(|(key, val)| (key, Self::consensus_from_nodes(&val)))
            .collect()
    }
    fn consensus_from_nodes(nodes: &[&Node]) -> Vec<u8> {
        use crate::assemble::ditch_graph::consensus;
        use crate::local_clustering::node_to_subchunks;
        let chunks: Vec<_> = nodes.iter().map(|n| node_to_subchunks(n, 100)).collect();
        let max_pos = chunks
            .iter()
            .filter_map(|cs| cs.iter().map(|c| c.pos).max())
            .max()
            .unwrap();
        let mut subseqs = vec![vec![]; max_pos + 1];
        for cs in chunks.iter() {
            for c in cs.iter() {
                subseqs[c.pos].push(c.seq.as_slice());
            }
        }
        subseqs.iter().flat_map(|seq| consensus(seq, 10)).collect()
    }
    fn edge_consensus(&self) -> HashMap<((u64, u64), (u64, u64)), EdgeSeq> {
        let mut slots: HashMap<_, Vec<_>> = HashMap::new();
        for read in self.reads.iter() {
            for (w, e) in read.nodes.windows(2).zip(read.edges.iter()) {
                if w[0].unit < w[1].unit {
                    let frm = (w[0].unit, w[0].cluster);
                    let to = (w[1].unit, w[1].cluster);
                    let seq = e.label.as_bytes().to_vec();
                    slots.entry((frm, to)).or_default().push((e.offset, seq));
                } else {
                    let frm = (w[1].unit, w[1].cluster);
                    let to = (w[0].unit, w[0].cluster);
                    let seq = bio_utils::revcmp(e.label.as_bytes());
                    slots.entry((frm, to)).or_default().push((e.offset, seq));
                };
            }
        }
        slots
            .into_iter()
            .map(|(key, val)| {
                let mean = val.iter().map(|x| x.0).sum::<i64>() / val.len() as i64;
                let cons = if mean <= 0 {
                    EdgeSeq::Offset(-mean)
                } else {
                    let longest = val
                        .into_iter()
                        .map(|x| x.1)
                        .max_by_key(|x| x.len())
                        .unwrap();
                    EdgeSeq::Seq(longest)
                };
                (key, cons)
            })
            .collect()
    }
    pub fn assemble(
        &mut self,
    ) -> (
        Vec<gfa::Segment>,
        Vec<gfa::Edge>,
        Vec<gfa::Group>,
        Vec<ContigSummary>,
    ) {
        // Simple path reduced.
        // We map each node to gfa::Segment,
        // each edge to gfa::Edge,
        // and each connected component into gfa::Group.
        let (group, cluster_num) = {
            let mut current_cluster = 0;
            let mut group: HashMap<_, _> = HashMap::new();
            for &id in self.nodes.keys() {
                if group.contains_key(&id) {
                    continue;
                }
                let mut stack = vec![id];
                'outer: while !stack.is_empty() {
                    let last = stack.last().unwrap();
                    if !group.contains_key(last) {
                        group.insert(*last, current_cluster);
                    }
                    for child in self.edges[last].iter().map(|e| e.to) {
                        if !group.contains_key(&child) {
                            stack.push(child);
                            continue 'outer;
                        }
                    }
                    stack.pop();
                }
                current_cluster += 1;
            }
            (group, current_cluster)
        };
        let now = std::time::Instant::now();
        debug!("Taking consensus(node)...");
        let node_sequences = self.node_consensus();
        debug!("{:?}", std::time::Instant::now() - now);
        let now = std::time::Instant::now();
        debug!("Taking consensus(edge)..");
        debug!("{:?}", std::time::Instant::now() - now);
        let edge_sequences = self.edge_consensus();
        let mut contig_number: HashMap<_, _> = (0..group.len()).map(|x| (x, 0)).collect();
        let name_of_nodes: HashMap<u64, _> = self
            .nodes
            .keys()
            .map(|id| {
                let group = group.get(id).unwrap();
                let node = contig_number.get_mut(group).unwrap();
                let name = format!("tig_{:04}_{:04}", group, node);
                *node += 1;
                (*id, name)
            })
            .collect();
        let mut segment_length: HashMap<_, usize> = HashMap::new();
        let segments: Vec<_> = self
            .nodes
            .iter()
            .map(|(key, node)| {
                let name = name_of_nodes[&key].clone();
                // Consensus.
                let mut seq = vec![];
                for w in node.windows(2) {
                    let (from, from_dir, to, _) = match *w {
                        [(ref f, fd), (ref t, td)] => (f, fd, t, td),
                        _ => panic!(),
                    };
                    if from.is_forward == from_dir {
                        seq.extend(&node_sequences[&(from.unit, from.cluster)]);
                    } else {
                        let rev = bio_utils::revcmp(&node_sequences[&(from.unit, from.cluster)]);
                        seq.extend(rev);
                    };
                    // TODO: Something's going wrong in this section....!!!!
                    if from.unit < to.unit {
                        let key = ((from.unit, from.cluster), (to.unit, to.cluster));
                        match edge_sequences.get(&key) {
                            Some(EdgeSeq::Offset(len)) => {
                                let len = (*len).min(seq.len() as i64);
                                (0..len).for_each(|_| {
                                    seq.pop().unwrap();
                                });
                            }
                            Some(EdgeSeq::Seq(ref xs)) => seq.extend(xs),
                            None => {}
                        }
                    } else {
                        let key = ((to.unit, to.cluster), (from.unit, from.cluster));
                        match edge_sequences.get(&key) {
                            Some(EdgeSeq::Offset(len)) => {
                                let len = (*len).min(seq.len() as i64);
                                (0..len).for_each(|_| {
                                    seq.pop().unwrap();
                                })
                            }
                            Some(EdgeSeq::Seq(ref xs)) => seq.extend(bio_utils::revcmp(xs)),
                            None => {}
                        }
                    }
                }
                let &(ref last, last_dir) = node.last().unwrap();
                if last.is_forward == last_dir {
                    seq.extend(&node_sequences[&(last.unit, last.cluster)]);
                } else {
                    let rev = bio_utils::revcmp(&node_sequences[&(last.unit, last.cluster)]);
                    seq.extend(rev);
                }
                let len = seq.len();
                segment_length.insert(key, len);
                gfa::Segment::from(name, len, Some(String::from_utf8(seq).unwrap()))
            })
            .collect();
        let edges: Vec<_> = self
            .edges
            .iter()
            .map(|x| x.1)
            .flatten()
            .enumerate()
            .map(|(idx, edge)| {
                let edgeid = Some(format!("edge_{}", idx));
                let sid1 = gfa::RefID::from(&name_of_nodes[&edge.from], true);
                let sid2 = gfa::RefID::from(&name_of_nodes[&edge.to], true);
                let len1 = segment_length[&edge.from];
                let len2 = segment_length[&edge.to];
                let beg1 = match edge.from_position {
                    Position::Head => gfa::Position::from(0, false),
                    Position::Tail => gfa::Position::from(len1, true),
                };
                let beg2 = match edge.to_position {
                    Position::Head => gfa::Position::from(0, false),
                    Position::Tail => gfa::Position::from(len2, true),
                };
                gfa::Edge::from(edgeid, sid1, sid2, beg1, beg1, beg2, beg2, None)
            })
            .collect();
        let group: Vec<_> = (0..cluster_num)
            .map(|gl| {
                let uid = Some(format!("group_{}", gl));
                let ids: Vec<_> = group
                    .iter()
                    .filter(|(_, val)| val == &&gl)
                    .map(|(key, _)| name_of_nodes[key].clone())
                    .collect();
                gfa::Group::Set(gfa::UnorderedGroup { uid, ids })
            })
            .collect();
        let read_ids: Vec<u64> = self.reads.iter().map(|r| r.id).collect();
        let summaries: Vec<ContigSummary> = name_of_nodes
            .iter()
            .map(|(&id, name)| {
                let parent = self.parent.find(id as usize).unwrap();
                let reads: Vec<_> = read_ids
                    .iter()
                    .filter(|&&id| self.parent.find(id as usize) == Some(parent))
                    .copied()
                    .collect();
                ContigSummary {
                    id: name.clone(),
                    reads,
                }
            })
            .collect();
        (segments, edges, group, summaries)
    }
}

enum EdgeSeq {
    Offset(i64),
    Seq(Vec<u8>),
}

#[derive(Debug, Clone, Copy)]
pub enum AlignmentResult {
    HeadToTail(usize),
    TailToHead(usize),
    DiagonalMatch,
    Contained,
    Containing,
}

pub fn unit_alignment_naive(
    xs: &[(u64, u64)],
    ys: &[(u64, u64)],
    _c: &AssembleConfig,
) -> Option<AlignmentResult> {
    let mut dp = vec![vec![0; ys.len() + 1]; xs.len() + 1];
    for (i, (u1, c1)) in xs.iter().enumerate().map(|(i, &p)| (i + 1, p)) {
        for (j, (u2, c2)) in ys.iter().enumerate().map(|(j, &p)| (j + 1, p)) {
            let mat_score = match (u1 == u2, c1 == c2) {
                (true, true) => 1,
                (true, false) => -1,
                (false, _) => -1000,
            };
            dp[i][j] = (dp[i - 1][j - 1] + mat_score)
                .max(dp[i - 1][j] - 1)
                .max(dp[i][j - 1] - 1);
        }
    }
    let (opt_score, xs_end, ys_end) = (0..xs.len() + 1)
        .map(|i| (i, ys.len()))
        .chain((0..ys.len()).map(|j| (xs.len(), j)))
        .map(|(i, j)| (dp[i][j], i, j))
        .max_by_key(|x| x.0)
        .unwrap();
    let (mut i, mut j) = (xs_end, ys_end);
    let mut match_num = 0;
    while 0 < i && 0 < j {
        let (u1, c1) = xs[i - 1];
        let (u2, c2) = ys[j - 1];
        let mat_score = match (u1 == u2, c1 == c2) {
            (true, true) => 1,
            (true, false) => -1,
            (false, _) => -1000,
        };
        if dp[i][j] == dp[i - 1][j - 1] + mat_score {
            i -= 1;
            j -= 1;
            match_num += ((u1 == u2) && (c1 == c2)) as i32;
        } else if dp[i][j] == dp[i - 1][j] - 1 {
            i -= 1;
        } else if dp[i][j] == dp[i][j - 1] - 1 {
            j -= 1;
        } else {
            unreachable!()
        }
    }
    let (xs_start, ys_start) = (i, j);
    if match_num < 3 || xs_end == xs_start || ys_end == ys_start {
        return None;
    }
    let normed_score =
        opt_score as f64 / (((xs_end - xs_start) * (ys_end - ys_start)) as f64).sqrt();
    if normed_score < 0.8 {
        return None;
    }
    match (
        xs_start == 0,
        ys_start == 0,
        xs_end == xs.len(),
        ys_end == ys.len(),
    ) {
        (true, true, true, true) => Some(AlignmentResult::DiagonalMatch),
        (true, _, true, _) => Some(AlignmentResult::Contained),
        (_, true, _, true) => Some(AlignmentResult::Containing),
        (true, false, false, true) => Some(AlignmentResult::HeadToTail(xs_end)),
        (false, true, true, false) => Some(AlignmentResult::TailToHead(xs.len() - xs_start)),
        _ => None,
    }
}

// Ok(r_start_pos, q_end_pos, true, score) or ()
pub fn unit_alignment(
    xs: &[(u64, u64)],
    ys: &[(u64, u64)],
    c: &AssembleConfig,
    score: &HashMap<(u64, u64), f64>,
) -> Option<AlignmentResult> {
    // Filling DP cell (Locally, Overlapping.)
    let mut dp = vec![vec![0f64; ys.len() + 1]; xs.len() + 1];
    for (i, (u1, c1)) in xs.iter().enumerate().map(|(i, &p)| (i + 1, p)) {
        for (j, (u2, c2)) in ys.iter().enumerate().map(|(j, &p)| (j + 1, p)) {
            let mat_score = match (u1 == u2, c1 == c2) {
                (true, true) => score[&(u1, c1)],
                (true, false) => MISM_SCORE,
                (false, _) => MISM_UNIT,
            };
            dp[i][j] = (dp[i - 1][j - 1] + mat_score)
                .max(dp[i - 1][j] + GAP_SCORE)
                .max(dp[i][j - 1] + GAP_SCORE);
        }
    }
    let (opt_score, xs_end, ys_end) = (0..xs.len() + 1)
        .map(|i| (i, ys.len()))
        .chain((0..ys.len()).map(|j| (xs.len(), j)))
        .map(|(i, j)| (dp[i][j], i, j))
        .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
        .unwrap();
    let (mut i, mut j) = (xs_end, ys_end);
    let mut match_num = 0;
    while 0 < i && 0 < j {
        let (u1, c1) = xs[i - 1];
        let (u2, c2) = ys[j - 1];
        let mat_score = match (u1 == u2, c1 == c2) {
            (true, true) => score[&(u1, c1)],
            (true, false) => MISM_SCORE,
            (false, _) => MISM_UNIT,
        };
        if (dp[i][j] - (dp[i - 1][j - 1] + mat_score)).abs() < 0.00001 {
            i -= 1;
            j -= 1;
            match_num += ((u1 == u2) && (c1 == c2)) as i32;
        } else if (dp[i][j] - (dp[i - 1][j] + GAP_SCORE)).abs() < 0.00001 {
            i -= 1;
        } else if (dp[i][j] - (dp[i][j - 1] + GAP_SCORE)).abs() < 0.00001 {
            j -= 1;
        } else {
            unreachable!()
        }
    }
    let (xs_start, ys_start) = (i, j);
    if match_num < c.ovlp || xs_end == xs_start || ys_end == ys_start {
        return None;
    }
    let normed_score = opt_score / (((xs_end - xs_start) * (ys_end - ys_start)) as f64).sqrt();
    if normed_score < c.identity {
        return None;
    }
    match (
        xs_start == 0,
        ys_start == 0,
        xs_end == xs.len(),
        ys_end == ys.len(),
    ) {
        (true, true, true, true) => Some(AlignmentResult::DiagonalMatch),
        (true, _, true, _) => Some(AlignmentResult::Contained),
        (_, true, _, true) => Some(AlignmentResult::Containing),
        (true, false, false, true) => Some(AlignmentResult::HeadToTail(xs_end)),
        (false, true, true, false) => Some(AlignmentResult::TailToHead(xs.len() - xs_start)),
        _ => None,
    }
    // dp.iter()
    //     .enumerate()
    //     .flat_map(|(r_argmax, line)| {
    //         line.iter()
    //             .enumerate()
    //             .filter(|&(_, &score)| score == max)
    //             .filter_map(|(q_argmax, _)| {
    //                 validate_alignment(r_argmax, r, q_argmax, q, &dp, max, c.identity)
    //             })
    //             .collect::<Vec<_>>()
    //     })
    //     .collect()
}

// fn validate_alignment(
//     r_argmax: usize,
//     r: &[(u64, u64)],
//     q_argmax: usize,
//     q: &[(u64, u64)],
//     dp: &[Vec<i32>],
//     max: i32,
//     iden: f64,
// ) -> Option<AlignmentResult> {
//     let (mut rpos, mut qpos) = (r_argmax, q_argmax);
//     while 0 < rpos && 0 < qpos {
//         let mat = if r[rpos - 1] == q[qpos - 1] { 1 } else { -1 };
//         if dp[rpos][qpos] == mat + dp[rpos - 1][qpos - 1] {
//             rpos -= 1;
//             qpos -= 1;
//         } else if dp[rpos][qpos] == dp[rpos - 1][qpos] - 1 {
//             rpos -= 1;
//         } else if dp[rpos][qpos] == dp[rpos][qpos - 1] - 1 {
//             qpos -= 1;
//         } else {
//             assert_eq!(dp[rpos][qpos], 0);
//             break;
//         }
//     }
//     // First condition:
//     // If this alignment truly overlaps,
//     // the optimal alignment should starts from the start position of either of reads.
//     // Second condition:
//     // If this alignment truly overlaps,
//     // the optimal alignment should ends at either of reads.
//     // eprintln!("{}\t{}", r.len(), q.len());
//     let is_valid_overlap = (qpos == 0 || rpos == 0) && (q_argmax == q.len() || r_argmax == r.len());
//     // eprintln!(
//     //     "{}-{}\t{}-{}\t{}",
//     //     rpos, r_argmax, qpos, q_argmax, is_valid_overlap
//     // );
//     if !is_valid_overlap {
//         None
//     } else if rpos == 0 && qpos == 0 && q_argmax == q.len() && r_argmax == r.len() {
//         Some(AlignmentResult::DiagonalMatch)
//     } else if rpos == 0 && r_argmax == r.len() {
//         Some(AlignmentResult::Contained)
//     } else if qpos == 0 && q_argmax == q.len() {
//         Some(AlignmentResult::Containing)
//     } else if rpos == 0 {
//         assert!(r_argmax != r.len());
//         let (rlen, qlen) = (r_argmax, q.len() - qpos);
//         let identity = max as f64 / ((rlen * qlen) as f64).sqrt();
//         // eprintln!("{}\t{}\t{}", rlen, qlen, identity);
//         if iden < identity {
//             Some(AlignmentResult::HeadToTail(r_argmax))
//         } else {
//             None
//         }
//     } else {
//         assert!(qpos == 0 && q_argmax != q.len());
//         let (rlen, qlen) = (r.len() - rpos, q_argmax);
//         let identity = max as f64 / ((rlen * qlen) as f64).sqrt();
//         // eprintln!("{}\t{}\t{}", rlen, qlen, identity);
//         if iden < identity {
//             Some(AlignmentResult::TailToHead(r.len() - rpos))
//         } else {
//             None
//         }
//     }
// }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrEdge {
    from: u64,
    to: u64,
    from_position: Position,
    to_position: Position,
    offset: usize,
}

use crate::global_clustering::error_correction::CorrectedRead;
impl StrEdge {
    #[allow(dead_code)]
    fn from_corrected_reads(
        r: &CorrectedRead,
        q: &CorrectedRead,
        c: &AssembleConfig,
        score: &HashMap<(u64, u64), f64>,
    ) -> Result<Vec<Self>, bool> {
        let mut rseq: Vec<_> = r.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        let qseq: Vec<_> = q.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        use AlignmentResult::*;
        use Position::*;
        let mut edges = vec![];
        for aln in unit_alignment(&rseq, &qseq, c, score) {
            match aln {
                HeadToTail(ofs) => edges.push(StrEdge::new(r.id, q.id, Head, Tail, ofs)),
                TailToHead(ofs) => edges.push(StrEdge::new(r.id, q.id, Tail, Head, ofs)),
                Contained => return Err(true),
                Containing => return Err(false),
                DiagonalMatch if r.id < q.id => return Err(true),
                DiagonalMatch => return Err(false),
            }
        }
        rseq.reverse();
        for aln in unit_alignment(&rseq, &qseq, c, score) {
            match aln {
                HeadToTail(ofs) => edges.push(StrEdge::new(r.id, q.id, Tail, Tail, ofs)),
                TailToHead(ofs) => edges.push(StrEdge::new(r.id, q.id, Head, Head, ofs)),
                Contained => return Err(true),
                Containing => return Err(false),
                DiagonalMatch if r.id < q.id => return Err(true),
                DiagonalMatch => return Err(false),
            }
        }
        Ok(edges)
    }
    fn from_reads(
        r: &EncodedRead,
        q: &EncodedRead,
        c: &AssembleConfig,
        score: &HashMap<(u64, u64), f64>,
    ) -> Result<Self, bool> {
        let mut rseq: Vec<_> = r.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        let qseq: Vec<_> = q.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        use AlignmentResult::*;
        use Position::*;
        if let Some(aln) = unit_alignment(&rseq, &qseq, c, score) {
            // if let Some(aln) = unit_alignment_naive(&rseq, &qseq, c) {
            match aln {
                HeadToTail(ofs) => return Ok(StrEdge::new(r.id, q.id, Head, Tail, ofs)),
                TailToHead(ofs) => return Ok(StrEdge::new(r.id, q.id, Tail, Head, ofs)),
                Contained => return Err(true),
                Containing => return Err(false),
                DiagonalMatch if r.id < q.id => return Err(true),
                DiagonalMatch => return Err(false),
            }
        }
        rseq.reverse();
        if let Some(aln) = unit_alignment(&rseq, &qseq, c, score) {
            // if let Some(aln) = unit_alignment_naive(&rseq, &qseq, c) {
            match aln {
                HeadToTail(ofs) => Ok(StrEdge::new(r.id, q.id, Tail, Tail, ofs)),
                TailToHead(ofs) => Ok(StrEdge::new(r.id, q.id, Head, Head, ofs)),
                Contained => Err(true),
                Containing => Err(false),
                DiagonalMatch if r.id < q.id => Err(true),
                DiagonalMatch => Err(false),
            }
        } else {
            Err(false)
        }
        // Ok(edges)
    }
    fn new(
        from: u64,
        to: u64,
        from_position: Position,
        to_position: Position,
        offset: usize,
    ) -> Self {
        Self {
            from,
            to,
            from_position,
            to_position,
            offset,
        }
    }
    fn is_transitive(&self, beta: &Self, gamma: &Self) -> bool {
        let (node_a, node_b, node_c) = (self.from, self.to, beta.to);
        assert_eq!(node_a, beta.from);
        assert_eq!(node_b, gamma.from);
        assert_eq!(node_c, gamma.to);
        gamma.to_position != beta.to_position
            && gamma.from_position == self.to_position
            && self.from_position == beta.from_position
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Position {
    Head,
    Tail,
}

impl std::ops::Not for Position {
    type Output = Self;
    fn not(self) -> Self::Output {
        match self {
            Self::Head => Self::Tail,
            Self::Tail => Self::Head,
        }
    }
}

impl Position {}

#[derive(Debug, Clone)]
pub struct ContigSummary {
    pub id: String,
    /// ID of the reads consisting this contig.
    pub reads: Vec<u64>,
}

fn get_match_score(ds: &DataSet) -> HashMap<(u64, u64), f64> {
    let mut counts: HashMap<u64, HashMap<u64, u32>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *counts
                .entry(node.unit)
                .or_default()
                .entry(node.cluster)
                .or_default() += 1;
        }
    }
    counts
        .iter()
        .flat_map(|(&unit, val)| {
            let total: u32 = val.values().sum();
            val.iter()
                .map(|(&cluster, &count)| {
                    let frac = count as f64 / total as f64;
                    let score = (MATCH + MISMATCH * frac).ln() - frac.ln();
                    ((unit, cluster), score)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}
