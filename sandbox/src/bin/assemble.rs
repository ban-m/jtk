use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let start = std::time::Instant::now();
    let ds = get_input_file()?;
    // let graph = assemble(&ds.reads);
    let end = std::time::Instant::now();
    eprintln!("{:?}", end - start);
    let mut graph = StringGraph::new(&ds.encoded_reads);
    graph.transitive_edge_reduction();
    graph.simple_path_reduction();
    let result = graph.assemble_as_gfa();
    println!("{}", result);
    // let mut r: EncodedRead = EncodedRead::default();
    // r.nodes = (5..15)
    //     .map(|i| Node {
    //         position_from_start: 0,
    //         unit: i as u64,
    //         cluster: 0,
    //         seq: String::new(),
    //         is_forward: false,
    //         cigar: vec![],
    //     })
    //     .collect();
    // r.edges = r
    //     .nodes
    //     .windows(2)
    //     .map(|w| Edge {
    //         from: w[0].unit,
    //         to: w[1].unit,
    //         offset: 0,
    //         label: String::new(),
    //     })
    //     .collect();
    // let mut q: EncodedRead = EncodedRead::default();
    // q.nodes = (0..10)
    //     .rev()
    //     .map(|i| Node {
    //         position_from_start: 0,
    //         unit: i as u64,
    //         cluster: 0,
    //         seq: String::new(),
    //         is_forward: false,
    //         cigar: vec![],
    //     })
    //     .collect();
    // q.edges = q
    //     .nodes
    //     .windows(2)
    //     .map(|w| Edge {
    //         from: w[0].unit,
    //         to: w[1].unit,
    //         offset: 0,
    //         label: String::new(),
    //     })
    //     .collect();
    // let aln = StrEdge::from_reads(&q, &r);
    // eprintln!("{:?}", aln);
    Ok(())
}

fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

// TODO: Convert them into Vector type.
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone)]
pub struct StringGraph<'a> {
    pub edges: HashMap<u64, Vec<StrEdge>>,
    pub nodes: HashMap<u64, Vec<(u64, u64)>>,
    pub reads: &'a [EncodedRead],
}

impl<'a> std::fmt::Display for StringGraph<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let n = self.reads.len();
        let e = self.edges.len();
        write!(f, "N:{}\tE:{}", n, e)
    }
}

impl<'a> StringGraph<'a> {
    pub fn new(reads: &'a [EncodedRead]) -> Self {
        // If Some(aln), alignment, If None, it is contained to some read.
        let edges: HashMap<_, Option<Vec<_>>> = reads
            .iter()
            .map(|r| {
                let mut edges = vec![];
                for q in reads.iter().filter(|q| q.id != r.id) {
                    match StrEdge::from_reads(r, q) {
                        Ok(res) => edges.push(res),
                        Err(is_contained) if is_contained => return (r.id, None),
                        Err(_) => {}
                    }
                }
                (r.id, Some(edges))
            })
            .collect();
        eprintln!("{} Reads.", reads.len());
        let num_edges = edges
            .values()
            .filter_map(|edges| edges.as_ref().map(|x| x.len()))
            .sum::<usize>();
        eprintln!("{} Alignments.", num_edges);
        // Contained read removal.
        let is_contained: HashSet<u64> = edges
            .iter()
            .filter_map(|(k, v)| if v.is_none() { Some(*k) } else { None })
            .collect();
        eprintln!("{} Contained reads.", is_contained.len());
        let edges: HashMap<_, Vec<_>> = edges
            .into_iter()
            .filter_map(|(k, v)| {
                v.and_then(|mut edges| {
                    edges.retain(|aln| {
                        !is_contained.contains(&aln.from) && !is_contained.contains(&aln.to)
                    });
                    if edges.is_empty() {
                        None
                    } else {
                        assert!(edges.iter().all(|e| e.from == k));
                        edges.sort_by_key(|e| e.to);
                        Some((k, edges))
                    }
                })
            })
            .collect();
        let is_used: HashSet<u64> = {
            let mut is_used = HashSet::new();
            for edges in edges.values() {
                for edge in edges.iter() {
                    is_used.insert(edge.from);
                    is_used.insert(edge.to);
                }
            }
            is_used
        };
        eprintln!("{} reads are used.", is_used.len());
        let nodes: HashMap<u64, Vec<_>> = reads
            .iter()
            .filter(|r| !is_used.contains(&r.id))
            .map(|r| (r.id, r.nodes.iter().map(|n| (n.unit, n.cluster)).collect()))
            .collect();
        // We save reads anyway.
        let num_edges = edges.values().map(|x| x.len()).sum::<usize>();
        eprintln!("{} edges.", num_edges);
        Self {
            reads,
            edges,
            nodes,
        }
    }
    pub fn transitive_edge_reduction(&mut self) {
        let is_transitive_edge: HashMap<u64, Vec<bool>> = self
            .edges
            .iter()
            .map(|(&k, edges)| {
                let is_transitive: Vec<_> = edges
                    .iter()
                    .map(|edge| {
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
            if let Some(transitive_flags) = is_transitive_edge.get(k) {
                let mut idx = 0;
                v.retain(|_| {
                    idx += 1;
                    transitive_flags[idx - 1]
                });
            }
        });
        let after_edges = self.edges.values().map(|x| x.len()).sum::<usize>();
        eprintln!("{}=>{}", num_edges, after_edges);
        // Sanity check.
        for (_k, edges) in self.edges.iter() {
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
    pub fn assemble_as_gfa(&mut self) -> gfa::GFA {
        use gfa::*;
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let (nodes, edges, group) = self.assemble();
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
    pub fn simple_path_reduction(&mut self) {
        let ids: Vec<_> = self.nodes.keys().copied().collect();
        for id in ids {
            let len = self.edges[&id].len();
            let mut is_marged = vec![];
            for idx in 0..len {
                let edge = self.edges[&id][idx];
                let num_edges = self.edges[&edge.to]
                    .iter()
                    .filter(|e| e.from_position == edge.to_position)
                    .count();
                if num_edges != 1 || edge.to == edge.from {
                    is_marged.push(false);
                } else {
                    // Never panic.
                    {
                        let mut poped_seq = self.nodes.remove(&edge.to).unwrap();
                        let edge_seq = self.nodes.get_mut(&edge.from).unwrap();
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
                                edge_seq.extend(poped_seq.into_iter().rev());
                            }
                            (Head, Tail) => {
                                poped_seq.extend(edge_seq.iter().copied().skip(edge.offset));
                                edge_seq.clear();
                                edge_seq.extend(poped_seq);
                            }
                            (Head, Head) => {
                                let mut ps = poped_seq;
                                ps.reverse();
                                ps.extend(edge_seq.iter().copied().skip(edge.offset));
                                edge_seq.clear();
                                edge_seq.extend(ps);
                            }
                        }
                    }
                    // Remove reverse edge. Never panic.
                    let mut new_edges = self.edges.remove(&edge.to).unwrap();
                    new_edges.retain(|e| e.from_position != edge.to_position);
                    new_edges.iter_mut().for_each(|e| {
                        self.edges
                            .get_mut(&e.to)
                            .unwrap()
                            .iter_mut()
                            .filter(|f| f.to == e.from)
                            .for_each(|f| {
                                f.to = edge.from;
                                f.to_position = edge.from_position;
                            });
                        e.from = edge.from;
                        e.from_position = edge.from_position;
                    });
                    is_marged.push(true);
                }
            }
            if let Some(edges) = self.edges.get_mut(&id) {
                let mut idx = 0;
                edges.retain(|_| {
                    idx += 1;
                    is_marged[idx - 1]
                });
            }
        }
    }
    fn assemble(&mut self) -> (Vec<gfa::Segment>, Vec<gfa::Edge>, Vec<gfa::Group>) {
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
            eprintln!("{}", current_cluster);
            (group, current_cluster)
        };
        let mut contig_number: HashMap<_, _> = (0..group.len()).map(|x| (x, 0)).collect();
        let name_of_nodes: HashMap<_, _> = self
            .nodes
            .keys()
            .map(|id| {
                let group = group.get(id).unwrap();
                let node = contig_number.get_mut(group).unwrap();
                let name = format!("tig_{:04}_{:04}", group, node);
                *node += 1;
                (id, name)
            })
            .collect();
        let segments: Vec<_> = self
            .nodes
            .iter()
            .map(|(key, node)| {
                let name = name_of_nodes[&key].clone();
                gfa::Segment::from(name, 2000 * node.len(), None)
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
                let len1 = self.nodes[&edge.from].len() * 2000;
                let len2 = self.nodes[&edge.to].len() * 2000;
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
        (segments, edges, group)
    }
}

enum AlignmentResult {
    HeadToTail(usize),
    TailToHead(usize),
    DiagonalMatch,
    Contained,
    Containing,
    LowScore,
}

// Ok(r_start_pos, q_end_pos, true, score) or ()
fn unit_alignment(r: &[(u64, u64)], q: &[(u64, u64)], thr: i32) -> AlignmentResult {
    let mut dp = vec![vec![0; q.len() + 1]; r.len() + 1];
    // Filling DP cell
    for (i, r_base) in r.iter().enumerate() {
        for (j, q_base) in q.iter().enumerate() {
            let mat = if r_base == q_base { 1 } else { -1 };
            dp[i + 1][j + 1] = (dp[i][j] + mat).max(dp[i][j + 1] - 1).max(dp[i + 1][j] - 1);
        }
    }
    // Traceback
    let (q_argmax, q_max) = dp
        .last()
        .unwrap()
        .iter()
        .enumerate()
        .max_by_key(|x| x.1)
        .unwrap();
    let (r_argmax, r_max) = dp
        .iter()
        .filter_map(|line| line.last())
        .enumerate()
        .max_by_key(|x| x.1)
        .unwrap();
    // let rs: Vec<_> = r.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
    // let qs: Vec<_> = q.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
    // eprintln!("{}\n{}", rs.join(":"), qs.join(":"));
    // eprintln!("{}\t{}\t{}\t{}", q_argmax, q_max, r_argmax, r_max);
    let (mut rpos, mut qpos) = if r_max <= q_max {
        (r.len(), q_argmax)
    } else {
        (r_argmax, q.len())
    };
    while 0 < rpos && 0 < qpos {
        let mat = if r[rpos - 1] == q[qpos - 1] { 1 } else { -1 };
        if dp[rpos][qpos] == mat + dp[rpos - 1][qpos - 1] {
            rpos -= 1;
            qpos -= 1;
        } else if dp[rpos][qpos] == dp[rpos - 1][qpos] - 1 {
            rpos -= 1;
        } else {
            assert_eq!(dp[rpos][qpos], dp[rpos][qpos - 1] - 1);
            qpos -= 1;
        }
    }
    if r_max < &thr && q_max < &thr {
        AlignmentResult::LowScore
    } else if rpos == 0 && qpos == 0 && q_argmax == q.len() && r_argmax == r.len() {
        let rs: Vec<_> = r.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
        let qs: Vec<_> = q.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
        eprintln!("Diag:\n{}\n{}", rs.join(":"), qs.join(":"));
        AlignmentResult::DiagonalMatch
    } else if r_max <= q_max {
        if rpos == 0 {
            let rs: Vec<_> = r.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
            let qs: Vec<_> = q.iter().map(|(u, c)| format!("{}-{}", u, c)).collect();
            eprintln!("Contained:\n{}\n{}", rs.join(":"), qs.join(":"));
            AlignmentResult::Contained
        } else {
            AlignmentResult::TailToHead(r.len() - rpos)
        }
    } else {
        if qpos == 0 {
            AlignmentResult::Containing
        } else {
            AlignmentResult::HeadToTail(r_argmax)
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct StrEdge {
    from: u64,
    to: u64,
    from_position: Position,
    to_position: Position,
    offset: usize,
}

impl StrEdge {
    fn from_reads(r: &EncodedRead, q: &EncodedRead) -> Result<Self, bool> {
        let mut rseq: Vec<_> = r.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        let qseq: Vec<_> = q.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
        use AlignmentResult::*;
        use Position::*;
        match unit_alignment(&rseq, &qseq, 4) {
            HeadToTail(ofs) => Ok(StrEdge::new(r.id, q.id, Head, Tail, ofs)),
            TailToHead(ofs) => Ok(StrEdge::new(r.id, q.id, Tail, Head, ofs)),
            Contained => Err(true),
            Containing => Err(false),
            DiagonalMatch if r.id < q.id => Err(true),
            DiagonalMatch => Err(false),
            LowScore => {
                rseq.reverse();
                match unit_alignment(&rseq, &qseq, 4) {
                    HeadToTail(ofs) => Ok(StrEdge::new(r.id, q.id, Tail, Tail, ofs)),
                    TailToHead(ofs) => Ok(StrEdge::new(r.id, q.id, Head, Head, ofs)),
                    Contained => Err(true),
                    Containing => Err(false),
                    DiagonalMatch if r.id < q.id => Err(true),
                    DiagonalMatch => Err(false),
                    LowScore => Err(false),
                }
            }
        }
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
