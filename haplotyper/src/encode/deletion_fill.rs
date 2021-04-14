//! Filling deletion.
use definitions::*;
use rayon::prelude::*;
// use log::*;
use std::collections::HashMap;
/// Fill deletions
pub fn correct_unit_deletion(mut ds: DataSet) -> DataSet {
    let read_skeltons: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|read| read.nodes.len() > 1)
        .map(|read| ReadSkelton::from_rich_nodes(&read.nodes))
        .collect();
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.seq()))
        .collect();
    let selected_chunks = &ds.selected_chunks;
    ds.encoded_reads
        .par_iter_mut()
        .filter(|r| r.nodes.len() > 1)
        .for_each(|r| {
            let seq: Vec<_> = raw_seq[&r.id]
                .iter()
                .map(|x| x.to_ascii_uppercase())
                .collect();
            correct_deletion_error(r, selected_chunks, &read_skeltons, &seq);
        });
    ds
}

// Aligment offset. We align [s-offset..e+offset] region to the unit.
const OFFSET: usize = 200;
fn correct_deletion_error(
    read: &mut EncodedRead,
    units: &[Unit],
    reads: &[ReadSkelton],
    seq: &[u8],
) {
    // Insertion counts.
    let pileups = get_pileup(read, reads);
    let threshold = get_threshold(&pileups);
    let nodes = &read.nodes;
    let take_len = nodes.len();
    let inserts: Vec<_> = pileups
        .iter()
        .enumerate()
        .take(take_len)
        .skip(1)
        .filter_map(|(idx, pileup)| {
            let max_num = pileup.max_insertion();
            if max_num < threshold {
                return None;
            }
            // Check if we can record new units.
            // It never panics.
            let (prev_offset, uid, direction, after_offset) = pileup.max_ins_information().unwrap();
            let start_position =
                (nodes[idx - 1].position_from_start + nodes[idx - 1].query_length()) as isize;
            let start_position = start_position + prev_offset;
            let start_position = (start_position as usize).saturating_sub(OFFSET);
            let end_position = nodes[idx].position_from_start as isize - after_offset;
            let end_position = (end_position as usize + OFFSET).min(seq.len());
            // Never panic.
            assert!(start_position < end_position,);
            let unit = units.iter().find(|u| u.id == uid).unwrap();
            encode_node(&seq, start_position, end_position, direction, unit)
                .map(|new_node| (idx, new_node))
        })
        .collect();
    let mut accum_insertes = 0;
    for (idx, node) in inserts {
        read.nodes.insert(idx + accum_insertes, node);
        accum_insertes += 1;
    }
    read.edges = read
        .nodes
        .windows(2)
        .map(|w| Edge::from_nodes(w, &seq))
        .collect();
}

// 0.3 * unit.seq() distance is regarded as too far: corresponding 30% errors.
const ALIGN_LIMIT: f64 = 0.3;
// Try to Encode Node. Return Some(node) if the alignment is good.
fn encode_node(
    seq: &[u8],
    start: usize,
    end: usize,
    is_forward: bool,
    unit: &Unit,
) -> Option<Node> {
    let query = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    // TODO: Maybe more sophisticated, or dedicated alignmnet parameters should be used.
    // Or, maybe we should use affine gap panalty version of semi-global alignment.
    let (dist, ops) = kiley::bialignment::edit_dist_slow_ops_semiglobal(unit.seq(), &query);
    let dist_thr = (unit.seq().len() as f64 * ALIGN_LIMIT).floor() as u32;
    // Leading/Trailing bases are offsets.
    if dist_thr + 2 * (OFFSET as u32) < dist {
        return None;
    };
    let mut ops = super::compress_kiley_ops(&ops);
    // How many bases are regarded as insertion from `start` position.
    // If the sequence is revcmped, it is lengt of the unnesesarry trailing sequence.
    // Remove leading insertions.
    let removed_base_from_start = if let Some(&Op::Ins(l)) = ops.first() {
        ops.remove(0);
        l
    } else {
        0
    };
    // Remove trailing insertions.
    let removed_base_from_end = if let Some(&Op::Ins(l)) = ops.last() {
        ops.pop();
        l
    } else {
        0
    };
    // Check whether there is too large indels.
    if ops.iter().any(|op| match *op {
        Op::Ins(l) | Op::Del(l) => l > super::INDEL_THRESHOLD,
        _ => false,
    }) {
        return None;
    }
    // Now we know this alignment is valid.
    let (start, end) = if is_forward {
        (start + removed_base_from_start, end - removed_base_from_end)
    } else {
        (start + removed_base_from_end, end - removed_base_from_start)
    };
    let seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    unit.seq().iter().all(|x| x.is_ascii_uppercase());
    seq.iter().all(|x| x.is_ascii_uppercase());
    Some(Node {
        position_from_start: start,
        unit: unit.id,
        cluster: 0,
        is_forward,
        seq: String::from_utf8(seq).unwrap(),
        cigar: ops,
    })
}

// Align read skeltons to read, return the pileup sumamries.
fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    for query in reads.iter() {
        if let Some((aln, is_forward)) = alignment(read, query) {
            let query = if is_forward {
                query.clone()
            } else {
                query.rev()
            };
            // debug!("{:?}", query);
            // debug!("{:?}", ReadSkelton::from_rich_nodes(&read.nodes));
            // debug!("{:?}", aln);
            let (mut r_ptr, mut q_ptr) = (0, 0);
            for op in aln {
                // We only record insertions with the size of 1, to
                // remove ambiguity.
                match op {
                    Op::Ins(l) => {
                        if r_ptr < read.nodes.len() && 0 < r_ptr && l == 1 {
                            for i in q_ptr..q_ptr + l {
                                pileups[r_ptr].add(query.nodes[i].clone());
                            }
                        }
                        q_ptr += l;
                    }
                    Op::Del(l) => {
                        r_ptr += l;
                    }
                    Op::Match(l) => {
                        for i in r_ptr..r_ptr + l {
                            pileups[i].coverage += 1;
                        }
                        r_ptr += l;
                        q_ptr += l;
                    }
                }
            }
        }
    }
    pileups
}

// Maybe we should tune this.
// For example, is it ok to use these parameters to treat long repeats?
// Maybe OK, as we confirm these candidate by alignment.
// Minimum required units to be matched.
const MIN_MATCH: usize = 2;
// Minimum required alignment score.
const SCORE_THR: i32 = 1;
// Alignment to the best direction, return the cigar and if the query should be reversed.
// (*false* if needed).
fn alignment(read: &EncodedRead, query: &ReadSkelton) -> Option<(Vec<Op>, bool)> {
    let read = ReadSkelton::from_rich_nodes(&read.nodes);
    use std::collections::HashSet;
    let read_nodes: HashSet<_> = read.nodes.iter().map(|n| n.unit).collect();
    if query.nodes.iter().all(|n| !read_nodes.contains(&n.unit)) {
        return None;
    }
    let (f_score, f_ops) = pairwise_alignment(&read, query);
    let f_match = get_match_units(&f_ops);
    let query = query.rev();
    let (r_score, r_ops) = pairwise_alignment(&read, &query);
    let r_match = get_match_units(&r_ops);
    // Dovetails should be proper. Not Del->In or In -> Del transition should reside.
    if r_score <= f_score && MIN_MATCH < f_match && SCORE_THR < f_score && is_proper(&f_ops) {
        return Some((f_ops, true));
    }
    if f_score <= r_score && MIN_MATCH < r_match && SCORE_THR < r_score && is_proper(&r_ops) {
        return Some((r_ops, false));
    }
    None
}

// Return true if the alignment is proper dovetail.
fn is_proper(ops: &[Op]) -> bool {
    ops.windows(2).all(|xs| match xs {
        &[Op::Ins(_), Op::Del(_)] | &[Op::Del(_), Op::Ins(_)] => false,
        _ => true,
    })
}

// This should return overlapping alignment and its score.
fn pairwise_alignment(read: &ReadSkelton, query: &ReadSkelton) -> (i32, Vec<Op>) {
    let (read, query) = (&read.nodes, &query.nodes);
    // We do not allow any mismatch by restricting the mismatch score to read.len() + query.len(),
    // which is apparently an upperbound.
    let score = |x: &LightNode, y: &LightNode| -> i32 {
        if x.unit == y.unit && x.is_forward == y.is_forward {
            1
        } else {
            -1 * (read.len() + query.len()) as i32
        }
    };
    // Fill DP cells.
    let mut dp = vec![vec![0; query.len() + 1]; read.len() + 1];
    for (i, x) in read.iter().enumerate() {
        for (j, y) in query.iter().enumerate() {
            let (i, j) = (i + 1, j + 1);
            let match_score = dp[i - 1][j - 1] + score(x, y);
            let ins = dp[i][j - 1] - 1;
            let del = dp[i - 1][j] - 1;
            dp[i][j] = match_score.max(del).max(ins);
        }
    }
    let (mut r_pos, mut q_pos, dist) = (0..read.len() + 1)
        .map(|i| (i, query.len()))
        .chain((0..query.len() + 1).map(|j| (read.len(), j)))
        .map(|(i, j)| (i, j, dp[i][j]))
        .max_by_key(|x| x.2)
        .unwrap();
    // We encode alignment by silly coding, in which every operation has length 1 and
    // successive operation might be the same.
    let mut ops = vec![];
    if read.len() != r_pos {
        assert_eq!(q_pos, query.len());
        ops.push(Op::Del(read.len() - r_pos));
    }
    if query.len() != q_pos {
        assert_eq!(r_pos, read.len());
        ops.push(Op::Ins(query.len() - q_pos));
    }
    while 0 < r_pos && 0 < q_pos {
        let current_dist = dp[r_pos][q_pos];
        if current_dist == dp[r_pos - 1][q_pos] - 1 {
            ops.push(Op::Del(1));
            r_pos -= 1;
        } else if current_dist == dp[r_pos][q_pos - 1] - 1 {
            ops.push(Op::Ins(1));
            q_pos -= 1;
        } else {
            let match_score = dp[r_pos - 1][q_pos - 1] + score(&read[r_pos - 1], &query[q_pos - 1]);
            assert_eq!(match_score, current_dist);
            ops.push(Op::Match(1));
            q_pos -= 1;
            r_pos -= 1;
        }
    }
    assert!(r_pos == 0 || q_pos == 0);
    if r_pos != 0 {
        ops.push(Op::Del(r_pos));
    }
    if q_pos != 0 {
        ops.push(Op::Ins(q_pos));
    }
    ops.reverse();
    let ops = compress_operations(ops);
    (dist, ops)
}

fn compress_operations(ops: Vec<Op>) -> Vec<Op> {
    assert!(!ops.is_empty());
    let mut current_op = ops[0];
    let mut compressed = vec![];
    for &op in ops.iter().skip(1) {
        match (op, current_op) {
            (Op::Match(l), Op::Match(m)) => {
                current_op = Op::Match(l + m);
            }
            (Op::Ins(l), Op::Ins(m)) => {
                current_op = Op::Ins(l + m);
            }
            (Op::Del(l), Op::Del(m)) => {
                current_op = Op::Del(l + m);
            }
            (x, _) => {
                compressed.push(current_op);
                current_op = x;
            }
        }
    }
    compressed.push(current_op);
    compressed
}

fn get_match_units(ops: &[Op]) -> usize {
    ops.iter()
        .map(|op| match op {
            Op::Match(l) => *l,
            _ => 0,
        })
        .sum::<usize>()
}

// Get threshold. In other words, a position would be regarded as an insertion if the
// count for a inserted unit is more than the return value of this function.
fn get_threshold(pileups: &[Pileup]) -> usize {
    let totcov = pileups.iter().map(|p| p.coverage).sum::<usize>();
    // We need at least 3 insertions to confirm.
    (totcov / 2 / pileups.len()).max(3)
}

#[derive(Debug, Clone)]
pub struct Pileup {
    inserted: Vec<LightNode>,
    coverage: usize,
}

impl Pileup {
    // Return the maximum insertion from the same unit, the same direction.
    fn max_insertion(&self) -> usize {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        count.values().max().copied().unwrap_or(0)
    }
    fn new() -> Self {
        Self {
            inserted: vec![],
            coverage: 0,
        }
    }
    fn max_ins_information(&self) -> Option<(isize, u64, bool, isize)> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        let ((max_unit, max_dir), _) = count.into_iter().max_by_key(|x| x.1)?;
        let inserts = self
            .inserted
            .iter()
            .filter(|node| node.unit == max_unit && node.is_forward == max_dir);
        let prev_offset = {
            let (count, total) =
                inserts
                    .clone()
                    .fold((0, 0), |(count, total), node| match node.prev_offset {
                        Some(len) => (count + 1, total + len),
                        None => (count, total),
                    });
            total / count
        };
        let after_offset = {
            let (count, total) =
                inserts
                    .clone()
                    .fold((0, 0), |(count, total), node| match node.after_offset {
                        Some(len) => (count + 1, total + len),
                        None => (count, total),
                    });
            total / count
        };
        Some((prev_offset, max_unit, max_dir, after_offset))
    }
    fn add(&mut self, node: LightNode) {
        self.inserted.push(node);
    }
}

#[derive(Clone)]
pub struct ReadSkelton {
    nodes: Vec<LightNode>,
}

impl std::fmt::Debug for ReadSkelton {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for node in self.nodes.iter() {
            write!(f, "{}:", node.unit)?;
        }
        Ok(())
    }
}

impl ReadSkelton {
    fn from_rich_nodes(nodes: &[Node]) -> Self {
        // Convert the nodes into (start_position, end_position)s
        let summaries: Vec<_> = nodes
            .iter()
            .map(|node| {
                let start = node.position_from_start;
                let end = start + node.query_length();
                (start as isize, end as isize)
            })
            .collect();
        let nodes: Vec<_> = nodes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                let prev_end = if i == 0 { None } else { summaries.get(i - 1) };
                let prev_offset = prev_end.map(|x| summaries[i].0 as isize - x.1 as isize);
                let after_offset = summaries.get(i + 1).map(|x| x.0 - summaries[i].1);
                LightNode {
                    prev_offset,
                    unit: n.unit,
                    is_forward: n.is_forward,
                    after_offset,
                }
            })
            .collect();
        ReadSkelton { nodes }
    }
    fn rev(&self) -> Self {
        let nodes: Vec<_> = self.nodes.iter().rev().map(LightNode::rev).collect();
        Self { nodes }
    }
}

#[derive(Clone)]
pub struct LightNode {
    // How long should be add to th last position of the previous node to
    // get the start position of this node.
    // None if this is the first node.
    prev_offset: Option<isize>,
    unit: u64,
    is_forward: bool,
    // Almost the same as prev_offset. The distance between the last postion of this node to
    // the start position of the next node.
    // None if this is the last node.
    after_offset: Option<isize>,
}

impl std::fmt::Debug for LightNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_forward {
            write!(f, "{}(+)", self.unit)
        } else {
            write!(f, "{}(-)", self.unit)
        }
    }
}

impl LightNode {
    fn rev(
        &Self {
            prev_offset,
            unit,
            is_forward,
            after_offset,
        }: &Self,
    ) -> Self {
        Self {
            prev_offset: after_offset,
            unit,
            is_forward: !is_forward,
            after_offset: prev_offset,
        }
    }
}
