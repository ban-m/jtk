//! Filling deletion.
use definitions::*;
use rayon::prelude::*;
// use log::*;
use std::collections::{HashMap, HashSet};

// The second argument is the vector of (index,unit_id) of the previous failed trials.
// for example, if failed_trials[i][0] = (j,id), then, we've already tried to encode the id-th unit after the j-th
// position of the i-th read, and failed it.
// If we can encode some position in the i-th read, the failed trials would be erased, as it change the
// condition of the read, making it possible to encode an unit previously failed to encode.
// sim_thr is the similarity threshold.
pub fn correct_unit_deletion(
    mut ds: DataSet,
    failed_trials: &mut [Vec<(usize, u64)>],
    sim_thr: f64,
) -> DataSet {
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
    assert_eq!(ds.encoded_reads.len(), failed_trials.len());
    let failed_trials: usize = ds
        .encoded_reads
        .par_iter_mut()
        .zip(failed_trials.par_iter_mut())
        .filter(|(r, _)| r.nodes.len() > 1)
        .map(|(r, failed_trial)| {
            let seq: Vec<_> = raw_seq[&r.id]
                .iter()
                .map(|x| x.to_ascii_uppercase())
                .collect();
            correct_deletion_error(
                r,
                failed_trial,
                selected_chunks,
                &read_skeltons,
                &seq,
                sim_thr,
            )
        })
        .sum();
    debug!("ENCODE\t{}", failed_trials);
    ds
}

/// Same as `correct unit deletion selected`. The only difference is that this function controls which read
/// would be corrected. Maybe boost the efficiency.
pub fn correct_unit_deletion_selected(ds: &mut DataSet, reads: &HashSet<u64>, sim_thr: f64) {
    // This skeltons contain all reads, as some of the unfocal reads
    // would be help to correct focal reads.
    let read_skeltons: Vec<_> = ds
        .encoded_reads
        .iter()
        .filter(|read| read.nodes.len() > 1)
        .map(|read| ReadSkelton::from_rich_nodes(&read.nodes))
        .collect();
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .filter(|r| reads.contains(&r.id))
        .map(|read| (read.id, read.seq()))
        .collect();
    let selected_chunks = &ds.selected_chunks;
    ds.encoded_reads
        .par_iter_mut()
        .filter(|r| r.nodes.len() > 1 && reads.contains(&r.id))
        .for_each(|r| {
            let seq: Vec<_> = raw_seq[&r.id]
                .iter()
                .map(|x| x.to_ascii_uppercase())
                .collect();
            correct_deletion_error(
                r,
                &mut Vec::new(),
                selected_chunks,
                &read_skeltons,
                &seq,
                sim_thr,
            );
        });
}

// Aligment offset. We align [s-offset..e+offset] region to the unit.
const OFFSET: usize = 300;
pub fn correct_deletion_error(
    read: &mut EncodedRead,
    failed_trials: &mut Vec<(usize, u64)>,
    units: &[Unit],
    reads: &[ReadSkelton],
    seq: &[u8],
    sim_thr: f64,
) -> usize {
    // Inserption counts.
    let pileups = get_pileup(read, reads);
    // let threshold = get_threshold(&pileups);
    let threshold = 3;
    let nodes = &read.nodes;
    let take_len = nodes.len();
    let mut inserts = vec![];
    let abs = |x: usize, y: usize| match x.cmp(&y) {
        std::cmp::Ordering::Greater => x - y,
        _ => y - x,
    };
    let mut failed_number = 0;
    for (idx, pileup) in pileups.iter().enumerate().take(take_len).skip(1) {
        let head_node = pileup
            .check_insertion_head(nodes, threshold, idx)
            .filter(|&(_, _, id)| !failed_trials.contains(&(idx, id)));
        if let Some((start_position, direction, uid)) = head_node {
            let unit = units.iter().find(|u| u.id == uid).unwrap();
            let end_position = (start_position + unit.seq().len() + 2 * OFFSET).min(seq.len());
            let is_the_same_encode = nodes[idx].unit == uid
                && abs(nodes[idx].position_from_start, start_position) < unit.seq().len();
            if start_position < end_position && !is_the_same_encode {
                match encode_node(seq, start_position, end_position, direction, unit, sim_thr) {
                    Some(head_node) => inserts.push((idx, head_node)),
                    None => {
                        failed_number += 1;
                        failed_trials.push((idx, uid));
                    }
                }
            }
        }
        let tail_node = pileup
            .check_insertion_tail(nodes, threshold, idx)
            .filter(|&(_, _, id)| !failed_trials.contains(&(idx, id)));
        if let Some((end_position, direction, uid)) = tail_node {
            let unit = units.iter().find(|u| u.id == uid).unwrap();
            let end_position = end_position.min(seq.len());
            let start_position = end_position.saturating_sub(unit.seq().len() + 2 * OFFSET);
            if start_position < end_position {
                match encode_node(seq, start_position, end_position, direction, unit, sim_thr) {
                    Some(tail_node) => inserts.push((idx, tail_node)),
                    None => {
                        failed_number += 1;
                        failed_trials.push((idx, uid));
                    }
                }
            }
        }
    }
    if !inserts.is_empty() {
        failed_trials.clear();
    }
    for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
        read.nodes.insert(idx + accum_inserts, node);
    }
    if !read.nodes.is_empty() {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use super::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| n.unit);
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
    failed_number
}

// 0.35 * unit.seq() distance is regarded as too far: corresponding 30% errors.
// pub const ALIGN_LIMIT: f64 = 0.35;
// Try to Encode Node. Return Some(node) if the alignment is good.
fn encode_node(
    seq: &[u8],
    start: usize,
    end: usize,
    is_forward: bool,
    unit: &Unit,
    sim_thr: f64,
) -> Option<Node> {
    let mut query = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    query.iter_mut().for_each(|x| x.make_ascii_uppercase());
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    // Note that unit.seq would be smaller than query! So the operations should be reversed.
    let unitseq = unit.seq();
    let alignment = edlib_sys::edlib_align(unitseq, &query, mode, task);
    let dist_thr = (unitseq.len() as f64 * sim_thr).floor() as u32;
    let locations = alignment.locations.unwrap();
    let (aln_start, aln_end) = locations[0];
    let seq = query[aln_start..=aln_end].to_vec();
    // Let's try to align the read once more.
    let band = (seq.len() / 10).max(20);
    let (_, ops) = kiley::bialignment::global_banded(unitseq, &seq, 2, -6, -5, -1, band);
    let aln_dist = {
        let (mut upos, mut spos) = (0, 0);
        ops.iter()
            .map(|op| match op {
                kiley::bialignment::Op::Del => {
                    upos += 1;
                    1
                }
                kiley::bialignment::Op::Ins => {
                    spos += 1;
                    1
                }
                kiley::bialignment::Op::Mat => {
                    spos += 1;
                    upos += 1;
                    (unitseq[upos - 1] != seq[spos - 1]) as u32
                }
            })
            .sum::<u32>()
    };
    let ops = super::compress_kiley_ops(&ops);
    let max_indel = ops
        .iter()
        .map(|&op| match op {
            Op::Ins(l) | Op::Del(l) => l,
            _ => 0,
        })
        .max()
        .unwrap();
    debug!("TRY\t{}\t{}\t{}", dist_thr, aln_dist, max_indel);
    if dist_thr < aln_dist {
        return None;
    };
    let position_from_start = if is_forward {
        start + aln_start
    } else {
        start + query.len() - aln_end - 1
    };
    assert!(seq.iter().all(|x| x.is_ascii_uppercase()));
    let has_large_indel = ops.iter().any(|&op| match op {
        Op::Ins(l) | Op::Del(l) => l < super::INDEL_THRESHOLD,
        _ => false,
    });
    if has_large_indel {
        return None;
    }
    let node = Node {
        position_from_start,
        unit: unit.id,
        cluster: 0,
        is_forward,
        seq: String::from_utf8(seq).unwrap(),
        cigar: ops,
    };
    Some(node)
}

// Align read skeltons to read, return the pileup sumamries.
fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    // let read_nodes: HashSet<_> = read.nodes.iter().map(|n| n.unit).collect();
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    for query in reads
        .iter()
        .filter(|q| read.nodes.iter().any(|n| q.sets.contains(&n.unit)))
    {
        if let Some((aln, is_forward)) = alignment(read, query) {
            let query = if is_forward {
                query.clone()
            } else {
                query.rev()
            };
            let (mut r_ptr, mut q_ptr) = (0, 0);
            for op in aln {
                match op {
                    Op::Ins(l) => {
                        if r_ptr < read.nodes.len() && 0 < r_ptr {
                            pileups[r_ptr].add_head(query.nodes[q_ptr].clone());
                            if 1 < l {
                                let last_insertion = q_ptr + l - 1;
                                pileups[r_ptr].add_tail(query.nodes[last_insertion].clone());
                            }
                        }
                        q_ptr += l;
                    }
                    Op::Del(l) => {
                        r_ptr += l;
                    }
                    Op::Match(l) => {
                        for pileup in pileups.iter_mut().skip(r_ptr).take(l) {
                            pileup.coverage += 1;
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
    let (f_score, f_ops) = pairwise_alignment_gotoh(&read, query);
    let f_match = get_match_units(&f_ops);
    let query = query.rev();
    let (r_score, r_ops) = pairwise_alignment_gotoh(&read, &query);
    let r_match = get_match_units(&r_ops);
    let score_thr = if read.nodes.len() == 2 { 1 } else { SCORE_THR };
    // Dovetails should be proper. Not Del->In or In -> Del transition should reside.
    if r_score <= f_score && MIN_MATCH <= f_match && score_thr <= f_score {
        is_proper(&f_ops).then(|| (f_ops, true))
    } else if f_score <= r_score && MIN_MATCH <= r_match && score_thr <= r_score {
        is_proper(&r_ops).then(|| (r_ops, false))
    } else {
        None
    }
}

// Return true if the alignment is proper dovetail.
fn is_proper(ops: &[Op]) -> bool {
    ops.windows(2)
        .all(|xs| !matches!(xs, &[Op::Ins(_), Op::Del(_)] | &[Op::Del(_), Op::Ins(_)]))
}

const MIN_ALN: i32 = -10000000;
fn score(x: &LightNode, y: &LightNode) -> i32 {
    if x.unit == y.unit && x.is_forward == y.is_forward {
        1
    } else {
        MIN_ALN
    }
}

// This should return overlapping alignment and its score.
#[allow(dead_code)]
fn pairwise_alignment(read: &ReadSkelton, query: &ReadSkelton) -> (i32, Vec<Op>) {
    let (read, query) = (&read.nodes, &query.nodes);
    // We do not allow any mismatch by restricting the mismatch score to read.len() + query.len(),
    // which is apparently an upperbound.
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
        ops.push(Op::Del(read.len() - r_pos));
    }
    if query.len() != q_pos {
        ops.push(Op::Ins(query.len() - q_pos));
    }
    while 0 < r_pos && 0 < q_pos {
        let current_dist = dp[r_pos][q_pos];
        let match_score = dp[r_pos - 1][q_pos - 1] + score(&read[r_pos - 1], &query[q_pos - 1]);
        if current_dist == match_score {
            ops.push(Op::Match(1));
            q_pos -= 1;
            r_pos -= 1;
        } else if current_dist == dp[r_pos - 1][q_pos] - 1 {
            ops.push(Op::Del(1));
            r_pos -= 1;
        } else {
            assert_eq!(current_dist, dp[r_pos][q_pos - 1] - 1);
            ops.push(Op::Ins(1));
            q_pos -= 1;
        }
    }
    if r_pos != 0 {
        ops.push(Op::Del(r_pos));
    }
    if q_pos != 0 {
        ops.push(Op::Ins(q_pos));
    }
    ops.reverse();
    (dist, compress_operations(ops))
}

fn pairwise_alignment_gotoh(read: &ReadSkelton, query: &ReadSkelton) -> (i32, Vec<Op>) {
    let (read, query) = (&read.nodes, &query.nodes);
    // Mat,Ins,Del
    let mut dp = vec![vec![vec![0; query.len() + 1]; read.len() + 1]; 3];
    // Initialize.
    for i in 0..read.len() + 1 {
        dp[0][i][0] = MIN_ALN;
        dp[1][i][0] = MIN_ALN;
    }
    for j in 0..query.len() + 1 {
        dp[0][0][j] = MIN_ALN;
        dp[2][0][j] = MIN_ALN;
    }
    dp[0][0][0] = 0;
    // Filling DP Table.
    for (i, x) in read.iter().enumerate() {
        for (j, y) in query.iter().enumerate() {
            let (i, j) = (i + 1, j + 1);
            dp[0][i][j] = dp[0][i - 1][j - 1]
                .max(dp[1][i - 1][j - 1])
                .max(dp[2][i - 1][j - 1])
                + score(x, y);
            dp[1][i][j] = (dp[0][i][j - 1] - 1).max(dp[1][i][j - 1]);
            dp[2][i][j] = (dp[0][i - 1][j] - 1).max(dp[2][i - 1][j]);
        }
    }
    let (mut state, mut r_pos, mut q_pos, dist) = (0..read.len() + 1)
        .map(|i| (i, query.len()))
        .chain((0..query.len() + 1).map(|j| (read.len(), j)))
        .flat_map(|(i, j)| vec![(0, i, j), (1, i, j), (2, i, j)])
        .map(|(s, i, j)| (s, i, j, dp[s][i][j]))
        .max_by_key(|x| x.3)
        .unwrap();
    let mut ops = vec![];
    if read.len() != r_pos {
        ops.push(Op::Del(read.len() - r_pos));
    }
    if query.len() != q_pos {
        ops.push(Op::Ins(query.len() - q_pos));
    }
    while 0 < r_pos && 0 < q_pos {
        let current_dist = dp[state][r_pos][q_pos];
        if state == 0 {
            let dist = current_dist - score(&read[r_pos - 1], &query[q_pos - 1]);
            state = match dist {
                x if x == dp[0][r_pos - 1][q_pos - 1] => 0,
                x if x == dp[1][r_pos - 1][q_pos - 1] => 1,
                x if x == dp[2][r_pos - 1][q_pos - 1] => 2,
                _ => panic!("{},{}", r_pos, q_pos),
            };
            ops.push(Op::Match(1));
            r_pos -= 1;
            q_pos -= 1;
        } else if state == 1 {
            state = match current_dist {
                x if x == dp[0][r_pos][q_pos - 1] - 1 => 0,
                x if x == dp[1][r_pos][q_pos - 1] => 1,
                _ => unreachable!(),
            };
            ops.push(Op::Ins(1));
            q_pos -= 1;
        } else {
            state = match current_dist {
                x if x == dp[0][r_pos - 1][q_pos] - 1 => 0,
                x if x == dp[2][r_pos - 1][q_pos] => 2,
                _ => unreachable!(),
            };
            ops.push(Op::Del(1));
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
// fn get_threshold(pileups: &[Pileup]) -> usize {
//     let totcov = pileups.iter().map(|p| p.coverage).sum::<usize>();
//     // We need at least 3 insertions to confirm.
//     (totcov / 3 / pileups.len()).max(3)
// }

#[derive(Debug, Clone)]
pub struct Pileup {
    // insertion at the beggining of this node
    head_inserted: Vec<LightNode>,
    // insertion at the last of this node
    tail_inserted: Vec<LightNode>,
    coverage: usize,
}

impl Pileup {
    // Return the maximum insertion from the same unit, the same direction.
    fn max_insertion_head(&self) -> Option<(usize, u64, bool)> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.head_inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        count
            .iter()
            .max_by_key(|x| x.1)
            .map(|(&(a, b), &y)| (y, a, b))
    }
    fn max_insertion_tail(&self) -> Option<(usize, u64, bool)> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.tail_inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        count
            .iter()
            .max_by_key(|x| x.1)
            .map(|(&(a, b), &y)| (y, a, b))
    }
    fn new() -> Self {
        Self {
            head_inserted: vec![],
            tail_inserted: vec![],
            coverage: 0,
        }
    }
    fn information_head(&self, unit: u64, is_forward: bool) -> (isize, isize) {
        let inserts = self
            .head_inserted
            .iter()
            .filter(|node| node.unit == unit && node.is_forward == is_forward);
        Self::summarize(inserts)
    }
    fn information_tail(&self, unit: u64, is_forward: bool) -> (isize, isize) {
        let inserts = self
            .tail_inserted
            .iter()
            .filter(|node| node.unit == unit && node.is_forward == is_forward);
        Self::summarize(inserts)
    }
    fn summarize<'a>(
        inserts: impl std::iter::Iterator<Item = &'a LightNode> + Clone,
    ) -> (isize, isize) {
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
                inserts.fold((0, 0), |(count, total), node| match node.after_offset {
                    Some(len) => (count + 1, total + len),
                    None => (count, total),
                });
            total / count
        };
        (prev_offset, after_offset)
    }
    fn add_head(&mut self, node: LightNode) {
        self.head_inserted.push(node);
    }
    fn add_tail(&mut self, node: LightNode) {
        self.tail_inserted.push(node);
    }
    fn check_insertion_head(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> Option<(usize, bool, u64)> {
        let (max_num, max_unit, max_dir) = self.max_insertion_head()?;
        (threshold <= max_num).then(|| {
            let (uid, direction) = (max_unit, max_dir);
            let (prev_offset, _) = self.information_head(max_unit, max_dir);
            let start_position =
                (nodes[idx - 1].position_from_start + nodes[idx - 1].query_length()) as isize;
            let start_position = start_position + prev_offset;
            let start_position = (start_position as usize).saturating_sub(OFFSET);
            (start_position, direction, uid)
        })
    }
    fn check_insertion_tail(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> Option<(usize, bool, u64)> {
        let end_position = nodes.get(idx)?.position_from_start as isize;
        let (max_num, max_unit, max_dir) = self.max_insertion_tail()?;
        (threshold <= max_num).then(|| {
            let (uid, direction) = (max_unit, max_dir);
            let (_, after_offset) = self.information_tail(max_unit, max_dir);
            let end_position = (end_position + after_offset) as usize + OFFSET;
            // TODO: Validate end position.
            (end_position, direction, uid)
        })
    }
}

#[derive(Clone)]
pub struct ReadSkelton {
    nodes: Vec<LightNode>,
    sets: HashSet<u64>,
}

impl std::fmt::Debug for ReadSkelton {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for node in self.nodes.iter() {
            write!(f, "{:?}:", node)?;
        }
        Ok(())
    }
}

impl ReadSkelton {
    pub fn from_rich_nodes(nodes: &[Node]) -> Self {
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
        let sets: HashSet<_> = nodes.iter().map(|n| n.unit).collect();
        ReadSkelton { nodes, sets }
    }
    fn rev(&self) -> Self {
        let nodes: Vec<_> = self.nodes.iter().rev().map(LightNode::rev).collect();
        let sets = self.sets.clone();
        Self { nodes, sets }
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

#[cfg(test)]
mod deletion_fill {
    use super::*;
    #[test]
    fn aln_test() {
        let nodes: Vec<_> = vec![69, 148, 318, 0]
            .into_iter()
            .zip(vec![false, false, true, true])
            .map(|(unit, is_forward)| LightNode {
                prev_offset: None,
                unit,
                is_forward,
                after_offset: None,
            })
            .collect();
        let read = ReadSkelton { nodes };
        let nodes: Vec<_> = vec![69, 221, 286, 148, 318]
            .into_iter()
            .zip(vec![false, true, true, false, true])
            .map(|(unit, is_forward)| LightNode {
                prev_offset: None,
                unit,
                is_forward,
                after_offset: None,
            })
            .collect();
        let query = ReadSkelton { nodes };
        let (score, ops) = pairwise_alignment(&read, &query);
        assert_eq!(score, 1, "{:?}", ops);
    }
    #[test]
    fn aln_test_gotoh() {
        let into_reads = |nodes: Vec<u64>| {
            let nodes: Vec<_> = nodes
                .into_iter()
                .map(|unit| LightNode {
                    prev_offset: None,
                    unit,
                    is_forward: true,
                    after_offset: None,
                })
                .collect();
            ReadSkelton { nodes }
        };
        let read = into_reads(vec![69, 148, 318, 0]);
        let query = into_reads(vec![69, 221, 286, 148, 318]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 2, "{:?}", ops);
        let read = into_reads(vec![0]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 1, "{:?}", ops);
        let read = into_reads(vec![0, 4]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 1, "{:?}", ops);
        let (score, ops) = pairwise_alignment_gotoh(&query, &read);
        assert_eq!(score, 1, "{:?}", ops);
        let read = into_reads(vec![0, 1]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 2, "{:?}", ops);
        let (score, ops) = pairwise_alignment_gotoh(&query, &read);
        assert_eq!(score, 2, "{:?}", ops);
    }
}
