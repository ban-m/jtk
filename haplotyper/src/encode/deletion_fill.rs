//! Filling deletion.
use definitions::*;
use rayon::prelude::*;
// use log::*;
use std::collections::{HashMap, HashSet};
/// The second argument is the vector of (index,unit_id) of the previous failed trials.
/// for example, if failed_trials[i][0] = (j,id), then, we've already tried to encode the id-th unit after the j-th
/// position of the i-th read, and failed it.
/// If we can encode some position in the i-th read, the failed trials would be erased, as it change the
/// condition of the read, making it possible to encode an unit previously failed to encode.
/// sim_thr is the similarity threshold.
/// This function corrects "unit-deletions". To do that,
/// it first align other reads onto a read to be corrected in unit resolution, detecting putative insertions.
/// Then, it tries to encode these putative insertions in base-pair resolution.
/// Note that, in the first - unit resolution - alignment, there's no distinction between clusters.
/// However, in the second alignment, it tries to encode the putative region by each cluster's representative.
/// Of course, if there's only one cluster for a unit, then, it just tries to encode by that unit.
pub fn correct_unit_deletion(mut ds: DataSet, sim_thr: f64) -> DataSet {
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.seq()))
        .collect();
    let mut current: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
    'outer: for _ in 0..3 {
        let representative = take_consensus_sequence(&ds);
        // i->vector of failed index and units.
        let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
        for i in 0..15 {
            let read_skeltons: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
            ds.encoded_reads
                .par_iter_mut()
                .zip(failed_trials.par_iter_mut())
                .filter(|(r, _)| r.nodes.len() > 1)
                .for_each(|(r, fails)| {
                    let mut seq: Vec<_> = raw_seq[&r.id].to_vec();
                    seq.iter_mut().for_each(u8::make_ascii_uppercase);
                    correct_deletion_error(
                        (r, &seq),
                        fails,
                        &representative,
                        &read_skeltons,
                        sim_thr,
                    );
                });
            let after: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
            debug!("Filled:{}\t{}", current, after);
            if after <= current {
                if i == 0 {
                    debug!("Filled\tBREAK\tOuter");
                    break 'outer;
                } else {
                    debug!("Filled\tBREAK\tInner");
                    break;
                }
            }
            current = after;
        }
    }
    ds
}

// /// Same as `correct unit deletion selected`. The only difference is that this function controls which read
// /// would be corrected. Maybe boost the efficiency.
// pub fn correct_unit_deletion_selected(ds: &mut DataSet, reads: &HashSet<u64>, sim_thr: f64) {
//     // This skeltons contain all reads, as some of the unfocal reads
//     // would be help to correct focal reads.
//     let read_skeltons: Vec<_> = ds
//         .encoded_reads
//         .iter()
//         .filter(|read| read.nodes.len() > 1)
//         .map(|read| ReadSkelton::from_rich_nodes(&read.nodes))
//         .collect();
//     let raw_seq: HashMap<_, _> = ds
//         .raw_reads
//         .iter()
//         .filter(|r| reads.contains(&r.id))
//         .map(|read| (read.id, read.seq()))
//         .collect();
//     let selected_chunks = &ds.selected_chunks;
//     ds.encoded_reads
//         .par_iter_mut()
//         .filter(|r| r.nodes.len() > 1 && reads.contains(&r.id))
//         .for_each(|r| {
//             correct_deletion_error(
//                 r,
//                 &mut Vec::new(),
//                 selected_chunks,
//                 &read_skeltons,
//                 raw_seq[&r.id],
//                 sim_thr,
//             );
//         });
// }

// Take consensus of each cluster of each unit, return the consensus seuqneces.
// UnitID->(clsuterID, its consensus).
#[allow(dead_code)]
fn take_consensus_sequence(ds: &DataSet) -> HashMap<u64, Vec<(u64, Vec<u8>)>> {
    fn polish(xs: &[&[u8]], unit: &Unit) -> Vec<u8> {
        // use kiley::gphmm::*;
        // let mut hmm = kiley::gphmm::GPHMM::<Cond>::clr();
        // let band_width = 100;
        // TODO:Maybe we need to polish the unit with GPHMM, why not?
        // Maybe the computational cost matters....?
        kiley::bialignment::polish_until_converge_banded(unit.seq(), xs, 100)
    }
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u)).collect();
    let mut bucket: HashMap<u64, Vec<_>> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        bucket
            .entry(node.unit)
            .or_default()
            .push((node.cluster, node.seq()));
    }
    bucket
        .par_iter()
        .filter(|&(_, xs)| !xs.is_empty())
        .map(|(&unit, bucket)| {
            let ref_unit = &ref_units[&unit];
            let mut clusters: HashMap<_, Vec<&[u8]>> = HashMap::new();
            for (cl, seq) in bucket {
                clusters.entry(*cl).or_default().push(seq);
            }
            assert!(!clusters.is_empty());
            let representative: Vec<_> = if clusters.len() == 1 {
                clusters
                    .iter()
                    .map(|(&cl, _)| (cl, ref_unit.seq().to_vec()))
                    .collect()
            } else {
                clusters
                    .iter()
                    .filter(|(_, xs)| !xs.is_empty())
                    .map(|(&cl, xs)| (cl, polish(xs, ref_unit)))
                    .collect()
            };
            (unit, representative)
        })
        .collect()
}

#[inline]
fn abs(x: usize, y: usize) -> usize {
    x.max(y) - x.min(y)
}

// Aligment offset. We align [s-offset..e+offset] region to the unit.
const OFFSET: usize = 300;
pub fn correct_deletion_error(
    (read, seq): (&mut EncodedRead, &[u8]),
    failed_trials: &mut Vec<(usize, u64)>,
    units: &HashMap<u64, Vec<(u64, Vec<u8>)>>,
    reads: &[ReadSkelton],
    sim_thr: f64,
) {
    // Inserption counts.
    let pileups = get_pileup(read, reads);
    // TODO:Parametrize here.
    let threshold = 3;
    let nodes = &read.nodes;
    let take_len = nodes.len();
    let mut inserts = vec![];
    let seq: Vec<_> = seq.iter().map(|x| x.to_ascii_uppercase()).collect();
    let seq = &seq;
    for (idx, pileup) in pileups.iter().enumerate().take(take_len).skip(1) {
        let head_cand = pileup.check_insertion_head(nodes, threshold, idx);
        let head_best = head_cand
            .iter()
            .filter(|&&(_, _, id)| !failed_trials.contains(&(idx, id)))
            .filter_map(|&(start_position, direction, uid)| {
                units
                    .get(&uid)?
                    .iter()
                    .filter_map(|&(cluster, ref unit)| {
                        let end_position =
                            (start_position + unit.len() + 2 * OFFSET).min(seq.len());
                        let is_the_same_encode = nodes[idx].unit == uid
                            && abs(nodes[idx].position_from_start, start_position) < unit.len();
                        if start_position < end_position && !is_the_same_encode {
                            let position = (start_position, end_position, direction);
                            let unit_info = (unit.as_slice(), uid, cluster);
                            encode_node(seq, position, unit_info, sim_thr)
                        } else {
                            None
                        }
                    })
                    .max_by_key(|x| x.1)
            })
            .max_by_key(|x| x.1);
        match head_best {
            Some((head_node, _)) => inserts.push((idx, head_node)),
            None => failed_trials.extend(head_cand.iter().map(|&(_, _, uid)| (idx, uid))),
        }
        let tail_cand = pileup.check_insertion_tail(nodes, threshold, idx);
        let tail_best = tail_cand
            .iter()
            .filter(|&&(_, _, id)| !failed_trials.contains(&(idx, id)))
            .filter_map(|&(end_position, direction, uid)| {
                units
                    .get(&uid)?
                    .iter()
                    .filter_map(|&(cluster, ref unit)| {
                        let end_position = end_position.min(seq.len());
                        let start_position = end_position.saturating_sub(unit.len() + 2 * OFFSET);
                        if start_position < end_position {
                            let positions = (start_position, end_position, direction);
                            let unit_info = (unit.as_slice(), uid, cluster);
                            encode_node(seq, positions, unit_info, sim_thr)
                        } else {
                            None
                        }
                    })
                    .max_by_key(|x| x.1)
            })
            .max_by_key(|x| x.1);
        match tail_best {
            Some((tail_node, _)) => inserts.push((idx, tail_node)),
            None => failed_trials.extend(tail_cand.iter().map(|&(_, _, uid)| (idx, uid))),
        }
    }
    if !inserts.is_empty() {
        failed_trials.clear();
        for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
            read.nodes.insert(idx + accum_inserts, node);
        }
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use super::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| (n.unit, n.position_from_start));
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
}

// Try to Encode Node. Return Some(node) if the alignment is good.
// Return also the alignment score of the encoding.
// The match score is 2, mism is -6, gap open is -5, and gap ext is -1.
fn encode_node(
    seq: &[u8],
    (start, end, is_forward): (usize, usize, bool),
    (unitseq, uid, cluster): (&[u8], u64, u64),
    sim_thr: f64,
) -> Option<(Node, i32)> {
    let mut query = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    query.iter_mut().for_each(|x| x.make_ascii_uppercase());
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    // Note that unit.seq would be smaller than query! So the operations should be reversed.
    // let unitseq = unit.seq();
    let alignment = edlib_sys::edlib_align(unitseq, &query, mode, task);
    let dist_thr = (unitseq.len() as f64 * sim_thr).floor() as u32;
    let locations = alignment.locations.unwrap();
    let (aln_start, aln_end) = locations[0];
    let seq = query[aln_start..=aln_end].to_vec();
    let band = (seq.len() / 10).max(20);
    let (score, ops) = kiley::bialignment::global_banded(unitseq, &seq, 2, -6, -5, -1, band);
    let aln_dist = ops
        .iter()
        .filter(|&&op| op != kiley::bialignment::Op::Mat)
        .count() as u32;
    // Mat=>-1, Other->1
    let indel_mism = ops
        .iter()
        .map(|&op| 1 - 2 * (op == kiley::bialignment::Op::Mat) as i32);
    let max_indel = super::max_region(indel_mism).max(0) as usize;
    let has_large_indel = max_indel > super::INDEL_THRESHOLD;
    let position_from_start = if is_forward {
        start + aln_start
    } else {
        start + query.len() - aln_end - 1
    };
    assert!(seq.iter().all(|x| x.is_ascii_uppercase()));
    if dist_thr < aln_dist || has_large_indel {
        if log_enabled!(log::Level::Trace) {
            trace!("{}\t{}\t{}\t{}\tNG", uid, cluster, max_indel, aln_dist);
            let (xr, ar, yr) = kiley::bialignment::recover(unitseq, &seq, &ops);
            for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
                eprintln!("{}", String::from_utf8_lossy(xr));
                eprintln!("{}", String::from_utf8_lossy(ar));
                eprintln!("{}\n", String::from_utf8_lossy(yr));
            }
        }
        return None;
    }
    trace!("{}\t{}\t{}\t{}\tOK", uid, cluster, max_indel, aln_dist);
    let ops = super::compress_kiley_ops(&ops);
    let node = Node {
        position_from_start,
        unit: uid,
        cluster,
        is_forward,
        seq: String::from_utf8(seq).unwrap(),
        cigar: ops,
    };
    Some((node, score))
}

// Align read skeltons to read, return the pileup sumamries.
fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    // let read_nodes: HashSet<_> = read.nodes.iter().map(|n| n.unit).collect();
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    for query in reads
        .iter()
        .filter(|q| read.nodes.iter().any(|n| q.sets.contains(&n.unit)))
    {
        let (aln, is_forward) = match alignment(read, query) {
            Some(res) => res,
            None => continue,
        };
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
    fn insertion_head(&self) -> impl std::iter::Iterator<Item = (usize, u64, bool)> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.head_inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        count.into_iter().map(|((a, b), y)| (y, a, b))
    }
    fn insertion_tail(&self) -> impl std::iter::Iterator<Item = (usize, u64, bool)> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.tail_inserted.iter() {
            *count.entry((node.unit, node.is_forward)).or_default() += 1;
        }
        count.into_iter().map(|((a, b), y)| (y, a, b))
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
    ) -> Vec<(usize, bool, u64)> {
        //Option<(usize, bool, u64)> {
        //let (max_num, max_unit, max_dir) = self.max_insertion_head()?;
        self.insertion_head()
            .filter(|&(num, _, _)| threshold <= num)
            .map(|(_, uid, direction)| {
                let (prev_offset, _) = self.information_head(uid, direction);
                let start_position =
                    (nodes[idx - 1].position_from_start + nodes[idx - 1].query_length()) as isize;
                let start_position = start_position + prev_offset;
                let start_position = (start_position as usize).saturating_sub(OFFSET);
                (start_position, direction, uid)
            })
            .collect()

        // (threshold <= max_num).then(|| {
        //     let (uid, direction) = (max_unit, max_dir);
        //     let (prev_offset, _) = self.information_head(max_unit, max_dir);
        //     let start_position =
        //         (nodes[idx - 1].position_from_start + nodes[idx - 1].query_length()) as isize;
        //     let start_position = start_position + prev_offset;
        //     let start_position = (start_position as usize).saturating_sub(OFFSET);
        //     (start_position, direction, uid)
        // })
    }
    fn check_insertion_tail(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> Vec<(usize, bool, u64)> {
        let end_position = match nodes.get(idx) {
            Some(res) => res.position_from_start as isize,
            None => return Vec::new(),
        };
        self.insertion_tail()
            .filter(|&(num, _, _)| threshold <= num)
            .map(|(_, uid, direction)| {
                let (_, after_offset) = self.information_tail(uid, direction);
                let end_position = (end_position + after_offset) as usize + OFFSET;
                // TODO: Validate end position.
                (end_position, direction, uid)
            })
            .collect()
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
    pub fn new(read: &EncodedRead) -> Self {
        Self::from_rich_nodes(&read.nodes)
    }
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
        let sets: HashSet<_> = nodes.iter().map(|x| x.unit).collect();
        let read = ReadSkelton { nodes, sets };
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
        let sets: HashSet<_> = nodes.iter().map(|x| x.unit).collect();
        let query = ReadSkelton { nodes, sets };
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
            let sets: HashSet<_> = nodes.iter().map(|x| x.unit).collect();
            ReadSkelton { nodes, sets }
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
