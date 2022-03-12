//! Filling deletion.
// TODO:Curretly the alignment can not handle to re-encode missing tips. It reduce some of the reads which should be encoded otherwise...
use crate::ALN_PARAMETER;
use definitions::*;
use rayon::prelude::*;
// use log::*;
use std::collections::{HashMap, HashSet};
// identity would be increased by this value when evaluating the edges.
const EDGE_BONUS: f64 = 0.05;
// Evaluating length of each side.
const EDGE_LEN: usize = 500;
#[derive(Debug, Clone)]
pub struct CorrectDeletionConfig {
    re_clustering: bool,
    sim_thr: f64,
}

impl CorrectDeletionConfig {
    pub fn new(re_clustering: bool, sim_thr: f64) -> Self {
        Self {
            re_clustering,
            sim_thr,
        }
    }
}

pub trait CorrectDeletion {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig);
}

impl CorrectDeletion for DataSet {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig) {
        let find_new_units = correct_unit_deletion(self, config.sim_thr);
        if config.re_clustering {
            // Log original assignments.
            let original_assignments = log_original_assignments(self);
            // Log previous copy number
            let prev_copy_numbers: HashMap<u64, _> = self
                .selected_chunks
                .iter()
                .map(|c| (c.id, c.copy_num))
                .collect();
            // Re-estimate copy number
            // -> Erase cluster
            self.encoded_reads
                .iter_mut()
                .flat_map(|r| r.nodes.iter_mut())
                .for_each(|n| n.cluster = 0);
            use crate::multiplicity_estimation::*;
            let seed = (231043290490.0 * config.sim_thr).round() as u64;
            let config = MultiplicityEstimationConfig::new(1, seed, None);
            self.estimate_multiplicity(&config);
            // Retain all the units changed their copy numbers.
            let chainged_units: HashSet<_> = self
                .selected_chunks
                .iter()
                .filter_map(|c| match prev_copy_numbers.get(&c.id) {
                    None => None,
                    Some(&prev) if c.copy_num == prev => None,
                    Some(_) => Some(c.id),
                })
                .collect();
            // Merge these two.
            let selection: HashSet<_> = find_new_units.union(&chainged_units).copied().collect();
            // Recover the original assignments on the retained units.
            self.encoded_reads
                .iter_mut()
                .zip(original_assignments)
                .for_each(|(read, (id, log))| {
                    assert_eq!(read.id, id);
                    recover_original_assignments(read, &log, &selection);
                });
            // Reclustering.
            crate::local_clustering::local_clustering_selected(self, &selection);
        }
    }
}
// Logging the original assignment into a vector.
fn log_original_assignments(ds: &DataSet) -> Vec<(u64, Vec<(u64, u64)>)> {
    ds.encoded_reads
        .iter()
        .map(|r| {
            let xs: Vec<_> = r.nodes.iter().map(|u| (u.unit, u.cluster)).collect();
            (r.id, xs)
        })
        .collect()
}

// Recover the previous clustering. Note that sometimes the node is added
// so that the length of the read is different from the logged one.
// But it does not removed!
fn recover_original_assignments(read: &mut EncodedRead, log: &[(u64, u64)], except: &HashSet<u64>) {
    let mut read = read.nodes.iter_mut();
    for &(unit, cluster) in log {
        for node in &mut read {
            if node.unit == unit {
                if !except.contains(&node.unit) {
                    node.cluster = cluster;
                }
                break;
            }
        }
    }
}

// pub fn correct_unit_deletion(ds: &mut DataSet, sim_thr: f64) -> HashSet<u64> {
//     let raw_seq: HashMap<_, _> = ds
//         .raw_reads
//         .iter()
//         .map(|read| (read.id, read.seq()))
//         .collect();
//     let mut current: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
//     const OUTER_LOOP: usize = 3;
//     const INNER_LOOP: usize = 15;
//     let mut find_new_node = HashSet::new();
//     'outer: for t in 0..OUTER_LOOP {
//         let representative = take_consensus_sequence(ds);
//         let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
//         // i->vector of failed index and units.
//         let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
//         for i in 0..INNER_LOOP {
//             let read_skeltons: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
//             let failed_trials = failed_trials.par_iter_mut();
//             let reads = ds.encoded_reads.par_iter_mut().zip(failed_trials);
//             let filtered_reads = reads.filter(|(r, _)| r.nodes.len() > 1);
//             let newly_encoded_units: Vec<_> = filtered_reads
//                 .flat_map(|(r, fails)| {
//                     let mut seq: Vec<_> = raw_seq[&r.id].to_vec();
//                     seq.iter_mut().for_each(u8::make_ascii_uppercase);
//                     let r_sk = &read_skeltons;
//                     let r = (r, seq.as_slice());
//                     correct_deletion_error(r, fails, &representative, &units, &r_sk, sim_thr)
//                 })
//                 .collect();
//             find_new_node.extend(newly_encoded_units);
//             let after: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
//             debug!("Filled:{}\t{}", current, after);
//             if after <= current && i == 0 {
//                 debug!("Filled\tBREAK\tOuter\t{}", t);
//                 break 'outer;
//             } else if after <= current {
//                 debug!("Filled\tBREAK\tInner");
//                 break;
//             }
//             current = after;
//         }
//     }
//     find_new_node
// }

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
/// Auto-tune the similarity threshold.
pub fn correct_unit_deletion(ds: &mut DataSet, fallback: f64) -> HashSet<u64> {
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.seq()))
        .collect();
    let mut current: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
    let read_error_rate: Vec<_> = estimate_upper_error_rate(ds, fallback);
    const OUTER_LOOP: usize = 3;
    const INNER_LOOP: usize = 15;
    let mut find_new_node = HashSet::new();
    'outer: for t in 0..OUTER_LOOP {
        let representative = take_consensus_sequence(ds);
        // i->vector of failed index and units.
        let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
        let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
        for i in 0..INNER_LOOP {
            let read_skeltons: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
            let failed_trials = failed_trials.par_iter_mut();
            let read_error_rate = read_error_rate.par_iter();
            let reads = ds
                .encoded_reads
                .par_iter_mut()
                .zip(failed_trials)
                .zip(read_error_rate);
            let filtered_reads = reads.filter(|((r, _), _)| r.nodes.len() > 1);
            let newly_encoded_units: Vec<_> = filtered_reads
                .flat_map(|((r, fails), &sim_thr)| {
                    let mut seq: Vec<_> = raw_seq[&r.id].to_vec();
                    seq.iter_mut().for_each(u8::make_ascii_uppercase);
                    let r_sk = &read_skeltons;
                    let r = (r, seq.as_slice());
                    correct_deletion_error(r, fails, &representative, &units, r_sk, sim_thr)
                })
                .collect();
            find_new_node.extend(newly_encoded_units);
            let after: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
            debug!("Filled:{}\t{}", current, after);
            if after <= current && i == 0 {
                debug!("Filled\tBREAK\tOuter\t{}", t);
                break 'outer;
            } else if after <= current {
                debug!("Filled\tBREAK\tInner");
                break;
            }
            current = after;
        }
    }
    find_new_node
}

fn estimate_upper_error_rate(ds: &DataSet, fallback: f64) -> Vec<f64> {
    let ref_chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let error_rates: Vec<Vec<_>> = ds
        .encoded_reads
        .par_iter()
        .map(|r| {
            r.nodes
                .iter()
                .map(|n| {
                    let (_, aln, _) = n.recover(ref_chunks[&n.unit]);
                    let error = aln.iter().filter(|&&b| b != b'|').count();
                    error as f64 / aln.len() as f64
                })
                .collect()
        })
        .collect();
    let (sd_sum, sd_num) = error_rates
        .iter()
        .filter(|errs| errs.len() > 2)
        .map(|errs| {
            let (sum, sumsq) = errs
                .iter()
                .fold((0f64, 0f64), |(sum, sumsq), x| (sum + x, sumsq + x * x));
            let mean = sum / errs.len() as f64;
            let var = sumsq / errs.len() as f64 - mean * mean;
            assert!(0f64 <= var);
            var.sqrt()
        })
        .fold((0f64, 0), |(sdsum, num), sd| (sdsum + sd, num + 1));
    let sd_mean = sd_sum / sd_num as f64;
    debug!("MEAN of SD\t{:.3}", sd_mean);
    error_rates
        .into_par_iter()
        .map(|mut errors| match errors.len() {
            0..=2 => fallback,
            x => {
                let median = errors
                    .select_nth_unstable_by(x / 2, |x, y| x.partial_cmp(y).unwrap())
                    .1;
                *median + 4f64 * sd_mean
            }
        })
        .collect()
}

// Take consensus of each cluster of each unit, return the consensus seuqneces.
// UnitID->(clsuterID, its consensus).
fn take_consensus_sequence(ds: &DataSet) -> HashMap<u64, Vec<(u64, Vec<u8>)>> {
    fn polish(xs: &[&[u8]], unit: &Unit, band: usize) -> Vec<u8> {
        kiley::bialignment::guided::polish_until_converge(unit.seq(), xs, band)
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
            let band = ds.read_type.band_width(ref_unit.seq().len());
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
                    .map(|(&cl, xs)| (cl, polish(xs, ref_unit, band)))
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
// returns the ids of the units newly encoded.
pub fn correct_deletion_error(
    (read, seq): (&mut EncodedRead, &[u8]),
    failed_trials: &mut Vec<(usize, u64)>,
    consensi: &HashMap<u64, Vec<(u64, Vec<u8>)>>,
    units: &HashMap<u64, &Unit>,
    reads: &[ReadSkelton],
    sim_thr: f64,
) -> Vec<u64> {
    let pileups = get_pileup(read, reads);
    // TODO:Parametrize here.
    let threshold = 3;
    let nodes = &read.nodes;
    let mut inserts = vec![];
    assert!(seq.iter().all(|x| x.is_ascii_uppercase()));
    for (idx, pileup) in pileups.iter().enumerate().skip(1) {
        let mut head_cand = pileup.check_insertion_head(nodes, threshold, idx);
        head_cand.retain(|&(_, _, uid)| !failed_trials.contains(&(idx, uid)));
        let head_best = try_encoding_head(nodes, &head_cand, idx, consensi, units, seq, sim_thr);
        match head_best {
            Some((head_node, _)) => inserts.push((idx, head_node)),
            None => failed_trials.extend(head_cand.into_iter().map(|x| (idx, x.2))),
        }
        let mut tail_cand = pileup.check_insertion_tail(nodes, threshold, idx);
        tail_cand.retain(|&(_, _, uid)| !failed_trials.contains(&(idx, uid)));
        let tail_best = try_encoding_tail(nodes, &tail_cand, idx, consensi, units, seq, sim_thr);
        match tail_best {
            Some((tail_node, _)) => inserts.push((idx, tail_node)),
            None => failed_trials.extend(tail_cand.into_iter().map(|x| (idx, x.2))),
        }
    }
    let new_inserts: Vec<_> = inserts.iter().map(|(_, n)| n.unit).collect();
    if !inserts.is_empty() {
        failed_trials.clear();
        for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
            read.nodes.insert(idx + accum_inserts, node);
        }
        let mut nodes = Vec::with_capacity(read.nodes.len());
        nodes.append(&mut read.nodes);
        use super::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| (n.unit, n.position_from_start));
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
    new_inserts
}

fn try_encoding_head(
    nodes: &[Node],
    head_cand: &[(usize, bool, u64)],
    idx: usize,
    consensi: &HashMap<u64, Vec<(u64, Vec<u8>)>>,
    units: &HashMap<u64, &Unit>,
    seq: &[u8],
    sim_thr: f64,
) -> Option<(Node, i32)> {
    head_cand
        .iter()
        .filter(|&&(start_position, _, _)| start_position < seq.len())
        .filter_map(|&(start_position, direction, uid)| {
            let unit = *units.get(&uid)?;
            let consensi = consensi.get(&uid)?.iter();
            consensi
                .filter_map(|&(cluster, ref cons)| {
                    let end_position = (start_position + cons.len() + 2 * OFFSET).min(seq.len());
                    let is_the_same_encode = match nodes.get(idx) {
                        Some(node) => {
                            node.unit == uid
                                && abs(node.position_from_start, start_position) < cons.len()
                        }
                        None => false,
                    };
                    assert!(start_position < end_position);
                    if is_the_same_encode {
                        return None;
                    }
                    let position = (start_position, end_position, direction);
                    let unit_info = (unit, cluster, cons.as_slice());
                    encode_node(seq, position, unit_info, sim_thr)
                })
                .max_by_key(|x| x.1)
        })
        .max_by_key(|x| x.1)
}

fn try_encoding_tail(
    nodes: &[Node],
    tail_cand: &[(usize, bool, u64)],
    idx: usize,
    consensi: &HashMap<u64, Vec<(u64, Vec<u8>)>>,
    units: &HashMap<u64, &Unit>,
    seq: &[u8],
    sim_thr: f64,
) -> Option<(Node, i32)> {
    tail_cand
        .iter()
        .filter_map(|&(end_position, direction, uid)| {
            let unit = *units.get(&uid)?;
            consensi
                .get(&uid)?
                .iter()
                .filter_map(|&(cluster, ref cons)| {
                    let end_position = end_position.min(seq.len());
                    let start_position = end_position.saturating_sub(cons.len() + 2 * OFFSET);
                    assert!(start_position < end_position);
                    let is_the_same_encode = match nodes.get(idx) {
                        Some(node) => {
                            node.unit == uid
                                && abs(node.position_from_start, start_position) < cons.len()
                        }
                        None => false,
                    };
                    if is_the_same_encode {
                        return None;
                    }
                    assert!(start_position < end_position);
                    let positions = (start_position, end_position, direction);
                    let unit_info = (unit, cluster, cons.as_slice());
                    encode_node(seq, positions, unit_info, sim_thr)
                })
                .max_by_key(|x| x.1)
        })
        .max_by_key(|x| x.1)
}

// Try to Encode Node. Return Some(node) if the alignment is good.
// Return also the alignment score of the encoding.
// The match score is 2, mism is -6, gap open is -5, and gap ext is -1.
fn encode_node(
    query: &[u8],
    (start, end, is_forward): (usize, usize, bool),
    (unit, cluster, unitseq): (&Unit, u64, &[u8]),
    sim_thr: f64,
) -> Option<(Node, i32)> {
    // Initial filter.
    // If the query is shorter than the unitseq,
    // at least we need the edit operations to fill the gaps.
    // This is lower bound of the sequence identity.
    let edit_dist_lower_bound = unitseq.len().saturating_sub(end - start);
    let diff_lower_bound = edit_dist_lower_bound as f64 / unitseq.len() as f64;
    if sim_thr < diff_lower_bound {
        return None;
    }
    // Tune the query...
    let mut query = if is_forward {
        query[start..end].to_vec()
    } else {
        bio_utils::revcmp(&query[start..end])
    };
    query.iter_mut().for_each(u8::make_ascii_uppercase);
    let (seq, aln_start, aln_end, kops, score) =
        fine_mapping(&query, (unit, cluster, unitseq), sim_thr)?;
    let ops = super::compress_kiley_ops(&kops);
    let cl = unit.cluster_num;
    let position_from_start = match is_forward {
        true => start + aln_start,
        false => start + query.len() - aln_end,
    };
    // I think we should NOT make likelihood gain to some biased value,
    // as 1. if the alignment gives the certaintly, then we can impute the clustering by the alignment,
    // 2. if `cluster` assignment is just by chance,
    // then we just should not introduce any bias into the likelihood gain.
    let mut node = Node::new(unit.id, is_forward, seq, ops, position_from_start, cl);
    node.cluster = cluster;
    Some((node, score))
}

type FineMapping<'a> = (&'a [u8], usize, usize, Vec<kiley::Op>, i32);
fn fine_mapping<'a>(
    query: &'a [u8],
    (unit, cluster, unitseq): (&Unit, u64, &[u8]),
    sim_thr: f64,
) -> Option<FineMapping<'a>> {
    //use kiley::bialignment::guided::global_guided;
    use kiley::bialignment::guided::infix_guided;
    fn edlib_op_to_kiley_op(ops: &[u8]) -> Vec<kiley::Op> {
        use kiley::Op::*;
        ops.iter()
            .map(|&op| [Match, Ins, Del, Mismatch][op as usize])
            .collect()
    }
    let origlen = query.len();
    let (query, aln_start, aln_end, ops, band) = {
        let mode = edlib_sys::AlignMode::Infix;
        let task = edlib_sys::AlignTask::Alignment;
        // Note that unit.seq would be smaller than query! So the operations should be reversed.
        let alignment = edlib_sys::edlib_align(unitseq, query, mode, task);
        let (aln_start, aln_end) = alignment.locations.unwrap()[0];
        let aln_end = aln_end + 1;
        let band = (((aln_end - aln_start + 1) as f64 * sim_thr * 0.3).ceil() as usize).max(10);
        let mut ops = Vec::with_capacity(3 * query.len() / 2);
        ops.extend(std::iter::repeat(kiley::Op::Del).take(aln_start));
        ops.extend(edlib_op_to_kiley_op(&alignment.operations.unwrap()));
        ops.extend(std::iter::repeat(kiley::Op::Del).take(query.len() - aln_end));
        let (_, mut ops) = infix_guided(query, unitseq, &ops, band, ALN_PARAMETER);
        // Reverse ops
        for op in ops.iter_mut() {
            *op = match *op {
                kiley::Op::Ins => kiley::Op::Del,
                kiley::Op::Del => kiley::Op::Ins,
                x => x,
            }
        }
        let (start, end) = trim_head_tail_insertion(&mut ops);
        let query = &query[start..query.len() - end];
        (query, start, query.len() - end, ops, band)
    };
    let (below_dissim, info) = {
        let mat_num = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
        let identity = mat_num as f64 / ops.len() as f64;
        let (head_identity, tail_identity) = edge_identity(unitseq, query, &ops, EDGE_LEN);
        let (rlen, qlen) = (unitseq.len(), query.len());
        let id = unit.id;
        let info = format!(
            "{}\t{}\t{:.2}\t{}\t{}\t{}",
            id, cluster, identity, rlen, qlen, origlen,
        );
        let iden_bound = 1f64 - sim_thr;
        let identity = identity
            .min(head_identity + EDGE_BONUS)
            .min(tail_identity + EDGE_BONUS);
        (iden_bound < identity, info)
    };
    if below_dissim {
        trace!("FILLDEL\t{}\tOK", info);
        let (score, new_ops) = infix_guided(unit.seq(), query, &ops, band, ALN_PARAMETER);
        Some((query, aln_start, aln_end, new_ops, score))
    } else {
        trace!("FILLDEL\t{}\tNG", info);
        if log_enabled!(log::Level::Trace) {
            let (xr, ar, yr) = kiley::recover(unitseq, query, &ops);
            for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
                eprintln!("ALN\t{}", String::from_utf8_lossy(xr));
                eprintln!("ALN\t{}", String::from_utf8_lossy(ar));
                eprintln!("ALN\t{}\n", String::from_utf8_lossy(yr));
            }
        }
        None
    }
}

fn edge_identity(unit: &[u8], _: &[u8], ops: &[kiley::Op], len: usize) -> (f64, f64) {
    let (mut head_aln_len, mut head_match) = (0, 0);
    let (mut tail_aln_len, mut tail_match) = (0, 0);
    let head_eval_end = len.min(unit.len());
    let tail_eval_start = unit.len().saturating_sub(len);
    let mut rpos = 0;
    for &op in ops {
        match op {
            kiley::Op::Mismatch | kiley::Op::Match => rpos += 1,
            kiley::Op::Ins => {}
            kiley::Op::Del => rpos += 1,
        }
        if rpos < head_eval_end {
            head_aln_len += 1;
            head_match += (kiley::Op::Match == op) as usize;
        }
        if tail_eval_start <= rpos {
            tail_aln_len += 1;
            tail_match += (kiley::Op::Match == op) as usize;
        }
    }
    let head_identity = head_match as f64 / head_aln_len as f64;
    let tail_identity = tail_match as f64 / tail_aln_len as f64;
    (head_identity, tail_identity)
}

// Triming the head/tail insertion, re-calculate the start and end position.
fn trim_head_tail_insertion(ops: &mut Vec<kiley::Op>) -> (usize, usize) {
    let mut tail_ins = 0;
    while ops.last() == Some(&kiley::Op::Ins) {
        ops.pop();
        tail_ins += 1;
    }
    ops.reverse();
    let mut head_ins = 0;
    while ops.last() == Some(&kiley::Op::Ins) {
        ops.pop();
        head_ins += 1;
    }
    ops.reverse();
    (head_ins, tail_ins)
}

// Align read skeltons to read, return the pileup sumamries.
// i-> insertions before the i-th nodes.
fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    assert!(!read.nodes.is_empty());
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    let units_in_read: HashSet<_> = read.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
    let filtered_read = reads
        .iter()
        .filter(|q| q.nodes.iter().any(|n| units_in_read.contains(&n.key())));
    let skelton = ReadSkelton::new(read);
    for query in filtered_read {
        let (aln, is_forward) = match alignment(&skelton, query) {
            Some(res) => res,
            None => continue,
        };
        let mut q_ptr = SkeltonIter::new(query, is_forward);
        let mut pileups = pileups.iter_mut();
        let mut current_pu = pileups.next().unwrap();
        let mut position = 0;
        for op in aln {
            // These unwraps are safe.
            match op {
                Op::Ins(l) if position == 0 => {
                    // Retain only the last insertion...
                    current_pu.add_tail(q_ptr.nth(l - 1).unwrap());
                }
                Op::Ins(l) if position == pileups.len() - 1 => {
                    // Retain only the first insertion...
                    current_pu.add_head(q_ptr.next().unwrap());
                    for _ in 0..l - 1 {
                        q_ptr.next().unwrap();
                    }
                }
                Op::Ins(l) => {
                    current_pu.add_head(q_ptr.next().unwrap());
                    for _ in 0..l - 1 {
                        current_pu.add_tail(q_ptr.next().unwrap());
                    }
                }
                Op::Del(l) => {
                    current_pu = pileups.nth(l - 1).unwrap();
                    position += l;
                }
                Op::Match(l) => {
                    q_ptr.nth(l - 1);
                    position += l;
                    for _ in 0..l {
                        current_pu.coverage += 1;
                        current_pu = pileups.next().unwrap();
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
// fn alignment(read: &EncodedRead, query: &ReadSkelton) -> Option<(Vec<Op>, bool)> {
fn alignment(read: &ReadSkelton, query: &ReadSkelton) -> Option<(Vec<Op>, bool)> {
    let (f_score, f_ops) = pairwise_alignment_gotoh(read, query);
    let f_match = get_match_units(&f_ops);
    let query = query.rev();
    let (r_score, r_ops) = pairwise_alignment_gotoh(read, &query);
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
    let is_same_dir = x.is_forward == y.is_forward;
    match (x.unit == y.unit, x.cluster == y.cluster, is_same_dir) {
        (_, _, false) | (false, _, _) => MIN_ALN,
        (true, true, true) => 1,
        (true, false, true) => -1,
    }
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
    fn information_head(&self, unit: u64, is_forward: bool) -> (Option<isize>, Option<isize>) {
        let inserts = self
            .head_inserted
            .iter()
            .filter(|node| node.unit == unit && node.is_forward == is_forward);
        Self::summarize(inserts)
    }
    fn information_tail(&self, unit: u64, is_forward: bool) -> (Option<isize>, Option<isize>) {
        let inserts = self
            .tail_inserted
            .iter()
            .filter(|node| node.unit == unit && node.is_forward == is_forward);
        Self::summarize(inserts)
    }
    fn summarize<'a>(
        inserts: impl std::iter::Iterator<Item = &'a LightNode>,
    ) -> (Option<isize>, Option<isize>) {
        let (mut prev_count, mut prev_total) = (0, 0);
        let (mut after_count, mut after_total) = (0, 0);
        for node in inserts {
            if let Some(len) = node.prev_offset {
                prev_count += 1;
                prev_total += len;
            }
            if let Some(len) = node.after_offset {
                after_count += 1;
                after_total += len;
            }
        }
        let prev_offset = (prev_count != 0).then(|| prev_total / prev_count);
        let after_offset = (after_count != 0).then(|| after_total / after_count);
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
        self.insertion_head()
            .filter(|&(num, _, _)| threshold <= num)
            .filter_map(|(_, uid, direction)| {
                let (prev_offset, _) = self.information_head(uid, direction);
                let start_position =
                    (nodes[idx - 1].position_from_start + nodes[idx - 1].query_length()) as isize;
                let start_position = start_position + prev_offset?;
                let start_position = (start_position as usize).saturating_sub(OFFSET);
                Some((start_position, direction, uid))
            })
            .collect()
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
            .filter_map(|(_, uid, direction)| {
                let (_, after_offset) = self.information_tail(uid, direction);
                let end_position = (end_position - after_offset?).max(0) as usize + OFFSET;
                Some((end_position, direction, uid))
            })
            .collect()
    }
}

#[derive(Clone)]
pub struct ReadSkelton {
    nodes: Vec<LightNode>,
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
                    cluster: n.cluster,
                    is_forward: n.is_forward,
                    after_offset,
                }
            })
            .collect();
        // let sets: HashSet<_> = nodes.iter().map(|n| n.unit).collect();
        ReadSkelton { nodes }
    }
    fn rev(&self) -> Self {
        let nodes: Vec<_> = self.nodes.iter().rev().map(LightNode::rev).collect();
        Self { nodes }
    }
}

#[derive(Clone)]
pub struct LightNode {
    // How long should be add to the last position of the previous node to
    // get the start position of this node.
    // None if this is the first node.
    prev_offset: Option<isize>,
    unit: u64,
    cluster: u64,
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
    fn key(&self) -> (u64, u64) {
        (self.unit, self.cluster)
    }
    fn rev(
        &Self {
            prev_offset,
            unit,
            cluster,
            is_forward,
            after_offset,
        }: &Self,
    ) -> Self {
        Self {
            prev_offset: after_offset,
            unit,
            cluster,
            is_forward: !is_forward,
            after_offset: prev_offset,
        }
    }
}

struct SkeltonIter<'a> {
    inner: &'a ReadSkelton,
    index: usize,
    is_forward: bool,
}

impl<'a> SkeltonIter<'a> {
    fn new(read: &'a ReadSkelton, is_forward: bool) -> Self {
        let mut it = Self {
            inner: read,
            is_forward,
            index: 0,
        };
        if !is_forward {
            it.index = read.nodes.len();
        }
        it
    }
}

impl<'a> std::iter::Iterator for SkeltonIter<'a> {
    type Item = LightNode;
    fn next(&mut self) -> Option<Self::Item> {
        if self.is_forward {
            self.index += 1;
            self.inner.nodes.get(self.index - 1).cloned()
        } else if 0 < self.index {
            self.index -= 1;
            self.inner.nodes.get(self.index).map(LightNode::rev)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod deletion_fill {
    use super::*;
    #[test]
    fn aln_test_gotoh() {
        let into_reads = |nodes: Vec<u64>| {
            let nodes: Vec<_> = nodes
                .into_iter()
                .map(|unit| LightNode {
                    prev_offset: None,
                    unit,
                    cluster: 0,
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
