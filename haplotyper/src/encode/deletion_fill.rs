//! Filling deletion.
use crate::ALN_PARAMETER;
use definitions::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
// identity would be increased by this value when evaluating the edges.
const EDGE_BOUND: f64 = 0.5;
// Evaluating length of each side.
const EDGE_LEN: usize = 100;
const INS_THR: usize = 2;
#[derive(Debug, Clone)]
pub struct CorrectDeletionConfig {
    re_clustering: bool,
    sim_thr: Option<f64>,
    stddev_of_error: Option<f64>,
}

impl CorrectDeletionConfig {
    /// If sim_thr is None, it is automatically estimated by `haplotyper::determine_chunks::calc_sim_thr`.
    pub fn new(re_clustering: bool, sim_thr: Option<f64>, stddev_of_error: Option<f64>) -> Self {
        Self {
            re_clustering,
            sim_thr,
            stddev_of_error,
        }
    }
}

pub trait CorrectDeletion {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig);
}
const SEED_CONST: usize = 49823094830;
impl CorrectDeletion for DataSet {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig) {
        let mut find_new_chunks = correct_chunk_deletion(self, config);
        // If half of the coverage supports large deletion, remove them.
        const OCCUPY_FRACTION: f64 = 0.5;
        use crate::purge_diverged::*;
        let p_config = PurgeLargeDelConfig::new(crate::MAX_ALLOWED_GAP, OCCUPY_FRACTION, false);
        find_new_chunks.extend(self.purge_largeindel(&p_config));
        let p_config = PurgeLargeDelConfig::new(crate::MAX_ALLOWED_GAP, OCCUPY_FRACTION, true);
        find_new_chunks.extend(self.purge_largeindel(&p_config));
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
            // Erase cluster
            self.encoded_reads
                .iter_mut()
                .flat_map(|r| r.nodes.iter_mut())
                .for_each(|n| n.cluster = 0);
            use crate::multiplicity_estimation::*;
            let seed = (SEED_CONST * self.encoded_reads.len()) as u64;
            let config = MultiplicityEstimationConfig::new(seed, None);
            self.estimate_multiplicity(&config);
            // Retain all the chunks changed their copy numbers.
            let changed_chunks: HashSet<_> = self
                .selected_chunks
                .iter()
                .filter_map(|c| match prev_copy_numbers.get(&c.id) {
                    None => None,
                    Some(&prev) if c.copy_num == prev => None,
                    Some(_) => Some(c.id),
                })
                .collect();
            // Merge these two.
            let selection: HashSet<_> = find_new_chunks.union(&changed_chunks).copied().collect();
            // Recover the original assignments on the retained chunks.
            self.encoded_reads
                .iter_mut()
                .zip(original_assignments)
                .for_each(|(read, (id, log))| {
                    assert_eq!(read.id, id);
                    recover_original_assignments(read, &log, &selection);
                });
            // Reclustering.
            use crate::local_clustering::LocalClustering;
            self.local_clustering_selected(&selection);
            // By the way, removing zero-copy chunks. Give the upper bound a very large value.
            self.purge_multiplicity(10000000);
        }
    }
}
// Logging the original assignment into a vector.
fn log_original_assignments(ds: &DataSet) -> Vec<(u64, Vec<(u64, u64)>)> {
    ds.encoded_reads
        .iter()
        .map(|r| {
            let xs: Vec<_> = r.nodes.iter().map(|u| (u.chunk, u.cluster)).collect();
            (r.id, xs)
        })
        .collect()
}

// Recover the previous clustering. Note that sometimes the node is added
// so that the length of the read is different from the logged one.
// But it does not removed!
fn recover_original_assignments(read: &mut EncodedRead, log: &[(u64, u64)], except: &HashSet<u64>) {
    let mut read = read.nodes.iter_mut();
    for &(chunk, cluster) in log {
        for node in &mut read {
            if node.chunk == chunk {
                if !except.contains(&node.chunk) {
                    node.cluster = cluster;
                }
                break;
            }
        }
    }
}

/**
The second argument is the vector of (index, chunk_id) of the previous failed trials.
for example, if failed_trials[i][0] = (j,id), then, we've already tried to encode the id-th chunk after the j-th
position of the i-th read, and failed it.
If we can encode some position in the i-th read, the failed trials would be erased, as it change the
condition of the read, making it possible to encode an chunk previously failed to encode.
sim_thr is the similarity threshold.
This function corrects "chunk-deletions". To do that,
it first align other reads onto a read to be corrected in chunk resolution, detecting putative insertions.
Then, it tries to encode these putative insertions in base-pair resolution.
Note that, in the first - chunk resolution - alignment, there's no distinction between clusters.
However, in the second alignment, it tries to encode the putative region by each cluster's representative.
Of course, if there's only one cluster for a chunk, then, it just tries to encode by that chunk.
Auto-tune the similarity threshold.
 */

pub fn correct_chunk_deletion(ds: &mut DataSet, config: &CorrectDeletionConfig) -> HashSet<u64> {
    const OUTER_LOOP: usize = 3;
    let mut find_new_node = HashSet::new();
    let consensi = take_consensus_sequence(ds);
    let fallback = config.sim_thr.unwrap_or_else(|| {
        crate::determine_chunks::calc_sim_thr(ds, crate::determine_chunks::TAKE_THR)
    });
    use crate::estimate_error_rate::estimate_error_rate;
    let errors = estimate_error_rate(ds, fallback);
    let standard_dev = config.stddev_of_error.unwrap_or(errors.median_of_sqrt_err);
    let mut failed_trials: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|r| FailedUpdates::new(r.id))
        .collect();
    for t in 0..OUTER_LOOP {
        ds.encoded_reads.retain(|r| !r.nodes.is_empty());
        debug!("ErrorRateSTDDev\t{}\t{}", t, standard_dev);
        let ft = &mut failed_trials;
        let (new_nodes, is_updated) = filling_until(ds, ft, &consensi, &errors, standard_dev);
        find_new_node.extend(new_nodes);
        {
            let mut idx = 0;
            failed_trials.retain(|_| {
                idx += 1;
                !ds.encoded_reads[idx - 1].nodes.is_empty()
            });
            ds.encoded_reads.retain(|r| !r.nodes.is_empty());
        }
        if !is_updated {
            break;
        }
    }
    find_new_node
}

const INNER_LOOP: usize = 12;
fn filling_until(
    ds: &mut DataSet,
    failed_trials: &mut Vec<FailedUpdates>,
    consensi: &HashMap<(u64, u64), Vec<u8>>,
    error_rates: &crate::estimate_error_rate::ErrorRate,
    stddev: f64,
) -> (HashSet<u64>, bool) {
    ds.encoded_reads.retain(|r| !r.nodes.is_empty());
    failed_trials.iter_mut().for_each(|ft| ft.revive());
    let raw_seq: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let mut find_new_node = HashSet::new();
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    let mut read_skeltons: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
    for i in 0..INNER_LOOP {
        let prev: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
        // let alive = failed_trials.iter().filter(|r| r.is_alive).count();
        //      debug!("Reads\t{}\t{}", alive, failed_trials.len());
        assert_eq!(ds.encoded_reads.len(), failed_trials.len());
        assert_eq!(ds.encoded_reads.len(), read_skeltons.len());
        let new_nodes = ds
            .encoded_reads
            .par_iter_mut()
            .zip(failed_trials.par_iter_mut())
            .filter(|(r, t)| !r.nodes.is_empty() && t.is_alive)
            .flat_map(|(read, fails)| {
                let seq = raw_seq[&read.id];
                let read = (read, seq);
                let chunks = (&chunks, error_rates, consensi);
                correct_deletion_error(read, fails, chunks, stddev, &read_skeltons)
            });
        find_new_node.par_extend(new_nodes);
        updates_reads(&mut read_skeltons, &ds.encoded_reads, failed_trials);
        let after: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
        // debug!("Filled\t{i}\t{prev}\t{after}");
        if after == prev && i == 0 {
            return (find_new_node, false);
        } else if after == prev {
            break;
        }
    }
    (find_new_node, true)
}

fn updates_reads(
    skeltons: &mut [ReadSkelton],
    reads: &[EncodedRead],
    failed_updates: &[FailedUpdates],
) {
    for ((ft, r), sk) in failed_updates.iter().zip(reads.iter()).zip(skeltons.iter()) {
        assert_eq!(ft.readid, r.id);
        assert_eq!(ft.readid, sk.id);
    }
    skeltons
        .iter_mut()
        .zip(reads.iter())
        .zip(failed_updates.iter())
        .filter(|&(_, t)| t.is_alive)
        .for_each(|((skelton, read), _)| *skelton = ReadSkelton::new(read));
}

#[derive(Debug, Clone)]
struct FailedUpdates {
    readid: u64,
    is_alive: bool,
    failed_trials: Vec<(usize, LightNode)>,
}

impl FailedUpdates {
    fn revive(&mut self) {
        self.is_alive = true;
        self.failed_trials.clear();
    }
    fn new(id: u64) -> Self {
        Self {
            readid: id,
            is_alive: true,
            failed_trials: vec![],
        }
    }
    fn extend<I: std::iter::Iterator<Item = (usize, LightNode)>>(&mut self, iter: I) {
        self.failed_trials.extend(iter);
    }
}

// Take consensus of each cluster of each chunk, return the consensus seuqneces.
// UnitID->(clsuterID, its consensus).
// fn take_consensus_sequence(ds: &DataSet) -> HashMap<u64, Vec<(u64, Vec<u8>)>> {
fn take_consensus_sequence(ds: &DataSet) -> HashMap<(u64, u64), Vec<u8>> {
    fn polish(xs: &[&[u8]], chunk: &Chunk, band: usize) -> Vec<u8> {
        kiley::bialignment::guided::polish_until_converge(chunk.seq(), xs, band)
    }
    let ref_chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u)).collect();
    let mut bucket: HashMap<_, Vec<_>> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        bucket
            .entry((node.chunk, node.cluster))
            .or_default()
            .push(node.seq());
    }
    bucket
        .par_iter()
        .filter(|(_, seq)| !seq.is_empty())
        .map(|(key, seqs)| {
            let ref_chunk = &ref_chunks[&key.0];
            let band = ds.read_type.band_width(ref_chunk.seq().len());
            let representative: Vec<_> = match key.1 {
                0 => ref_chunk.seq().to_vec(),
                _ => polish(seqs, ref_chunk, band),
            };
            (*key, representative)
        })
        .collect()
}

#[inline]
fn abs(x: usize, y: usize) -> usize {
    x.max(y) - x.min(y)
}

// Aligment offset. We align [s-offset..e+offset] region to the chunk.
const OFFSET_FACTOR: f64 = 0.1;
// returns the ids of the chunks newly encoded.
// Maybe each (chunk, cluster) should corresponds to a key...?
type UnitInfo<'a> = (
    &'a HashMap<u64, &'a Chunk>,
    &'a crate::estimate_error_rate::ErrorRate,
    &'a HashMap<(u64, u64), Vec<u8>>,
);
fn correct_deletion_error(
    (read, seq): (&mut EncodedRead, &[u8]),
    ft: &mut FailedUpdates,
    chunkinfo: UnitInfo,
    stddev: f64,
    reads: &[ReadSkelton],
) -> Vec<u64> {
    let read_error = chunkinfo.1.read(read.id);
    let pileups = get_pileup(read, reads);
    let nodes = &read.nodes;
    let mut inserts = vec![];
    let ins_thr = mean_cov(&pileups)
        .map(|x| (x / 5).min(INS_THR))
        .unwrap_or(INS_THR);
    for (idx, pileup) in pileups.iter().enumerate() {
        let mut head_cand = pileup.check_insertion_head(nodes, ins_thr, idx);
        head_cand.retain(|node, _| !ft.failed_trials.contains(&(idx, *node)));
        let head_best =
            try_encoding_head(nodes, &head_cand, idx, chunkinfo, seq, read_error, stddev);
        match head_best {
            Some((head_node, _)) => inserts.push((idx, head_node)),
            None => ft.extend(head_cand.keys().map(|&n| (idx, n))),
        }
        let mut tail_cand = pileup.check_insertion_tail(nodes, ins_thr, idx);
        tail_cand.retain(|node, _| !ft.failed_trials.contains(&(idx, *node)));
        let tail_best =
            try_encoding_tail(nodes, &tail_cand, idx, chunkinfo, seq, read_error, stddev);
        // tries.extend(tail_cand.iter().map(|(n, p)| (n.chunk, *p, idx, 1)));
        // oks.extend(tail_best.iter().map(|(n, p)| (n.chunk, *p, idx, 1)));
        match tail_best {
            Some((tail_node, _)) => inserts.push((idx, tail_node)),
            None => ft.extend(tail_cand.into_iter().map(|x| (idx, x.0))),
        }
    }
    let new_inserts: Vec<_> = inserts.iter().map(|(_, n)| n.chunk).collect();
    ft.is_alive = !inserts.is_empty();
    if !inserts.is_empty() {
        ft.revive();
        read.nodes.extend(inserts.into_iter().map(|x| x.1));
    }
    if ft.is_alive && !read.nodes.is_empty() {
        let mut nodes = Vec::with_capacity(read.nodes.len());
        nodes.append(&mut read.nodes);
        use super::{nodes_to_encoded_read, remove_overlapping_encoding, remove_slippy_alignment};
        nodes.sort_by_key(|n| (n.chunk, n.position_from_start));
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        nodes = remove_overlapping_encoding(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
    new_inserts
}

const THR: f64 = 10f64;
fn try_encoding_head(
    nodes: &[Node],
    head_cand: &HashMap<LightNode, usize>,
    idx: usize,
    (chunks, chunk_error_rate, consensi): UnitInfo,
    seq: &[u8],
    read_error: f64,
    stddev: f64,
) -> Option<(Node, i32)> {
    head_cand
        .iter()
        .filter(|&(_, &start_position)| start_position < seq.len())
        .filter_map(|(node, &start_position)| {
            let (uid, cluster) = (node.chunk, node.cluster);
            let chunk = *chunks.get(&uid)?;
            let cons = consensi.get(&(uid, cluster))?;
            let offset = (OFFSET_FACTOR * cons.len() as f64).ceil() as usize;
            let end_position = (start_position + cons.len() + offset).min(seq.len());
            let start_position = start_position.saturating_sub(offset);
            let is_the_same_encode = match nodes.get(idx) {
                Some(node) => {
                    node.chunk == uid && abs(node.position_from_start, start_position) < cons.len()
                }
                None => false,
            };
            assert!(start_position < end_position);
            if is_the_same_encode {
                return None;
            }
            let expected = read_error + chunk_error_rate.chunk((uid, cluster));
            let error_rate_bound = expected + THR * stddev;
            let position = (start_position, end_position, node.is_forward);
            let chunk_info = (chunk, cluster, cons.as_slice());
            encode_node(seq, position, chunk_info, error_rate_bound)
        })
        .max_by_key(|x| x.1)
}

fn try_encoding_tail(
    nodes: &[Node],
    tail_cand: &HashMap<LightNode, usize>,
    idx: usize,
    (chunks, chunk_error_rate, consensi): UnitInfo,
    seq: &[u8],
    read_error: f64,
    stddev: f64,
) -> Option<(Node, i32)> {
    tail_cand
        .iter()
        .filter_map(|(node, &end_position)| {
            let (uid, cluster) = (node.chunk, node.cluster);
            let chunk = *chunks.get(&uid)?;
            let cons = consensi.get(&(uid, cluster))?;
            let offset = (OFFSET_FACTOR * cons.len() as f64).ceil() as usize;
            let start_position = end_position
                .min(seq.len())
                .saturating_sub(offset + cons.len());
            let end_position = (end_position + offset).min(seq.len());
            assert!(start_position < end_position);
            let is_the_same_encode = match nodes.get(idx) {
                Some(node) => {
                    node.chunk == uid && abs(node.position_from_start, start_position) < cons.len()
                }
                None => false,
            };
            if is_the_same_encode {
                return None;
            }
            assert!(start_position < end_position);
            let positions = (start_position, end_position, node.is_forward);
            let chunk_info = (chunk, cluster, cons.as_slice());
            let expected = read_error + chunk_error_rate.chunk((uid, cluster));
            let error_rate_bound = expected + THR * stddev;
            encode_node(seq, positions, chunk_info, error_rate_bound)
        })
        .max_by_key(|x| x.1)
}

// Try to Encode Node. Return Some(node) if the alignment is good.
// Return also the alignment score of the encoding.
// The match score is 2, mism is -6, gap open is -5, and gap ext is -1.
fn encode_node(
    query: &[u8],
    (start, end, is_forward): (usize, usize, bool),
    (chunk, cluster, chunkseq): (&Chunk, u64, &[u8]),
    sim_thr: f64,
) -> Option<(Node, i32)> {
    // Initial filter.
    // If the query is shorter than the chunkseq,
    // at least we need the edit operations to fill the gaps.
    // This is lower bound of the sequence identity.
    let edit_dist_lower_bound = chunkseq.len().saturating_sub(end - start);
    let diff_lower_bound = edit_dist_lower_bound as f64 / chunkseq.len() as f64;
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
    let (seq, trim_head, trim_tail, kops, score) =
        fine_mapping(&query, (chunk, cluster, chunkseq), sim_thr)?;
    let ops = crate::misc::kiley_op_to_ops(&kops).0;
    let cl = chunk.cluster_num;
    let position_from_start = match is_forward {
        true => start + trim_head,
        false => start + trim_tail,
    };
    // I think we should NOT make likelihood gain to some biased value,
    // as 1. if the alignment gives the certaintly, then we can impute the clustering by the alignment,
    // 2. if `cluster` assignment is just by chance,
    // then we just should not introduce any bias into the likelihood gain.
    let seq = seq.to_vec();
    let mut node = Node::new(chunk.id, is_forward, seq, ops, position_from_start, cl);
    node.cluster = cluster;
    Some((node, score))
}

// Sequence, trimed base from the head, trimed base from the tail, ops, score.
type FineMapping<'a> = (&'a [u8], usize, usize, Vec<kiley::Op>, i32);
// const EDLIB_OFS: f64 = 0.10;
fn fine_mapping<'a>(
    orig_query: &'a [u8],
    (chunk, cluster, chunkseq): (&Chunk, u64, &[u8]),
    sim_thr: f64,
) -> Option<FineMapping<'a>> {
    let (query, trim_head, trim_tail, ops, band) =
        fit_query_by_edlib(chunkseq, orig_query, sim_thr)?;
    let mat_num = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
    let identity = mat_num as f64 / ops.len() as f64;
    let (head_identity, tail_identity) = edge_identity(chunkseq, query, &ops, EDGE_LEN);
    let iden_bound = 1f64 - sim_thr;
    let below_dissim = iden_bound < identity && EDGE_BOUND < head_identity.min(tail_identity);
    {
        let (rlen, qlen) = (chunkseq.len(), query.len());
        let id = chunk.id;
        let orig_len = orig_query.len();
        let info = format!(
            "{id}\t{cluster}\t{identity:.2}\t{rlen}\t{qlen}\t{orig_len}\t{trim_head}\t{trim_tail}"
        );
        trace!("FILLDEL\t{}\t{}", info, ["NG", "OK"][below_dissim as usize]);
    }
    if log_enabled!(log::Level::Trace) && !below_dissim {
        let (xr, ar, yr) = kiley::op::recover(chunkseq, query, &ops);
        for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
            eprintln!("ALN\t{}", String::from_utf8_lossy(xr));
            eprintln!("ALN\t{}", String::from_utf8_lossy(ar));
            eprintln!("ALN\t{}", String::from_utf8_lossy(yr));
        }
    }
    if !below_dissim {
        return None;
    }
    let (score, new_ops) = infix_guided(chunk.seq(), query, &ops, band, ALN_PARAMETER);
    Some((query, trim_head, trim_tail, new_ops, score))
}

use kiley::bialignment::guided::infix_guided;
fn edlib_op_to_kiley_op(ops: &[u8]) -> Vec<kiley::Op> {
    use kiley::Op::*;
    ops.iter()
        .map(|&op| [Match, Ins, Del, Mismatch][op as usize])
        .collect()
}

type FitQuery<'a> = (&'a [u8], usize, usize, Vec<kiley::op::Op>, usize);
fn fit_query_by_edlib<'a>(
    chunkseq: &[u8],
    orig_query: &'a [u8],
    sim_thr: f64,
) -> Option<FitQuery<'a>> {
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    let alignment = edlib_sys::align(chunkseq, orig_query, mode, task);
    let (start, end) = alignment.location().unwrap();
    let mut ops = vec![kiley::Op::Del; start];
    ops.extend(edlib_op_to_kiley_op(alignment.operations().unwrap()));
    ops.extend(std::iter::repeat(kiley::Op::Del).take(orig_query.len() - end - 1));
    // Align twice, to get an accurate alignment.
    let band = ((orig_query.len() as f64 * sim_thr * 0.3).ceil() as usize).max(10);
    let (_, ops) = infix_guided(orig_query, chunkseq, &ops, band, ALN_PARAMETER);
    let (_, mut ops) = infix_guided(orig_query, chunkseq, &ops, band, ALN_PARAMETER);
    // Reverse ops
    for op in ops.iter_mut() {
        *op = match *op {
            kiley::Op::Ins => kiley::Op::Del,
            kiley::Op::Del => kiley::Op::Ins,
            x => x,
        }
    }
    let (trim_head, trim_tail) = trim_head_tail_insertion(&mut ops);
    let query = &orig_query[trim_head..orig_query.len() - trim_tail];
    Some((query, trim_head, trim_tail, ops, band))
}

fn edge_identity(chunk: &[u8], _: &[u8], ops: &[kiley::Op], len: usize) -> (f64, f64) {
    let (mut head_aln_len, mut head_match) = (0, 0);
    let (mut tail_aln_len, mut tail_match) = (0, 0);
    let head_eval_end = len.min(chunk.len());
    let tail_eval_start = chunk.len().saturating_sub(len);
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

fn check_alignment_by_chunkmatch(chunks: &[(u64, u64, bool)], query: &ReadSkelton) -> Option<bool> {
    fn count_match(chunks: &[(u64, u64, bool)], query: &[(u64, u64, bool)]) -> usize {
        let mut r_ptr = chunks.iter().peekable();
        let mut q_ptr = query.iter().peekable();
        let mut match_num = 0;
        while r_ptr.peek().is_some() && q_ptr.peek().is_some() {
            match r_ptr.peek().unwrap().cmp(q_ptr.peek().unwrap()) {
                std::cmp::Ordering::Less => r_ptr.next(),
                std::cmp::Ordering::Equal => {
                    match_num += 1;
                    r_ptr.next();
                    q_ptr.next()
                }
                std::cmp::Ordering::Greater => q_ptr.next(),
            };
        }
        match_num
    }
    let mut keys: Vec<_> = query.nodes.iter().map(|n| n.key()).collect();
    keys.sort_unstable();
    let forward_match = count_match(chunks, &keys);
    keys.iter_mut().for_each(|x| x.2 = !x.2);
    keys.sort_unstable();
    let reverse_match = count_match(chunks, &keys);
    let min_match = MIN_MATCH.min(chunks.len());
    (min_match <= forward_match.max(reverse_match)).then_some(reverse_match <= forward_match)
}

// Align read skeltons to read, return the pileup sumamries.
// i-> insertions before the i-th nodes.
// The coverage of the last slot is always zero.
fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    assert!(!read.nodes.is_empty());
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    let skelton = ReadSkelton::new(read);
    let mut chunks_in_read: Vec<_> = skelton.nodes.iter().map(|n| n.key()).collect();
    chunks_in_read.sort_unstable();
    for query in reads.iter() {
        let is_forward = match check_alignment_by_chunkmatch(&chunks_in_read, query) {
            Some(is_forward) => is_forward,
            None => continue,
        };
        let id = read.id;
        let aln = match alignment(id, &skelton, query, is_forward) {
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
                    if 2 <= l {
                        current_pu.add_tail(q_ptr.nth(l - 2).unwrap());
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
// Minimum required chunks to be matched.
const MIN_MATCH: usize = 2;
// Minimum required alignment score.
const SCORE_THR: i32 = 1;
fn alignment(_: u64, read: &ReadSkelton, query: &ReadSkelton, dir: bool) -> Option<Vec<Op>> {
    // let mut query = query.clone();
    let (score, ops) = match dir {
        true => pairwise_alignment_gotoh(read, query),
        false => {
            let query = query.rev();
            pairwise_alignment_gotoh(read, &query)
        }
    };
    let match_num = get_match_chunks(&ops);
    let min_match = MIN_MATCH.min(read.nodes.len()).min(query.nodes.len());
    (min_match <= match_num && SCORE_THR <= score && is_proper(&ops)).then_some(ops)
}

// Return true if the alignment is proper dovetail.
fn is_proper(ops: &[Op]) -> bool {
    ops.windows(2)
        .all(|xs| !matches!(xs, &[Op::Ins(_), Op::Del(_)] | &[Op::Del(_), Op::Ins(_)]))
}

const MIN_ALN: i32 = -10000000;
fn score(x: &LightNode, y: &LightNode) -> i32 {
    if x.chunk != y.chunk || x.is_forward != y.is_forward {
        MIN_ALN
    } else if x.cluster == y.cluster {
        1
    } else {
        -1
    }
}

fn pairwise_alignment_gotoh(read: &ReadSkelton, query: &ReadSkelton) -> (i32, Vec<Op>) {
    let (read, query) = (&read.nodes, &query.nodes);
    let (row_num, col_num) = (read.len() + 1, query.len() + 1);
    let mut dp = vec![0; row_num * col_num * 3];
    let read_row = col_num * 3;
    // Initialize.
    for i in 0..read.len() + 1 {
        dp[read_row * i] = MIN_ALN;
        dp[read_row * i + 1] = MIN_ALN;
    }
    for j in 0..query.len() + 1 {
        dp[3 * j] = MIN_ALN;
        dp[3 * j + 2] = MIN_ALN;
    }
    dp[0] = 0;
    // Filling DP Table.
    for (i, x) in read.iter().enumerate() {
        for (j, y) in query.iter().enumerate() {
            let (i, j) = (i + 1, j + 1);
            let fill_pos = read_row * i + 3 * j;
            let prev_match = read_row * (i - 1) + 3 * (j - 1);
            dp[fill_pos] = dp[prev_match..prev_match + 3].iter().max().unwrap() + score(x, y);
            let prev_ins = read_row * i + 3 * (j - 1);
            dp[fill_pos + 1] = (dp[prev_ins] - 1).max(dp[prev_ins + 1]);
            let prev_del = read_row * (i - 1) + 3 * j;
            dp[fill_pos + 2] = (dp[prev_del] - 1).max(dp[prev_del + 2]);
        }
    }
    let (mut r_pos, mut q_pos, mut state, dist) = (0..read.len() + 1)
        .map(|i| (i, query.len()))
        .chain((0..query.len() + 1).map(|j| (read.len(), j)))
        .filter_map(|(i, j)| {
            let position = read_row * i + 3 * j;
            dp[position..position + 3]
                .iter()
                .enumerate()
                .max_by_key(|x| x.1)
                .map(|(state, &score)| (i, j, state, score))
        })
        .max_by_key(|x| x.3)
        .unwrap();
    let mut ops = Vec::with_capacity(row_num.max(col_num) + 2);
    if read.len() != r_pos {
        ops.push(Op::Del(read.len() - r_pos));
    }
    if query.len() != q_pos {
        ops.push(Op::Ins(query.len() - q_pos));
    }
    while 0 < r_pos && 0 < q_pos {
        let current_pos = read_row * r_pos + 3 * q_pos + state;
        let current_dist = dp[current_pos];
        if state == 0 {
            let dist = current_dist - score(&read[r_pos - 1], &query[q_pos - 1]);
            let prev_pos = read_row * (r_pos - 1) + 3 * (q_pos - 1);
            let (new_state, _) = dp[prev_pos..prev_pos + 3]
                .iter()
                .enumerate()
                .find(|&(_, &score)| score == dist)
                .unwrap();
            state = new_state;
            ops.push(Op::Match(1));
            r_pos -= 1;
            q_pos -= 1;
        } else if state == 1 {
            let prev_pos = read_row * r_pos + 3 * (q_pos - 1);
            state = (current_dist != dp[prev_pos] - 1) as usize;
            ops.push(Op::Ins(1));
            q_pos -= 1;
        } else {
            let prev_pos = read_row * (r_pos - 1) + 3 * q_pos;
            state = if current_dist == dp[prev_pos] - 1 {
                0
            } else {
                2
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
    let mut compressed = Vec::with_capacity(ops.len());
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

fn get_match_chunks(ops: &[Op]) -> usize {
    ops.iter()
        .map(|op| match op {
            Op::Match(l) => *l,
            _ => 0,
        })
        .sum::<usize>()
}

#[derive(Debug, Clone)]
struct Pileup {
    // insertion at the beggining of this node
    head_inserted: Vec<LightNode>,
    // insertion at the last of this node
    tail_inserted: Vec<LightNode>,
    coverage: usize,
}

fn mean_cov(pileups: &[Pileup]) -> Option<usize> {
    let len = pileups.len();
    let sum: usize = pileups.iter().map(|p| p.coverage).sum();
    match len {
        0 => None,
        _ => Some(sum / len),
    }
}

impl Pileup {
    // Return the maximum insertion from the same chunk, the same direction.
    fn insertion_head(&self) -> HashMap<LightNode, usize> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.head_inserted.iter() {
            *count.entry(*node).or_default() += 1;
        }
        count
    }
    fn insertion_tail(&self) -> HashMap<LightNode, usize> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.tail_inserted.iter() {
            *count.entry(*node).or_default() += 1;
        }
        count
    }
    fn new() -> Self {
        Self {
            head_inserted: Vec::with_capacity(5),
            tail_inserted: Vec::with_capacity(5),
            coverage: 0,
        }
    }
    fn information_head(&self, node: &LightNode) -> Option<isize> {
        Self::summarize(&self.head_inserted, node).0
    }
    fn information_tail(&self, node: &LightNode) -> Option<isize> {
        Self::summarize(&self.tail_inserted, node).1
    }
    fn summarize(inserts: &[LightNode], target: &LightNode) -> (Option<isize>, Option<isize>) {
        let inserts = inserts.iter().filter(|&node| node == target);
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
        let prev_offset = match prev_count {
            0 => None,
            _ => Some(prev_total / prev_count),
        };
        let after_offset = match after_count {
            0 => None,
            _ => Some(after_total / after_count),
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
    ) -> HashMap<LightNode, usize> {
        let mut inserts = self.insertion_head();
        inserts.retain(|_, num| threshold <= *num);
        inserts.retain(|node, num| {
            let prev_offset = self.information_head(node);
            let start_position = nodes[idx - 1].position_from_start + nodes[idx - 1].query_length();
            match prev_offset {
                Some(x) => {
                    *num = (start_position as isize + x) as usize;
                    true
                }
                None => false,
            }
        });
        inserts
    }
    fn check_insertion_tail(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> HashMap<LightNode, usize> {
        let end_position = match nodes.get(idx) {
            Some(res) => res.position_from_start as isize,
            None => return HashMap::new(),
        };
        let mut inserts = self.insertion_tail();
        inserts.retain(|_, num| threshold <= *num);
        inserts.retain(|node, num| match self.information_tail(node) {
            Some(after_offset) => {
                *num = (end_position - after_offset).max(0) as usize;
                true
            }
            None => false,
        });
        inserts
    }
}

#[derive(Clone)]
struct ReadSkelton {
    id: u64,
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
    fn new(read: &EncodedRead) -> Self {
        Self::from_rich_nodes(read.id, &read.nodes)
    }
    fn from_rich_nodes(id: u64, nodes: &[Node]) -> Self {
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
                let prev_offset = prev_end.map(|x| summaries[i].0 - x.1);
                let after_offset = summaries.get(i + 1).map(|x| x.0 - summaries[i].1);
                LightNode {
                    prev_offset,
                    chunk: n.chunk,
                    cluster: n.cluster,
                    is_forward: n.is_forward,
                    after_offset,
                }
            })
            .collect();
        ReadSkelton { id, nodes }
    }
    fn rev(&self) -> Self {
        let nodes: Vec<_> = self.nodes.iter().rev().map(LightNode::rev).collect();
        Self { id: self.id, nodes }
    }
}

#[derive(Clone, Copy)]
struct LightNode {
    // How long should be add to the last position of the previous node to
    // get the start position of this node.
    // None if this is the first node.
    prev_offset: Option<isize>,
    chunk: u64,
    cluster: u64,
    is_forward: bool,
    // Almost the same as prev_offset. The distance between the last postion of this node to
    // the start position of the next node.
    // None if this is the last node.
    after_offset: Option<isize>,
}

impl std::cmp::PartialEq for LightNode {
    fn eq(&self, other: &Self) -> bool {
        self.chunk == other.chunk
            && self.cluster == other.cluster
            && self.is_forward == other.is_forward
    }
}

impl std::cmp::Eq for LightNode {}

impl std::hash::Hash for LightNode {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.chunk.hash(state);
        self.cluster.hash(state);
        self.is_forward.hash(state);
    }
}

impl std::fmt::Debug for LightNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (u, c) = (self.chunk, self.cluster);
        let dir = if self.is_forward { '+' } else { '-' };
        let prev = self.prev_offset.unwrap_or(-1);
        let after = self.after_offset.unwrap_or(-1);
        write!(f, "{u}-{c}({dir},{prev},{after})")
    }
}

impl LightNode {
    fn key(&self) -> (u64, u64, bool) {
        (self.chunk, self.cluster, self.is_forward)
    }
    fn rev(
        &Self {
            prev_offset,
            chunk,
            cluster,
            is_forward,
            after_offset,
        }: &Self,
    ) -> Self {
        Self {
            prev_offset: after_offset,
            chunk,
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
mod deletion_fill_test {
    use super::*;
    #[test]
    fn aln_test_gotoh() {
        let into_reads = |nodes: Vec<u64>| {
            let nodes: Vec<_> = nodes
                .into_iter()
                .map(|chunk| LightNode {
                    prev_offset: None,
                    chunk,
                    cluster: 0,
                    is_forward: true,
                    after_offset: None,
                })
                .collect();
            ReadSkelton { id: 0, nodes }
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
