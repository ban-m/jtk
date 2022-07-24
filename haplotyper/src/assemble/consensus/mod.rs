use crate::assemble::ditch_graph::ContigEncoding;
use crate::assemble::ditch_graph::UnitAlignmentInfo;
use crate::model_tune::get_model;
use definitions::*;
use gfa::Segment;
use kiley::hmm::guided::PairHiddenMarkovModel;
use kiley::Op;
use rand::prelude::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;

pub trait Polish: private::Sealed {
    fn polish_segment(
        &self,
        segments: &[Segment],
        encs: &[ContigEncoding],
        config: &PolishConfig,
    ) -> Vec<Segment>;
    fn distribute_to_contig(
        &self,
        segments: &[Segment],
        encs: &[ContigEncoding],
        config: &PolishConfig,
    ) -> BTreeMap<String, Vec<Alignment>>;
}

use std::collections::BTreeMap;
mod private {
    pub trait Sealed {}
    impl Sealed for definitions::DataSet {}
}

#[derive(Debug, Clone)]
pub struct PolishConfig {
    seed: u64,
    min_coverage: usize,
    window_size: usize,
    radius: usize,
    round_num: usize,
}
impl PolishConfig {
    pub fn new(
        seed: u64,
        min_coverage: usize,
        window_size: usize,
        radius: usize,
        round_num: usize,
    ) -> PolishConfig {
        Self {
            seed,
            min_coverage,
            window_size,
            radius,
            round_num,
        }
    }
}

impl std::default::Default for PolishConfig {
    fn default() -> Self {
        Self {
            seed: Default::default(),
            min_coverage: 3,
            window_size: 2000,
            radius: 100,
            round_num: 2,
        }
    }
}

use rayon::prelude::*;
impl Polish for DataSet {
    fn polish_segment(
        &self,
        segments: &[Segment],
        encs: &[ContigEncoding],
        config: &PolishConfig,
    ) -> Vec<Segment> {
        let mut alignments_on_contigs = self.distribute_to_contig(segments, encs, config);
        let hmm = get_model(self).unwrap();
        let polished = alignments_on_contigs.iter_mut().map(|(sid, alignments)| {
            log::debug!("POLISH\t{sid}");
            let seg = segments.iter().find(|seg| &seg.sid == sid).unwrap();
            let seg = seg.sequence.as_ref().unwrap().as_bytes();
            let seq = polish(sid, seg, alignments, &hmm, config);
            (sid.clone(), seq)
        });
        polished
            .into_iter()
            .map(|(sid, seq)| {
                let slen = seq.len();
                let sequence = String::from_utf8(seq).ok();
                Segment::from(sid, slen, sequence)
            })
            .collect()
    }

    fn distribute_to_contig(
        &self,
        segments: &[Segment],
        encs: &[ContigEncoding],
        config: &PolishConfig,
    ) -> BTreeMap<String, Vec<Alignment>> {
        let alignments: Vec<_> = self
            .encoded_reads
            .par_iter()
            .enumerate()
            .map(|(i, read)| {
                let seed = config.seed + i as u64;
                let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
                let mut raw_read = read.recover_raw_read();
                raw_read.iter_mut().for_each(u8::make_ascii_uppercase);
                align_to_contigs(read, &raw_read, encs, segments, &mut rng)
            })
            .collect();
        let mut alignments_on_contigs: BTreeMap<_, Vec<_>> = BTreeMap::new();
        for aln in alignments.into_iter().flatten() {
            alignments_on_contigs
                .entry(aln.contig.clone())
                .or_default()
                .push(aln);
        }
        alignments_on_contigs
    }
}

fn log_identity(sid: &str, alignments: &[Alignment]) {
    if log_enabled!(log::Level::Debug) {
        let (mut mat, mut len) = (0, 0);
        for aln in alignments.iter() {
            let match_num = aln.ops.iter().filter(|&&op| op == Op::Match).count();
            let aln_len = aln.ops.len();
            mat += match_num;
            len += aln_len;
            assert!(len < 2 * mat);
        }
        log::debug!("ALN\t{sid}\t{mat}\t{len}\t{}", mat as f64 / len as f64);
    }
}

fn polish(
    sid: &str,
    draft: &[u8],
    alignments: &mut [Alignment],
    hmm: &PairHiddenMarkovModel,
    config: &PolishConfig,
) -> Vec<u8> {
    let window = config.window_size;
    let mut polished = draft.to_vec();
    log_identity(sid, alignments);
    for _ in 0..config.round_num {
        // Allocation.
        let num_slot = polished.len() / window + (polished.len() % window != 0) as usize;
        let mut pileup_seq = vec![vec![]; num_slot];
        let mut pileup_ops = vec![vec![]; num_slot];
        let allocated_positions: Vec<_> = alignments
            .iter()
            .map(|aln| {
                let mut put_position = vec![];
                let (start, chunks, end) = split(aln, window, num_slot, polished.len());
                for (pos, seq, ops) in chunks {
                    assert_eq!(pileup_seq[pos].len(), pileup_ops[pos].len());
                    put_position.push((pos, pileup_seq[pos].len()));
                    pileup_seq[pos].push(seq);
                    pileup_ops[pos].push(ops);
                }
                (start, put_position, end)
            })
            .collect();
        // for (i, pu) in pileup_seq.iter().enumerate() {
        //     debug!("PILEUP\t{i}\t{}", pu.len());
        // }
        //  Polish
        let polished_seg: Vec<_> = polished
            .par_chunks(window)
            .zip(pileup_seq.par_iter())
            .zip(pileup_ops.par_iter_mut())
            .map(
                |((draft, seqs), ops)| match seqs.len() < config.min_coverage {
                    true => draft.to_vec(),
                    false => polish_seg(hmm, draft, seqs, ops, config.radius),
                },
            )
            .collect();
        let (acc_len, _) = polished_seg.iter().fold((vec![0], 0), |(mut xs, len), x| {
            let len = len + x.len();
            xs.push(len);
            (xs, len)
        });
        assert_eq!(acc_len.len(), polished_seg.len() + 1);
        polished = polished_seg.iter().flatten().copied().collect();
        // Fix alignment.
        alignments
            .iter_mut()
            .zip(allocated_positions.iter())
            .for_each(|(aln, aloc_pos)| {
                fix_alignment(aln, aloc_pos, &polished, &acc_len, &pileup_ops)
            });
        log_identity(sid, alignments);
    }
    polished
}

fn polish_seg(
    hmm: &PairHiddenMarkovModel,
    draft: &[u8],
    seqs: &[&[u8]],
    ops: &mut [Vec<Op>],
    radius: usize,
) -> Vec<u8> {
    use kiley::bialignment::guided::polish_until_converge_with;
    let draft = polish_until_converge_with(draft, seqs, ops, radius);
    // Check validity.
    let has_large_del = seqs.iter().zip(ops.iter()).any(|(seq, op)| {
        let lk_pre = hmm.likelihood_guided(&draft, seq, op, radius);
        let lk_pos = hmm.likelihood_guided_post(&draft, seq, op, radius);
        lk_pre.is_nan() || lk_pos.is_nan()
    });
    match has_large_del {
        true => draft,
        false => hmm.polish_until_converge_with(&draft, seqs, ops, radius),
    }
}

fn fix_alignment(
    aln: &mut Alignment,
    aloc_pos: &(TipPos, Vec<(usize, usize)>, TipPos),
    polished: &[u8],
    acc_len: &[usize],
    pileup_ops: &[Vec<Vec<Op>>],
) {
    aln.ops.clear();
    let &((start, first_pos), ref allocated_pos, (end, last_pos)) = aloc_pos;
    if start == end && first_pos == last_pos {
        // Contained alignment.
        let (contig_start, contig_end) = (acc_len[first_pos], acc_len[first_pos + 1]);
        let (start, ops, end) = align_infix(&aln.query, &polished[contig_start..contig_end]);
        aln.ops.extend(ops);
        aln.contig_start = contig_start + start;
        aln.contig_end = contig_start + end;
        return;
    }
    if start != 0 {
        let first_pos_bp = acc_len[first_pos + 1];
        let (ops, contig_len) = align_leading(&aln.query[..start], &polished[..first_pos_bp]);
        aln.contig_start = first_pos_bp - contig_len;
        aln.ops.extend(ops);
    } else {
        aln.contig_start = acc_len[first_pos];
    }
    for &(pos, idx) in allocated_pos.iter() {
        aln.ops.extend(pileup_ops[pos][idx].iter());
    }
    if end != aln.query.len() {
        let last_pos_bp = acc_len[last_pos];
        let (ops, contig_len) = align_trailing(&aln.query[end..], &polished[last_pos_bp..]);
        aln.ops.extend(ops);
        aln.contig_end = last_pos_bp + contig_len;
    } else {
        aln.contig_end = acc_len[last_pos];
    }
}

fn align_infix(query: &[u8], seg: &[u8]) -> (usize, Vec<Op>, usize) {
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    let aln = edlib_sys::align(query, seg, mode, task);
    let (start, end) = aln.location().unwrap();
    let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
    (start, ops, end + 1)
}

// ------------> seg
//       ------> query
fn align_leading(query: &[u8], seg: &[u8]) -> (Vec<Op>, usize) {
    let len = query.len();
    let seg = &seg[seg.len().saturating_sub(2 * len)..];
    if query.is_empty() {
        (Vec::new(), 0)
    } else if seg.is_empty() {
        (vec![Op::Ins; query.len()], 0)
    } else {
        let mode = edlib_sys::AlignMode::Infix;
        let task = edlib_sys::AlignTask::Alignment;
        let aln = edlib_sys::align(query, seg, mode, task);
        let (_, end) = aln.location().unwrap();
        let end = end + 1;
        let mut ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
        let rem_len = seg.len() - end;
        ops.extend(std::iter::repeat(Op::Del).take(rem_len));
        let seg_len = ops.iter().filter(|&&op| op != Op::Ins).count();
        (ops, seg_len)
    }
}

// -----------------> seg
// ----->             query
fn align_trailing(query: &[u8], seg: &[u8]) -> (Vec<Op>, usize) {
    let len = query.len();
    let seg = &seg[..(2 * len).min(seg.len())];
    if query.is_empty() {
        (Vec::new(), 0)
    } else if seg.is_empty() {
        (vec![Op::Ins; query.len()], 0)
    } else {
        let mode = edlib_sys::AlignMode::Prefix;
        let task = edlib_sys::AlignTask::Alignment;
        let aln = edlib_sys::align(query, seg, mode, task);
        let (_, end) = aln.location().unwrap();
        let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
        (ops, end + 1)
    }
}

const EDGE: usize = 100;
// (bp position in the query, chunk id in the contig)
type TipPos = (usize, usize);
type Chunk<'a> = (usize, &'a [u8], Vec<Op>);
fn split(
    alignment: &Alignment,
    window: usize,
    window_num: usize,
    contig_len: usize,
) -> (TipPos, Vec<Chunk>, TipPos) {
    let len = alignment.ops.iter().filter(|&&op| op != Op::Del).count();
    assert_eq!(len, alignment.query.len());
    let refr = alignment.ops.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(refr, alignment.contig_end - alignment.contig_start);
    let start_chunk_id = alignment.contig_start / window;
    let start_pos_in_contig = match alignment.contig_start % window == 0 {
        true => start_chunk_id * window,
        false => (start_chunk_id + 1) * window,
    };
    let (mut qpos, mut cpos) = (0, alignment.contig_start);
    let mut ops = alignment.ops.iter();
    // Seek to the start position.
    while cpos < start_pos_in_contig {
        match ops.next() {
            Some(Op::Match) | Some(Op::Mismatch) => {
                qpos += 1;
                cpos += 1;
            }
            Some(Op::Ins) => qpos += 1,
            Some(Op::Del) => cpos += 1,
            None => break,
        }
    }
    let start_pos_in_query = qpos;
    if cpos < start_pos_in_contig {
        let start = (start_pos_in_query, start_chunk_id);
        return (start, Vec::new(), start);
    }
    let mut chunks = vec![];
    let mut current_chunk_id = start_pos_in_contig / window;
    let mut end_pos = qpos;
    'outer: loop {
        let start = qpos;
        let mut chunk_ops = vec![];
        let target = (current_chunk_id + 1) * window;
        while cpos < target {
            match ops.next() {
                Some(Op::Match) => {
                    cpos += 1;
                    qpos += 1;
                    chunk_ops.push(Op::Match);
                }
                Some(Op::Mismatch) => {
                    cpos += 1;
                    qpos += 1;
                    chunk_ops.push(Op::Mismatch);
                }
                Some(Op::Del) => {
                    cpos += 1;
                    chunk_ops.push(Op::Del);
                }
                Some(Op::Ins) => {
                    qpos += 1;
                    chunk_ops.push(Op::Ins);
                }
                None if current_chunk_id == window_num - 1 && (contig_len - cpos) < EDGE => {
                    chunks.push((current_chunk_id, &alignment.query[start..qpos], chunk_ops));
                    end_pos = qpos;
                    current_chunk_id += 1;
                    break 'outer;
                }
                None => break 'outer,
            }
        }
        chunks.push((current_chunk_id, &alignment.query[start..qpos], chunk_ops));
        end_pos = qpos;
        current_chunk_id += 1;
    }
    let start = (start_pos_in_query, start_chunk_id);
    let end = (end_pos, current_chunk_id);
    (start, chunks, end)
}

fn align_to_contigs<R: Rng>(
    read: &EncodedRead,
    seq: &[u8],
    encs: &[ContigEncoding],
    segs: &[gfa::Segment],
    rng: &mut R,
) -> Vec<Alignment> {
    let mut chains = enumerate_chain(read, encs);
    let mut alns = vec![];
    while !chains.is_empty() {
        let choises: Vec<_> = (0..chains.len()).collect();
        let max = chains
            .iter()
            .map(|chain| chain.apporox_score())
            .max()
            .unwrap();
        let picked = choises
            .choose_weighted(rng, |&idx| {
                ((chains[idx].apporox_score() - max) as f64).exp()
            })
            .unwrap();
        let chain = chains.remove(*picked);
        chains.retain(|c| c.overlap_frac(&chain) < 0.5);
        let seg = segs.iter().find(|seg| seg.sid == chain.id).unwrap();
        let enc = encs.iter().find(|enc| enc.id == chain.id).unwrap();
        alns.push(base_pair_alignment(read, seq, &chain, seg, enc));
    }
    alns
}

fn enumerate_chain(read: &EncodedRead, encs: &[ContigEncoding]) -> Vec<Chain> {
    let mut chains = vec![];
    let mut nodes_run: Vec<_> = read
        .nodes
        .iter()
        .map(|n| {
            let start = n.position_from_start;
            let end = start + n.query_length();
            LightNode::new((n.unit, n.cluster), n.is_forward, (start, end))
        })
        .collect();
    for enc in encs.iter() {
        chains.extend(enumerate_chain_norev(&nodes_run, enc, true));
    }
    // Reverse
    let len = read.original_length;
    nodes_run.reverse();
    nodes_run.iter_mut().for_each(|x| x.reverse(len));
    for enc in encs.iter() {
        chains.extend(enumerate_chain_norev(&nodes_run, enc, false));
    }
    chains
}

#[derive(Debug, Clone)]
struct LightNode {
    node: (u64, u64),
    is_forward: bool,
    range: (usize, usize),
}

impl LightNode {
    fn new(node: (u64, u64), is_forward: bool, range: (usize, usize)) -> Self {
        Self {
            node,
            is_forward,
            range,
        }
    }
    fn reverse(&mut self, len: usize) {
        self.is_forward = !self.is_forward;
        self.range = (len - self.range.1, len - self.range.0);
    }
}

fn enumerate_chain_norev(nodes: &[LightNode], enc: &ContigEncoding, direction: bool) -> Vec<Chain> {
    let mut chain_nodes = vec![];
    for (q_idx, node) in nodes.iter().enumerate() {
        for (r_idx, target) in enc.matches(node.node, node.is_forward) {
            chain_nodes.push(ChainNode::new(q_idx, node, r_idx, target));
        }
    }
    // log::debug!("EnumChain\t{}", chain_nodes.len());
    if chain_nodes.is_empty() {
        return Vec::new();
    }
    // We obtain the topological order by sorting both start positions.
    chain_nodes.sort_by_key(|c| (c.contig_start, c.read_start));
    let mut chains: Vec<_> = vec![];
    while !chain_nodes.is_empty() {
        let chain_indices = min_chain(&chain_nodes);
        let first = chain_indices.first().map(|&i| chain_nodes[i]).unwrap();
        let last = chain_indices.last().map(|&i| chain_nodes[i]).unwrap();
        chains.push(align_in_chunk_space(nodes, enc, direction, first, last));
        let mut idx = 0;
        chain_nodes.retain(|_| {
            idx += 1;
            !chain_indices.contains(&(idx - 1))
        });
    }
    chains
}
const CHAIN_MATCH: i64 = -4_000;
fn min_chain(chain_nodes: &[ChainNode]) -> Vec<usize> {
    let sentinel = chain_nodes.len();
    let mut parents = vec![sentinel; chain_nodes.len()];
    // Initialize.
    let mut min_dist = vec![CHAIN_MATCH; chain_nodes.len()];
    for (idx, c_node) in chain_nodes.iter().enumerate() {
        let arg_min = chain_nodes
            .iter()
            .enumerate()
            .take(idx)
            .filter_map(|(i, other)| {
                let gap = other.to(c_node)?;
                let penalty = (gap as f64).ln().ceil() as i64;
                Some((i, min_dist[i] + penalty))
            })
            .min_by_key(|x| x.1);
        if let Some((parent, min)) = arg_min {
            let dist = min + CHAIN_MATCH;
            if dist < min_dist[idx] {
                min_dist[idx] = dist;
                parents[idx] = parent;
            }
        }
    }
    let (mut idx, _min) = min_dist.iter().enumerate().min_by_key(|x| x.1).unwrap();
    // log::debug!("MinChain\t{_min}");
    let mut traceback = vec![idx];
    while parents[idx] != sentinel {
        idx = parents[idx];
        traceback.push(idx);
    }
    traceback.reverse();
    traceback
}

fn align_in_chunk_space(
    nodes: &[LightNode],
    enc: &ContigEncoding,
    is_forward: bool,
    first: ChainNode,
    last: ChainNode,
) -> Chain {
    let query = &nodes[first.read_index..last.read_index + 1];
    let refr = &enc.tiles()[first.contig_index..last.contig_index + 1];
    let mut ops = alignment(query, refr);
    // if nodes.iter().any(|n| n.node.0 == 20) {
    //     debug!("{first:?}\t{ops:?}\t{:?}", &nodes[first.read_index]);
    // }
    // Removing head/tail...
    let (start_position, end_position) = get_range(first, &mut ops);
    Chain::new(
        enc.id.clone(),
        start_position,
        end_position,
        nodes.len(),
        is_forward,
        ops,
    )
}

const MIS_CLUSTER: i64 = 2;
const MISM: i64 = 100;
const GAP: i64 = 4;
const MATCH: i64 = -10;
fn match_score(n: &LightNode, m: &UnitAlignmentInfo) -> (i64, Op) {
    let ((unit, cluster), dir) = m.unit_and_dir_info();
    let (n_unit, n_cluster) = n.node;
    if dir != n.is_forward || unit != n_unit {
        (MISM, Op::Mismatch)
    } else if n_cluster != cluster {
        (MIS_CLUSTER, Op::Match)
    } else {
        assert_eq!((n_unit, n_cluster, n.is_forward), (unit, cluster, dir));
        (MATCH, Op::Match)
    }
}
fn alignment(query: &[LightNode], refr: &[UnitAlignmentInfo]) -> Vec<Op> {
    let mut dp = vec![vec![0; refr.len() + 1]; query.len() + 1];
    for (i, row) in dp.iter_mut().enumerate() {
        row[0] = i as i64 * GAP;
    }
    for j in 0..refr.len() + 1 {
        dp[0][j] = j as i64 * GAP;
    }
    for (i, n) in query.iter().enumerate() {
        let i = i + 1;
        for (j, unit) in refr.iter().enumerate() {
            let j = j + 1;
            let (mat, _) = match_score(n, unit);
            dp[i][j] = (dp[i - 1][j] + GAP)
                .min(dp[i][j - 1] + GAP)
                .min(dp[i - 1][j - 1] + mat);
        }
    }
    let mut ops = vec![];
    let (mut qpos, mut rpos) = (query.len(), refr.len());
    while 0 < qpos && 0 < rpos {
        let current = dp[qpos][rpos];
        let (mat, op) = match_score(&query[qpos - 1], &refr[rpos - 1]);
        if current == dp[qpos - 1][rpos] + GAP {
            qpos -= 1;
            ops.push(Op::Ins);
        } else if current == dp[qpos][rpos - 1] + GAP {
            rpos -= 1;
            ops.push(Op::Del);
        } else {
            assert_eq!(dp[qpos - 1][rpos - 1] + mat, current);
            rpos -= 1;
            qpos -= 1;
            ops.push(op);
        }
    }
    ops.extend(std::iter::repeat(Op::Del).take(rpos));
    ops.extend(std::iter::repeat(Op::Ins).take(qpos));
    ops.reverse();
    ops
}

fn get_range(start: ChainNode, ops: &mut Vec<Op>) -> ((usize, usize), (usize, usize)) {
    // Remove useless Ins/Del.
    while matches!(ops.last(), Some(Op::Del) | Some(Op::Ins)) {
        assert!(matches!(ops.pop(), Some(Op::Del) | Some(Op::Ins)))
    }
    // Seek to the beginning.
    let (mut q_pos, mut r_pos) = (start.read_index, start.contig_index);
    ops.reverse();
    loop {
        match ops.last() {
            Some(Op::Del) => {
                ops.pop();
                r_pos += 1;
            }
            Some(Op::Ins) => {
                ops.pop();
                q_pos += 1;
            }
            _ => break,
        }
    }
    ops.reverse();
    let start_pair = (q_pos, r_pos);
    // Consume.
    for op in ops.iter() {
        match &op {
            Op::Mismatch | Op::Match => {
                r_pos += 1;
                q_pos += 1;
            }
            Op::Ins => q_pos += 1,
            Op::Del => r_pos += 1,
        }
    }
    let end_pair = (q_pos, r_pos);
    (start_pair, end_pair)
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
struct ChainNode {
    contig_index: usize,
    contig_start: usize,
    read_start: usize,
    read_index: usize,
}

impl ChainNode {
    fn new(
        read_index: usize,
        node: &LightNode,
        contig_index: usize,
        contig_elm: &UnitAlignmentInfo,
    ) -> Self {
        Self {
            contig_index,
            contig_start: contig_elm.contig_range().0,
            read_start: node.range.0,
            read_index,
        }
    }
    fn to(&self, to: &Self) -> Option<i64> {
        (self.read_start < to.read_start && self.contig_start < to.contig_start)
            .then(|| (to.read_start - self.read_start + to.contig_start - self.contig_start) as i64)
    }
}

#[derive(Debug, Clone)]
struct Chain {
    id: String,
    contig_start_idx: usize,
    #[allow(dead_code)]
    contig_end_idx: usize,
    // If false, query indices are after rev-comped.
    is_forward: bool,
    query_start_idx: usize,
    query_end_idx: usize,
    query_nodes_len: usize,
    ops: Vec<Op>,
}

impl Chain {
    fn apporox_score(&self) -> i64 {
        // TODO: ????
        self.ops
            .iter()
            .map(|&op| match op {
                Op::Mismatch => -10,
                Op::Match => 10,
                Op::Ins => -5,
                Op::Del => -5,
            })
            .sum()
    }
    fn coord(&self) -> (usize, usize) {
        match self.is_forward {
            true => (self.query_start_idx, self.query_end_idx),
            false => (
                self.query_nodes_len - self.query_end_idx,
                self.query_nodes_len - self.query_start_idx,
            ),
        }
    }
    fn overlap_frac(&self, other: &Self) -> f64 {
        // let len = self.query_end_idx - self.query_start_idx;
        let (start, end) = self.coord();
        let (others, othere) = other.coord();
        let len = end - start;
        let ovlp = end.min(othere).saturating_sub(start.max(others));
        assert!(ovlp <= len);
        ovlp as f64 / len as f64
    }
    fn new(
        id: String,
        (query_start_idx, contig_start_idx): (usize, usize),
        (query_end_idx, contig_end_idx): (usize, usize),
        query_nodes_len: usize,
        is_forward: bool,
        ops: Vec<Op>,
    ) -> Self {
        Self {
            id,
            contig_start_idx,
            contig_end_idx,
            is_forward,
            ops,
            query_start_idx,
            query_end_idx,
            query_nodes_len,
        }
    }
}

fn base_pair_alignment(
    read: &EncodedRead,
    seq: &[u8],
    chain: &Chain,
    seg: &gfa::Segment,
    encs: &crate::assemble::ditch_graph::ContigEncoding,
) -> Alignment {
    let seg_seq = seg.sequence.as_ref().unwrap().as_bytes();
    let tiles = convert_into_tiles(read, chain, encs);
    let (mut query, mut ops, tip_len) = match chain.contig_start_idx {
        0 => align_tip(seq, seg, chain, tiles.first().unwrap()),
        _ => (vec![], vec![], 0),
    };
    for w in tiles.windows(2) {
        append_range(&mut query, &mut ops, &w[0]);
        let seg_start = w[0].ctg_end;
        let seg_end = w[1].ctg_start;
        assert!(seg_start <= seg_end);
        let seg_bet = &seg_seq[seg_start..seg_end];
        if chain.is_forward {
            let r_start = w[0].read_end;
            let r_end = w[1].read_start;
            assert!(r_start <= r_end);
            let read_bet = &seq[r_start..r_end];
            extend_between(&mut query, &mut ops, read_bet, seg_bet);
        } else {
            let r_start = w[1].read_end;
            let r_end = w[0].read_start;
            assert!(r_start <= r_end);
            let read_bet = bio_utils::revcmp(&seq[r_start..r_end]);
            extend_between(&mut query, &mut ops, &read_bet, seg_bet);
        }
    }
    append_range(&mut query, &mut ops, tiles.last().unwrap());
    let (tail, tail_ops, tail_len) = match chain.contig_end_idx == encs.tiles().len() {
        true => align_tail(seq, seg, chain, tiles.last().unwrap()),
        false => (Vec::new(), Vec::new(), 0),
    };
    query.extend(tail);
    ops.extend(tail_ops);
    let contig_start = tiles.first().map(|t| t.ctg_start).unwrap() - tip_len;
    let contig_end = tiles.last().map(|t| t.ctg_end).unwrap() + tail_len;
    if log_enabled!(log::Level::Trace) {
        let op_len: usize = ops.iter().filter(|&&op| op != Op::Ins).count();
        let len = contig_end - contig_start;
        assert_eq!(len, op_len, "R\t{}\t{}", contig_start, contig_end);
        let len = query.len();
        let op_len = ops.iter().filter(|&&op| op != Op::Del).count();
        assert_eq!(len, op_len, "Q");
        check(&query, seq, chain.is_forward);
    }
    Alignment {
        read_id: read.id,
        contig: encs.id.to_string(),
        contig_start,
        contig_end,
        query,
        ops,
        is_forward: chain.is_forward,
    }
}

fn check(query: &[u8], seq: &[u8], is_forward: bool) {
    let task = edlib_sys::AlignTask::Alignment;
    let mode = edlib_sys::AlignMode::Infix;
    assert!(query.iter().all(u8::is_ascii_uppercase));
    assert!(seq.iter().all(u8::is_ascii_uppercase));
    if is_forward {
        let contains = seq.windows(query.len()).any(|w| w == query);
        if !contains {
            let aln = edlib_sys::align(query, seq, mode, task);
            let (start, end) = aln.location().unwrap();
            let ops: Vec<_> = crate::misc::edlib_to_kiley(aln.operations().unwrap());
            let refr = &seq[start..end + 1];
            let (r, a, q) = kiley::recover(refr, query, &ops);
            for ((r, a), q) in r.chunks(200).zip(a.chunks(200)).zip(q.chunks(200)) {
                eprintln!("{}", std::str::from_utf8(r).unwrap());
                eprintln!("{}", std::str::from_utf8(a).unwrap());
                eprintln!("{}", std::str::from_utf8(q).unwrap());
            }
        }
        assert!(contains);
    } else {
        let rev = bio_utils::revcmp(seq);
        let contains = rev.windows(query.len()).any(|w| w == query);
        if !contains {
            let aln = edlib_sys::align(query, &rev, mode, task);
            let (start, end) = aln.location().unwrap();
            let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
            let refr = &rev[start..end + 1];
            let (r, a, q) = kiley::recover(refr, query, &ops);
            for ((r, a), q) in r.chunks(200).zip(a.chunks(200)).zip(q.chunks(200)) {
                eprintln!("{}", std::str::from_utf8(r).unwrap());
                eprintln!("{}", std::str::from_utf8(a).unwrap());
                eprintln!("{}", std::str::from_utf8(q).unwrap());
            }
        }
        assert!(contains);
    }
}

fn align_tail(
    seq: &[u8],
    seg: &gfa::Segment,
    chain: &Chain,
    tile: &Tile,
) -> (Vec<u8>, Vec<Op>, usize) {
    let seg_start = tile.ctg_end;
    let seg = seg.sequence.as_ref().unwrap().as_bytes();
    if seg_start == seg.len() {
        return (Vec::new(), Vec::new(), 0);
    }
    let seg_seq = &seg[seg_start..];
    let (query, ops) = if chain.is_forward {
        let read_start = tile.read_end;
        let read_end = (read_start + 2 * seg_seq.len()).min(seq.len());
        let trailing_seq = &seq[read_start..read_end];
        align_tail_inner(seg_seq, trailing_seq)
    } else {
        let read_end = tile.read_start;
        let read_start = read_end.saturating_sub(seg_seq.len() * 2);
        let trailing_seq = bio_utils::revcmp(&seq[read_start..read_end]);
        align_tail_inner(seg_seq, &trailing_seq)
    };
    let len_seg = ops.iter().filter(|&&op| op != Op::Ins).count();
    (query, ops, len_seg)
}

fn align_tip(
    seq: &[u8],
    seg: &gfa::Segment,
    chain: &Chain,
    tile: &Tile,
) -> (Vec<u8>, Vec<Op>, usize) {
    let seg_end = tile.ctg_start;
    if seg_end == 0 {
        return (vec![], vec![], 0);
    }
    // log::debug!("SegEnd\t{seg_end}");
    // TODO:Maybe we should constraint the length ... ?
    let seg_seq = &seg.sequence.as_ref().unwrap().as_bytes()[..seg_end];
    let (query, ops) = if chain.is_forward {
        let read_end = tile.read_start;
        let read_start = read_end.saturating_sub(2 * seg_end);
        // log::debug!("QUERY\t?\t{read_end}");
        let leading_seq = &seq[read_start..read_end];
        align_tip_inner(seg_seq, leading_seq)
    } else {
        let read_start = tile.read_end;
        let read_end = (read_start + 2 * seg_end).min(seq.len());
        // log::debug!("QUERY\t{read_start}\t?");
        let leading_seq = bio_utils::revcmp(&seq[read_start..read_end]);
        align_tip_inner(seg_seq, &leading_seq)
    };
    let len_seg = ops.iter().filter(|&&op| op != Op::Ins).count();
    (query, ops, len_seg)
}

fn align_tip_inner(seg: &[u8], query: &[u8]) -> (Vec<u8>, Vec<Op>) {
    // log::debug!("{}VS{}", seg.len(), query.len());
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    if seg.len() < query.len() {
        // Seg:         --------->
        // Que:  ------------------->
        //              |         |
        let aln = edlib_sys::align(seg, query, mode, task);
        let (start, end) = aln.location().unwrap();
        let end = end + 1;
        let ops_table = [Op::Match, Op::Del, Op::Ins, Op::Mismatch];
        let ops = aln.operations().unwrap();
        let mut ops: Vec<_> = ops.iter().map(|&op| ops_table[op as usize]).collect();
        assert!(end <= query.len());
        ops.extend(std::iter::repeat(Op::Ins).take(query.len() - end));
        let query = query[start..].to_vec();
        (query, ops)
    } else {
        // Seg: ------------------->
        // Query:          -------->
        let aln = edlib_sys::align(query, seg, mode, task);
        let (_, end) = aln.location().unwrap();
        let end = end + 1;
        let mut ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
        ops.extend(std::iter::repeat(Op::Del).take(seg.len() - end));
        // Pop leading insertions.
        ops.reverse();
        let mut start = 0;
        while ops.last() == Some(&Op::Ins) {
            start += ops.pop().is_some() as usize;
        }
        ops.reverse();
        let query = query[start..].to_vec();
        (query, ops)
    }
}

fn align_tail_inner(seg: &[u8], query: &[u8]) -> (Vec<u8>, Vec<Op>) {
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    if seg.len() < query.len() {
        // Seg:  --------->
        // Que:  ------------------->
        //         |      |
        let aln = edlib_sys::align(seg, query, mode, task);
        let (start, end) = aln.location().unwrap();
        let end = end + 1;
        let ops_table = [Op::Match, Op::Del, Op::Ins, Op::Mismatch];
        let aln = aln.operations().unwrap();
        let mut ops = vec![Op::Ins; start];
        ops.extend(aln.iter().map(|&op| ops_table[op as usize]));
        let query = query[..end].to_vec();
        (query, ops)
    } else {
        //          |      |
        // Seg:   ------------------->
        // Query:   -------|->
        let aln = edlib_sys::align(query, seg, mode, task);
        let (start, _) = aln.location().unwrap();
        let mut ops = vec![Op::Del; start];
        ops.extend(crate::misc::edlib_to_kiley(aln.operations().unwrap()));
        let mut poped = 0;
        while ops.last() == Some(&Op::Ins) {
            poped += ops.pop().is_some() as usize;
        }
        let query = query[..query.len() - poped].to_vec();
        (query, ops)
    }
}

fn extend_between(query: &mut Vec<u8>, ops: &mut Vec<Op>, read: &[u8], seg: &[u8]) {
    if seg.is_empty() {
        query.extend(read);
        ops.extend(std::iter::repeat(Op::Ins).take(read.len()));
    } else if read.is_empty() {
        ops.extend(std::iter::repeat(Op::Del).take(seg.len()));
    } else {
        let mode = edlib_sys::AlignMode::Global;
        let task = edlib_sys::AlignTask::Alignment;
        let aln = edlib_sys::align(read, seg, mode, task);
        query.extend(read);
        ops.extend(crate::misc::edlib_to_kiley(aln.operations().unwrap()));
    }
}

struct Tile<'a, 'b> {
    ctg_start: usize,
    ctg_end: usize,
    read_start: usize,
    read_end: usize,
    node: &'a definitions::Node,
    encoding: &'b UnitAlignmentInfo,
}

impl<'a, 'b> Tile<'a, 'b> {
    fn new(
        ctg_start: usize,
        ctg_end: usize,
        read_start: usize,
        read_end: usize,
        node: &'a definitions::Node,
        encoding: &'b UnitAlignmentInfo,
    ) -> Self {
        Self {
            ctg_start,
            ctg_end,
            read_start,
            read_end,
            node,
            encoding,
        }
    }
}

fn convert_into_tiles<'a, 'b>(
    read: &'a EncodedRead,
    chain: &Chain,
    encs: &'b ContigEncoding,
) -> Vec<Tile<'a, 'b>> {
    let mut tiles = vec![];
    let (mut read_idx, mut seg_idx) = (chain.query_start_idx, chain.contig_start_idx);
    let mut prev_r_pos = None;
    for &op in chain.ops.iter() {
        match op {
            Op::Match => {
                let read_node = match chain.is_forward {
                    true => &read.nodes[read_idx],
                    false => &read.nodes[read.nodes.len() - read_idx - 1],
                };
                let node_start = read_node.position_from_start;
                let seg_node = &encs.tiles()[seg_idx];
                let (start, end) = convert_to_read_coordinate(read_node, seg_node);
                let (read_start, read_end) = (node_start + start, node_start + end);
                let (read_start, read_end) = match (chain.is_forward, prev_r_pos) {
                    (_, None) => (read_start, read_end),
                    (true, Some((_, e))) => (read_start.max(e), read_end.max(e)),
                    (false, Some((s, _))) => (read_start.min(s), read_end.min(s)),
                };
                assert!(
                    node_start <= read_start,
                    "{},{},{},{:?},{}",
                    node_start,
                    read_start,
                    start,
                    prev_r_pos.unwrap(),
                    chain.is_forward
                );
                prev_r_pos = Some((read_start, read_end));
                let (seg_start, seg_end) = seg_node.contig_range();
                let tile = Tile::new(
                    seg_start, seg_end, read_start, read_end, read_node, seg_node,
                );
                tiles.push(tile);
                read_idx += 1;
                seg_idx += 1;
            }
            Op::Mismatch => {
                read_idx += 1;
                seg_idx += 1;
            }
            Op::Del => seg_idx += 1,
            Op::Ins => read_idx += 1,
        }
    }
    tiles
}

fn convert_to_read_coordinate(node: &definitions::Node, seg: &UnitAlignmentInfo) -> (usize, usize) {
    let map_range_in_unit = match seg.unit_range() {
        (true, start, end) => (start, end),
        (false, start, end) => (seg.unit_len() - end, seg.unit_len() - start),
    };
    let (start, end) = convert_to_map_range(node, map_range_in_unit);
    match node.is_forward {
        true => (start, end),
        false => (node.seq().len() - end, node.seq().len() - start),
    }
}

fn convert_to_map_range(node: &definitions::Node, (start, end): (usize, usize)) -> (usize, usize) {
    let (mut readpos, mut unitpos) = (0, 0);
    let mut read_start = 0;
    for op in node.cigar.iter() {
        match op {
            definitions::Op::Match(l) if end <= unitpos + l => {
                readpos += end - unitpos;
                break;
            }
            definitions::Op::Match(l) if start <= unitpos => {
                readpos += l;
                unitpos += l;
            }
            definitions::Op::Match(l) if start <= unitpos + l => {
                assert_eq!(read_start, 0);
                read_start = readpos + start - unitpos;
                readpos += l;
                unitpos += l;
            }
            definitions::Op::Match(l) => {
                readpos += l;
                unitpos += l;
            }
            definitions::Op::Del(l) if end <= unitpos + l => {
                readpos += end - unitpos;
                break;
            }
            definitions::Op::Del(l) if start <= unitpos => unitpos += l,
            definitions::Op::Del(l) if start <= unitpos + l => {
                assert_eq!(read_start, 0);
                read_start = readpos + start - unitpos;
                unitpos += l;
            }
            definitions::Op::Del(l) => unitpos += l,
            definitions::Op::Ins(l) => readpos += l,
        }
    }
    // .min is to cap the readpos. Sometimes the reference is too large...
    (read_start, readpos.min(node.seq().len()))
}

// TODO:Make this as an iterator, not allocating a vector!
fn defop2kileyop(ops: &definitions::Ops) -> Vec<Op> {
    ops.iter()
        .flat_map(|&op| match op {
            definitions::Op::Match(l) => std::iter::repeat(Op::Match).take(l),
            definitions::Op::Del(l) => std::iter::repeat(Op::Del).take(l),
            definitions::Op::Ins(l) => std::iter::repeat(Op::Ins).take(l),
        })
        .collect()
}

fn append_range(query: &mut Vec<u8>, ops: &mut Vec<Op>, tile: &Tile) {
    let (r_start, r_end) = (tile.read_start, tile.read_end);
    let offset = tile.node.position_from_start;
    assert!(offset <= r_start);
    assert!(r_start <= r_end);
    // (r_start - offset, r_end - offset);
    let seqlen = tile.node.seq().len();
    let (r_start, r_end) = match tile.node.is_forward {
        true => (r_start - offset, r_end - offset),
        false => (seqlen + offset - r_end, seqlen + offset - r_start),
    };
    assert!(r_end <= tile.node.seq().len());
    let (start, end) = match tile.encoding.unit_range() {
        (true, start, end) => (start, end),
        (false, start, end) => (
            tile.encoding.unit_len() - end,
            tile.encoding.unit_len() - start,
        ),
    };
    let alignment = defop2kileyop(&tile.node.cigar);
    let r_len = alignment.iter().filter(|&&op| op != Op::Del).count();
    let u_len = alignment.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(u_len, tile.encoding.unit_len());
    assert_eq!(r_len, tile.node.seq().len());
    let (mut r_pos, mut u_pos) = (0, 0);
    // log::debug!("SEG\t{start}-{end}");
    // log::debug!("QUE\t{r_start}-{r_end}");
    let mut temp_o = vec![];
    for op in alignment {
        if (start..end).contains(&u_pos) && (r_start..r_end).contains(&r_pos) {
            temp_o.push(op);
        } else if (start..end).contains(&u_pos) && op != Op::Ins {
            temp_o.push(Op::Del);
        } else if (r_start..r_end).contains(&r_pos) && op != Op::Del {
            temp_o.push(Op::Ins)
        }
        match op {
            Op::Match | Op::Mismatch => {
                u_pos += 1;
                r_pos += 1;
            }
            Op::Ins => r_pos += 1,
            Op::Del => u_pos += 1,
        }
    }
    let aligned_seq = &tile.node.seq()[r_start..r_end];
    let r_len = temp_o.iter().filter(|&&op| op != Op::Del).count();
    let u_len = temp_o.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(r_len, r_end - r_start);
    assert_eq!(u_len, end - start);
    if tile.encoding.unit_range().0 {
        query.extend(aligned_seq);
        ops.extend(temp_o);
    } else {
        query.extend(crate::seq::DNAIter::new(aligned_seq, false));
        temp_o.reverse();
        ops.extend(temp_o);
    }
}

#[derive(Debug, Clone)]
pub struct Alignment {
    #[allow(dead_code)]
    read_id: u64,
    contig: String,
    contig_start: usize,
    contig_end: usize,
    query: Vec<u8>,
    ops: Vec<Op>,
    #[allow(dead_code)]
    is_forward: bool,
}

impl Alignment {}

#[cfg(test)]
mod align_test {
    use super::*;
    fn mock_chain_node(rstart: usize, cstart: usize) -> ChainNode {
        ChainNode {
            contig_index: 0,
            contig_start: cstart,
            read_start: rstart,
            read_index: 0,
        }
    }
    #[test]
    fn to_test() {
        let s1 = mock_chain_node(0, 0);
        let s2 = mock_chain_node(1, 1);
        let s3 = mock_chain_node(1, 2);
        let s4 = mock_chain_node(2, 1);
        let s5 = mock_chain_node(2, 2);
        assert_eq!(s1.to(&s2), Some(2));
        assert_eq!(s3.to(&s1), None);
        assert_eq!(s4.to(&s2), None);
        assert_eq!(s2.to(&s5), Some(2));
        assert_eq!(s5.to(&s2), None);
    }
    #[test]
    fn seek_ops() {
        use Op::*;
        let mut ops = vec![Del, Ins, Match, Ins, Del, Mismatch, Match, Ins, Del];
        let start = ChainNode {
            contig_index: 10,
            contig_start: 0,
            read_start: 0,
            read_index: 1,
        };
        let (start_range, end_range) = get_range(start, &mut ops);
        assert_eq!(ops, vec![Match, Ins, Del, Mismatch, Match]);
        assert_eq!(start_range, (2, 11));
        assert_eq!(end_range, (6, 15));
    }
    #[test]
    fn test_alignment() {
        let range = (0, 0);
        let query: Vec<_> = vec![
            LightNode::new((0, 0), true, range),
            LightNode::new((1, 0), true, range),
            LightNode::new((2, 0), true, range),
            LightNode::new((3, 1), true, range),
            LightNode::new((4, 0), true, range),
        ];
        let refr = vec![
            unit_aln_info((7, 0), true),
            unit_aln_info((1, 0), true),
            unit_aln_info((2, 1), true),
            unit_aln_info((3, 1), true),
            unit_aln_info((6, 0), true),
        ];
        let ops = alignment(&query, &refr);
        use Op::*;
        assert_eq!(ops, vec![Del, Ins, Match, Match, Match, Del, Ins]);
    }
    fn unit_aln_info((unit, cluster): (u64, u64), direction: bool) -> UnitAlignmentInfo {
        UnitAlignmentInfo::new((unit, cluster), (direction, 0, 0), (0, 0), 0)
    }
    fn chain_node(read_start: usize, contig_start: usize) -> ChainNode {
        ChainNode {
            contig_index: 0,
            contig_start,
            read_start,
            read_index: 0,
        }
    }
    #[test]
    fn min_chain_test() {
        let chain_nodes = vec![
            chain_node(0, 0),
            chain_node(1, 1),
            chain_node(1, 3),
            chain_node(2, 2),
        ];
        let answer = vec![0, 1, 3];
        let indices = min_chain(&chain_nodes);
        assert_eq!(indices, answer);
        let chain_nodes = vec![
            chain_node(0, 0),  // 0
            chain_node(0, 5),  // 1
            chain_node(1, 1),  // 2
            chain_node(1, 6),  // 3
            chain_node(1, 10), // 4
            chain_node(2, 2),  // 5
            chain_node(2, 10), // 6
            chain_node(3, 11), //7
        ];
        let answer = vec![1, 3, 6, 7];
        let indices = min_chain(&chain_nodes);
        assert_eq!(indices, answer);
    }
    #[test]
    fn leading_aln_test() {
        let seq = b"AAAAA";
        let leading = b"TTTTTTTTAAAAAC";
        let (ops, len) = align_leading(seq, leading);
        assert_eq!(ops, vec![vec![Op::Match; 5], vec![Op::Del]].concat());
        assert_eq!(len, 6);
    }
    #[test]
    fn trailing_aln_test() {
        let seq = b"AAAAA";
        let trailing = b"GAAAAATTTTTTTTTTTTC";
        let (ops, len) = align_trailing(seq, trailing);
        assert_eq!(ops, vec![vec![Op::Mismatch], vec![Op::Match; 4]].concat());
        assert_eq!(len, 5);
    }
}
