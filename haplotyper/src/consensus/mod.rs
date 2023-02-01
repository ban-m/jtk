use crate::assemble::ditch_graph::ContigEncoding;
use crate::assemble::ditch_graph::UnitAlignmentInfo;
use crate::model_tune::ModelFit;
use definitions::*;
use gfa::Segment;
use kiley::hmm::PairHiddenMarkovModelOnStrands;
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
        dump_path: &Option<String>,
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
    max_coverage: usize,
    window_size: usize,
    radius: usize,
    round_num: usize,
}
impl PolishConfig {
    pub fn new(
        seed: u64,
        min_coverage: usize,
        max_coverage: usize,
        window_size: usize,
        radius: usize,
        round_num: usize,
    ) -> PolishConfig {
        Self {
            seed,
            max_coverage,
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
            max_coverage: 40,
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
        dump_path: &Option<String>,
    ) -> Vec<Segment> {
        let mut alignments_on_contigs = self.distribute_to_contig(segments, encs, config);
        let models = self.fit_models_on_both_strands().unwrap();
        let polished: Vec<_> = alignments_on_contigs
            .iter_mut()
            .map(|(sid, alignments)| {
                log::debug!("POLISH\t{sid}");
                let seg = segments.iter().find(|seg| &seg.sid == sid).unwrap();
                let seg = seg.sequence.as_ref().unwrap().as_bytes();
                let seq = polish(sid, seg, alignments, &models, config);
                (sid.clone(), seq)
            })
            .map(|(sid, seq)| {
                let slen = seq.len();
                let sequence = String::from_utf8(seq).ok();
                Segment::from(sid, slen, sequence)
            })
            .collect();
        if let Some(path) = dump_path.as_ref() {
            if let Err(why) = dump_alignments(self, &alignments_on_contigs, &polished, path) {
                warn!("{why:?}");
            }
        }
        polished
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
            .flat_map(|(i, read)| {
                let seed = config.seed + i as u64;
                let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
                let mut raw_read = read.recover_raw_read();
                raw_read.iter_mut().for_each(u8::make_ascii_uppercase);
                align_to_contigs(read, &raw_read, encs, segments, &mut rng)
            })
            .collect();
        let mut alignments_on_contigs: BTreeMap<_, Vec<_>> = BTreeMap::new();
        for aln in alignments {
            alignments_on_contigs
                .entry(aln.contig.clone())
                .or_default()
                .push(aln);
        }
        alignments_on_contigs
    }
}

fn dump_alignments(
    ds: &DataSet,
    alignments_on_contigs: &BTreeMap<String, Vec<Alignment>>,
    polished: &[Segment],
    dump_path: &str,
) -> std::io::Result<()> {
    use std::io::*;
    let mut sam_file = std::fs::File::create(format!("{dump_path}.sam")).map(BufWriter::new)?;
    dump_sam_header(polished, &mut sam_file)?;
    dump_sam_records(ds, alignments_on_contigs, &mut sam_file)?;
    let mut coverage_file =
        std::fs::File::create(format!("{dump_path}.coverage.tsv")).map(BufWriter::new)?;
    dump_coverages(polished, alignments_on_contigs, &mut coverage_file)?;
    Ok(())
}

fn dump_sam_header<W: std::io::Write>(polished: &[Segment], wtr: &mut W) -> std::io::Result<()> {
    writeln!(wtr, "@HD\tVN:1.6")?;
    for seg in polished {
        writeln!(wtr, "@SQ\tSN:{}\tLN:{}", seg.sid, seg.slen)?;
    }
    writeln!(wtr, "@PG\tID:jtk\tPN:jtk\tVN:0.1")?;
    Ok(())
}
fn dump_sam_records<W: std::io::Write>(
    ds: &DataSet,
    alignments_on_contigs: &BTreeMap<String, Vec<Alignment>>,
    wtr: &mut W,
) -> std::io::Result<()> {
    use std::collections::HashSet;
    let mut is_secondary = HashSet::new();
    let raw_reads: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r)).collect();
    for (sid, alignments) in alignments_on_contigs {
        for aln in alignments {
            let read = match raw_reads.get(&aln.read_id) {
                Some(res) => res,
                _ => {
                    warn!("NOALN\t{sid},{}", aln.read_id,);
                    continue;
                }
            };
            let secondary = is_secondary.contains(&aln.read_id);
            is_secondary.insert(aln.read_id);
            writeln!(wtr, "{}", sam_record(sid, aln, read, secondary))?;
        }
    }
    Ok(())
}
fn sam_record(sid: &str, aln: &Alignment, read: &definitions::RawRead, secondary: bool) -> String {
    let mut records = String::new();
    records.push_str(&read.name);
    records.push('\t');
    let flag = match aln.is_forward {
        true => 0,
        false => 0x10,
    };
    let flag = match secondary {
        true => flag + 0x800,
        false => flag,
    };
    records.push_str(&format!("{flag}\t"));
    records.push_str(sid);
    records.push('\t');
    records.push_str(&format!("{}\t", aln.contig_start + 1));
    records.push_str("60\t");
    records.push_str(&aln.cigar());
    records.push('\t');
    records.push_str("*\t0\t0\t");
    records.push_str(std::str::from_utf8(&aln.query).unwrap());
    records.push_str("\t*");
    records
}

const SMOOTH_WINDOW: usize = 1_000;
fn dump_coverages<W: std::io::Write>(
    polished: &[Segment],
    alignments_on_contigs: &BTreeMap<String, Vec<Alignment>>,
    wtr: &mut W,
) -> std::io::Result<()> {
    let mut coverages: HashMap<_, _> = polished
        .iter()
        .map(|seg| {
            let len = (seg.slen as usize) / SMOOTH_WINDOW + 1;
            (seg.sid.clone(), vec![0; len])
        })
        .collect();
    for (sid, alns) in alignments_on_contigs.iter() {
        let coverage = coverages.get_mut(sid).unwrap();
        for aln in alns.iter() {
            let mut rpos = aln.contig_start;
            for op in aln.ops.iter() {
                match op {
                    Op::Mismatch | Op::Match | Op::Del => {
                        rpos += 1;
                        coverage[rpos / SMOOTH_WINDOW] += 1;
                    }
                    Op::Ins => {}
                }
            }
        }
    }
    writeln!(wtr, "contig\tposition\tcoverage")?;
    for (sid, coverages) in coverages.iter() {
        for (i, &cov) in coverages.iter().enumerate() {
            let pos = i * SMOOTH_WINDOW;
            let cov = cov as f64 / SMOOTH_WINDOW as f64;
            writeln!(wtr, "{sid}\t{pos}\t{cov:.3}")?
        }
    }
    Ok(())
}

fn log_identity(sid: &str, alignments: &[Alignment]) {
    if log_enabled!(log::Level::Debug) {
        let (mut mat, mut len) = (0, 0);
        for aln in alignments.iter() {
            let match_num = aln.ops.iter().filter(|&&op| op == Op::Match).count();
            let aln_len = aln.ops.len();
            mat += match_num;
            len += aln_len;
        }
        log::debug!("ALN\t{sid}\t{mat}\t{len}\t{}", mat as f64 / len as f64);
    }
}

// (alignment index, node index, direction, operations, seq)
type SeqOps<'a> = (usize, usize, bool, &'a [u8], Vec<Op>);
type SeqOpsOnWindow<'a> = Vec<SeqOps<'a>>;
type UsedRange = (TipPos, TipPos);
use std::collections::HashMap;
fn allocate_on_windows(
    alignments: &[Alignment],
    round: usize,
    window: usize,
    draft_len: usize,
) -> (Vec<SeqOpsOnWindow>, HashMap<usize, UsedRange>) {
    let num_slot = draft_len / window + (draft_len % window != 0) as usize;
    let mut slots = vec![vec![]; num_slot];
    let used_range: HashMap<_, _> = alignments
        .iter()
        .enumerate()
        .map(|(aln_idx, aln)| {
            let (start, chunks, end) = split(aln, window, num_slot, draft_len);
            if chunks.is_empty() {
                assert!(start.0 == end.0 || start.0 + 1 == end.0);
            }
            for (idx, (pos, seq, ops)) in chunks.into_iter().enumerate() {
                slots[pos].push((aln_idx, idx, aln.is_forward, seq, ops));
            }
            (aln_idx, (start, end))
        })
        .collect();
    if round != 0 {
        for pileup in slots.iter_mut() {
            pileup.sort_by_cached_key(|x| x.4.iter().filter(|&&op| op != Op::Match).count());
        }
    }
    (slots, used_range)
}

pub fn polish(
    sid: &str,
    draft: &[u8],
    alignments: &mut [Alignment],
    models: &PairHiddenMarkovModelOnStrands,
    config: &PolishConfig,
) -> Vec<u8> {
    let window = config.window_size;
    debug!("POLISH\tWINDOW\t{window}");
    let mut polished = draft.to_vec();
    log_identity(sid, alignments);
    let min_coverage = config.min_coverage;
    for round in 0..config.round_num {
        let to_refresh = 0 == round;
        let (mut pileups, used_ranges) =
            allocate_on_windows(alignments, round, window, polished.len());
        let polished_seg: Vec<_> = polished
            .par_chunks(window)
            .zip(pileups.par_iter_mut())
            .map(|(draft, pileup)| {
                let (mut seqs, mut opss, mut strands) = (vec![], vec![], vec![]);
                for (_, _, strand, seq, ops) in pileup.iter_mut() {
                    strands.push(*strand);
                    seqs.push(*seq);
                    opss.push(ops);
                }
                match seqs.len() < min_coverage {
                    true => draft.to_vec(),
                    false => polish_seg(models, draft, &strands, &seqs, opss, config, to_refresh),
                }
            })
            .collect();
        let (acc_len, _) =
            polished_seg
                .iter()
                .enumerate()
                .fold((vec![0], 0), |(mut xs, len), (i, x)| {
                    if window * 3 < x.len() {
                        debug!("{i}\t{len}\t{}", x.len());
                    }
                    let len = len + x.len();
                    xs.push(len);
                    (xs, len)
                });
        assert_eq!(acc_len.len(), polished_seg.len() + 1);
        polished = polished_seg.iter().flatten().copied().collect();
        // Fix alignment.
        let mut recovered: HashMap<_, Vec<_>> = HashMap::new();
        for pileup in pileups {
            for (idx, pos, _, _, ops) in pileup {
                recovered.entry(idx).or_default().push((pos, ops));
            }
        }
        recovered
            .values_mut()
            .for_each(|xs| xs.sort_by_key(|x| x.0));
        assert_eq!(alignments.len(), used_ranges.len());
        alignments
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, aln)| {
                let used_range = used_ranges[&idx];
                match recovered.get(&idx) {
                    Some(aloc_pos) => fix_alignment(aln, used_range, aloc_pos, &polished, &acc_len),
                    None => fix_alignment(aln, used_range, &[], &polished, &acc_len),
                }
            });
        log_identity(sid, alignments);
    }
    polished
}

// #[allow(dead_code)]
// fn get_threshold(pileups: &[SeqOpsOnWindow]) -> f64 {
//     let identity_rates: Vec<_> = pileups
//         .iter()
//         .filter(|pileup| !pileup.is_empty())
//         .map(|pileup| {
//             let matches: f64 = pileup
//                 .iter()
//                 .map(|(_, _, _, ops)| {
//                     let mat = ops.iter().filter(|&&op| op == Op::Match).count();
//                     mat as f64 / ops.len() as f64
//                 })
//                 .sum();
//             matches / pileup.len() as f64
//         })
//         .collect();
//     let identity_rate = identity_rates.iter().sum::<f64>() / identity_rates.len() as f64;
//     let var = identity_rates
//         .iter()
//         .map(|e| (e - identity_rate).powi(2))
//         .sum::<f64>()
//         / identity_rates.len() as f64;
//     let sd = var.sqrt().min(0.02);
//     let thr = identity_rate - 5f64 * sd;
//     debug!("TRAIN\t{identity_rate:.3}\t{sd:.3}\t{thr:.3}");
//     thr
// }

// #[allow(dead_code)]
// fn train_hmm_dev(
//     hmm: &mut PairHiddenMarkovModel,
//     draft: &[u8],
//     window: usize,
//     pileups: &[SeqOpsOnWindow],
// ) {
//     let identity_threshold = get_threshold(pileups);
//     let mut transitions = [[1f64; 3]; 3];
//     let mut mat_emit = [1f64; 16];
//     let mut ins_emit = [1f64; 20];
//     for (template, pileup) in draft.chunks_exact(window).zip(pileups.iter()) {
//         let filtered = pileup.iter().filter_map(|(_, _, seq, ops)| {
//             let mat = ops.iter().filter(|&&op| op == Op::Match).count();
//             (identity_threshold < (mat as f64 / ops.len() as f64)).then_some((seq, ops.as_slice()))
//         });
//         for (seq, ops) in filtered {
//             let params = (&mut transitions, &mut mat_emit, &mut ins_emit);
//             crate::model_tune::register_alignments(template, seq, ops, params);
//         }
//     }
//     let mat = (transitions[0][0], transitions[0][1], transitions[0][2]);
//     let ins = (transitions[1][0], transitions[1][1], transitions[1][2]);
//     let del = (transitions[2][0], transitions[2][1], transitions[2][2]);
//     *hmm = PairHiddenMarkovModel::new(mat, ins, del, &mat_emit, &ins_emit);
//     debug!("TRAIN\n{hmm}");
// }

// fn train_hmm(
//     hmm: &mut PairHiddenMarkovModel,
//     draft: &[u8],
//     window: usize,
//     pileups: &[SeqOpsOnWindow],
//     radius: usize,
// ) {
//     let mut coverages: Vec<_> = pileups.iter().map(|x| x.len()).collect();
//     if coverages.is_empty() {
//         return;
//     }
//     let idx = coverages.len() / 2;
//     let (_, &mut med_cov, _) = coverages.select_nth_unstable(idx);
//     let range = 2 * window / 3..4 * window / 3;
//     let cov_range = 2 * med_cov / 3..4 * med_cov / 3;
//     let iterator = draft
//         .chunks(window)
//         .zip(pileups.iter())
//         .map(|(template, pileup)| {
//             let (seqs, ops): (Vec<&[u8]>, Vec<&[Op]>) = pileup
//                 .iter()
//                 .filter_map(|(_, _, _, seq, ops)| {
//                     let indel = ops.iter().map(|op| match op {
//                         Op::Mismatch => 1,
//                         Op::Match => -1,
//                         Op::Ins => 1,
//                         Op::Del => 1,
//                     });
//                     if crate::misc::max_region(indel) < 30 {
//                         Some((seq, ops.as_slice()))
//                     } else {
//                         None
//                     }
//                 })
//                 .unzip();
//             (template, seqs, ops)
//         })
//         .filter(|(draft, _, _)| range.contains(&draft.len()))
//         .filter(|(_, _, ops)| cov_range.contains(&ops.len()))
//         .take(3);
//     for (template, seqs, ops) in iterator {
//         hmm.fit_naive_with_par(template, &seqs, ops.as_slice(), radius);
//     }
// }

fn length_median(seqs: &[&[u8]]) -> usize {
    let mut len: Vec<_> = seqs.iter().map(|x| x.len()).collect();
    let idx = seqs.len() / 2;
    *len.select_nth_unstable(idx).1
}

fn within_range(median: usize, len: usize, frac: f64) -> bool {
    let diff = median.max(len) - median.min(len);
    (diff as f64 / median as f64) < frac
}

fn remove_refernece(ops: Vec<&mut Vec<Op>>, seqs: &[&[u8]]) {
    for (seq, op) in std::iter::zip(seqs, ops) {
        let qlen = op.iter().filter(|&&op| op != Op::Del).count();
        assert_eq!(seq.len(), qlen);
        op.clear();
        op.extend(std::iter::repeat(Op::Ins).take(qlen));
    }
}

type SplitQuery<'a, 'b> = (
    Vec<bool>,
    Vec<&'a [u8]>,
    Vec<&'b mut Vec<Op>>,
    Vec<(usize, &'b mut Vec<Op>)>,
);
fn split_sequences<'a, 'b>(
    strands: &[bool],
    seqs: &[&'a [u8]],
    ops: Vec<&'b mut Vec<Op>>,
    length_median: usize,
) -> SplitQuery<'a, 'b> {
    let mut will_be_use = vec![];
    let mut update_indices = vec![];
    assert_eq!(seqs.len(), ops.len());
    for (idx, ((seq, ops), strand)) in seqs.iter().zip(ops).zip(strands).enumerate() {
        if within_range(length_median, seq.len(), 0.15) {
            will_be_use.push((*seq, ops, *strand));
        } else {
            update_indices.push((idx, ops));
        }
    }
    let (mut use_seqs, mut use_ops, mut use_strands) = (vec![], vec![], vec![]);
    for (seq, ops, strand) in will_be_use {
        use_seqs.push(seq);
        use_ops.push(ops);
        use_strands.push(strand);
    }
    // let (use_seqs, use_ops): (Vec<_>, Vec<_>) = will_be_use.into_iter().unzip();
    (use_strands, use_seqs, use_ops, update_indices)
}

fn global_align(query: &[u8], target: &[u8]) -> Vec<Op> {
    if query.is_empty() {
        vec![Op::Del; target.len()]
    } else if target.is_empty() {
        vec![Op::Ins; query.len()]
    } else {
        let mode = edlib_sys::AlignMode::Global;
        let task = edlib_sys::AlignTask::Alignment;
        let align = edlib_sys::align(query, target, mode, task);
        crate::misc::edlib_to_kiley(align.operations().unwrap())
    }
}

fn bootstrap_consensus(seqs: &[&[u8]], ops: &mut [Vec<Op>], radius: usize) -> Vec<u8> {
    // let draft = kiley::ternary_consensus_by_chunk(seqs, radius);
    let draft = seqs[0].to_vec();
    for (seq, ops) in std::iter::zip(seqs, ops.iter_mut()) {
        *ops = global_align(seq, &draft);
    }
    kiley::bialignment::guided::polish_until_converge_with(&draft, seqs, ops, radius)
}

fn polish_seg(
    models: &PairHiddenMarkovModelOnStrands,
    draft: &[u8],
    strands: &[bool],
    seqs: &[&[u8]],
    ops: Vec<&mut Vec<Op>>,
    config: &PolishConfig,
    refresh: bool,
) -> Vec<u8> {
    let radius = config.radius;
    let max_cov = config.max_coverage;
    let length_median = length_median(seqs);
    if length_median == 0 {
        remove_refernece(ops, seqs);
        return Vec::new();
    }
    for ops in ops.iter() {
        let reflen = ops.iter().filter(|&&op| op != Op::Ins).count();
        assert_eq!(reflen, draft.len());
    }
    let (use_strands, use_seqs, use_ops, update_indices) =
        split_sequences(strands, seqs, ops, length_median);
    assert_eq!(seqs.len(), use_seqs.len() + update_indices.len());
    let mut polished = draft.to_vec();
    let mut temp_ops: Vec<_> = use_ops.iter().map(|x| x.to_vec()).collect();
    use kiley::bialignment::guided::polish_until_converge_with;
    polished = match within_range(draft.len(), length_median, 0.2) {
        true if refresh => polish_until_converge_with(&polished, &use_seqs, &mut temp_ops, radius),
        true => polished,
        false => bootstrap_consensus(&use_seqs, &mut temp_ops, radius),
    };
    let config = kiley::hmm::guided::HMMConfig::new(radius, max_cov, 0);
    polished = models.polish_until_converge_with_conf(
        &polished,
        &use_seqs,
        &mut temp_ops,
        &use_strands,
        &config,
    );
    assert_eq!(temp_ops.len(), use_ops.len());
    for (old, new) in std::iter::zip(use_ops, temp_ops) {
        *old = new;
        let reflen = old.iter().filter(|&&op| op != Op::Ins).count();
        assert_eq!(reflen, polished.len());
    }
    for (idx, ops) in update_indices {
        *ops = global_align(seqs[idx], &polished);
        let reflen = ops.iter().filter(|&&op| op != Op::Ins).count();
        assert_eq!(reflen, polished.len());
    }
    polished
}

fn fix_alignment(
    aln: &mut Alignment,
    used_range: UsedRange,
    aloc_ops: &[(usize, Vec<Op>)],
    polished: &[u8],
    acc_len: &[usize],
) {
    aln.ops.clear();
    let ((start, first_pos), (end, last_pos)) = used_range;
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
        // let first_pos_bp = acc_len[first_pos + 1];
        let first_pos_bp = acc_len[first_pos];
        let (ops, contig_len) = align_leading(&aln.query[..start], &polished[..first_pos_bp]);
        aln.contig_start = first_pos_bp - contig_len;
        aln.ops.extend(ops);
    } else {
        aln.contig_start = acc_len[first_pos];
    }
    aln.ops.extend(aloc_ops.iter().flat_map(|(_, ops)| ops));
    if end != aln.query.len() {
        let last_pos_bp = acc_len[last_pos];
        let (ops, contig_len) = align_trailing(&aln.query[end..], &polished[last_pos_bp..]);
        aln.ops.extend(ops);
        aln.contig_end = last_pos_bp + contig_len;
    } else {
        aln.contig_end = acc_len[last_pos];
    }
    while aln.ops.last() == Some(&Op::Del) {
        aln.ops.pop();
        aln.contig_end -= 1;
    }
    aln.ops.reverse();
    while aln.ops.last() == Some(&Op::Del) {
        aln.ops.pop();
        aln.contig_start += 1;
    }
    aln.ops.reverse();
    let reflen = aln.ops.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(
        reflen,
        aln.contig_end - aln.contig_start,
        "{},{},{},{},{},{},{},{}",
        start,
        first_pos,
        end,
        last_pos,
        aln.query.len(),
        aln.contig_start,
        aln.contig_end,
        acc_len[first_pos]
    );
    let refr = &polished[aln.contig_start..aln.contig_end];
    use kiley::bialignment::guided;
    aln.ops = guided::global_guided(refr, &aln.query, &aln.ops, 10, crate::ALN_PARAMETER).1;
}

fn align_infix(query: &[u8], seg: &[u8]) -> (usize, Vec<Op>, usize) {
    assert!(!query.is_empty());
    assert!(!seg.is_empty());
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
    assert_eq!(len, alignment.query.len(), "{:?}", alignment);
    let refr = alignment.ops.iter().filter(|&&op| op != Op::Ins).count();
    let ref_len = alignment.contig_end - alignment.contig_start;
    assert_eq!(
        refr, ref_len,
        "{},{}",
        alignment.is_forward, alignment.read_id
    );
    //    let start_chunk_id = alignment.contig_start / window;
    let start_pos_in_contig = match alignment.contig_start % window == 0 {
        true => (alignment.contig_start / window) * window,
        false => (alignment.contig_start / window + 1) * window,
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
    let start_chunk_id = cpos / window;
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
                    assert_eq!(qpos, alignment.query.len());
                    chunk_ops.extend(std::iter::repeat(Op::Del).take(contig_len - cpos));
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
            .map(|chain| chain.score)
            // .map(|chain| chain.apporox_score())
            .max()
            .unwrap();
        //                ((chains[idx].apporox_score() - max) as f64).exp()
        let picked = choises
            .choose_weighted(rng, |&idx| ((chains[idx].score - max) as f64).exp())
            .unwrap();
        let chain = chains.remove(*picked);
        chains.retain(|c| c.overlap_frac(&chain) < 0.5);
        let seg = segs.iter().find(|seg| seg.sid == chain.id).unwrap();
        let enc = encs.iter().find(|enc| enc.id == chain.id).unwrap();
        alns.push(base_pair_alignment(read, seq, &chain, seg, enc, alns.len()));
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
    for c in chains.iter() {
        assert_eq!(c.query_nodes_len, read.nodes.len());
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
        let chain = align_in_chunk_space(nodes, enc, direction, first, last);
        if !chain.ops.is_empty() {
            chains.push(chain);
        }
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
    let (mut ops, score) = alignment(query, refr);
    let (start_position, end_position) = get_range(first, &mut ops);
    Chain::new(
        enc.id.clone(),
        start_position,
        end_position,
        nodes.len(),
        is_forward,
        ops,
        score,
    )
}

const MIS_CLUSTER: i64 = -2;
const MISM: i64 = -100;
const GAP: i64 = -4;
const MATCH: i64 = 10;
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
// TODO: this function should return alignment score.
// Note that the posterior probabiltiy of the node is not available anymore
// but if we have that, this program would be much correct. Also, the copy number of the clsuter should be considered as well....
fn alignment(query: &[LightNode], refr: &[UnitAlignmentInfo]) -> (Vec<Op>, i64) {
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
                .max(dp[i][j - 1] + GAP)
                .max(dp[i - 1][j - 1] + mat);
        }
    }
    let mut ops = vec![];
    let (mut qpos, mut rpos) = (query.len(), refr.len());
    let score = dp[qpos][rpos];
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
    (ops, score)
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
        if self.read_start < to.read_start && self.contig_start < to.contig_start {
            Some((to.read_start - self.read_start + to.contig_start - self.contig_start) as i64)
        } else {
            None
        }
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
    score: i64,
}

impl Chain {
    // fn apporox_score(&self) -> i64 {
    //     self.ops
    //         .iter()
    //         .map(|&op| match op {
    //             Op::Mismatch => -10,
    //             Op::Match => 10,
    //             Op::Ins => -5,
    //             Op::Del => -5,
    //         })
    //         .sum()
    // }
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
        score: i64,
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
            score,
        }
    }
}

fn base_pair_alignment(
    read: &EncodedRead,
    seq: &[u8],
    chain: &Chain,
    seg: &gfa::Segment,
    encs: &crate::assemble::ditch_graph::ContigEncoding,
    _id: usize,
) -> Alignment {
    let seg_seq = seg.sequence.as_ref().unwrap().as_bytes();
    let tiles = convert_into_tiles(read, chain, encs);
    if tiles.is_empty() {
        eprintln!("{chain:?}\n{read}\n{}", read.nodes.len());
    }
    let (mut query, mut ops, tip_len, head_clip) = align_tip(seq, seg, chain, &tiles[0], read.id);
    for (_, w) in tiles.windows(2).enumerate() {
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
    let (tail, tail_ops, tail_len) = align_tail(seq, seg, chain, tiles.last().unwrap());
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
    let contig_range = (contig_start, contig_end);
    let id = encs.id.to_string();
    let tail_clip = seq.len() - head_clip - query.len();
    let clips = (head_clip, tail_clip);
    let is_forward = chain.is_forward;
    let contig_len = contig_end - contig_start;
    let reflen = ops.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(contig_len, reflen);
    Alignment::new(read.id, id, contig_range, clips, query, ops, is_forward)
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
            let (r, a, q) = kiley::op::recover(refr, query, &ops);
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
            let (r, a, q) = kiley::op::recover(refr, query, &ops);
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
    let (query, ops) = if chain.is_forward {
        let read_start = tile.read_end;
        let trailing_seq = &seq[read_start..];
        let seg_end = (seg_start + trailing_seq.len()).min(seg.len());
        let seg_seq = &seg[seg_start..seg_end];
        align_tail_inner(seg_seq, trailing_seq)
    } else {
        let read_end = tile.read_start;
        let trailing_seq = bio_utils::revcmp(&seq[..read_end]);
        let seg_end = (seg_start + trailing_seq.len()).min(seg.len());
        let seg_seq = &seg[seg_start..seg_end];
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
    _id: u64,
) -> (Vec<u8>, Vec<Op>, usize, usize) {
    let seg_end = tile.ctg_start;
    if seg_end == 0 {
        return (vec![], vec![], 0, 0);
    }
    let (query, ops) = if chain.is_forward {
        let read_end = tile.read_start;
        let leading_seq = &seq[..read_end];
        let seg_start = seg_end.saturating_sub(leading_seq.len());
        let seg_seq = &seg.sequence.as_ref().unwrap().as_bytes()[seg_start..seg_end];
        align_tip_inner(seg_seq, leading_seq)
    } else {
        let read_start = tile.read_end;
        let leading_seq = bio_utils::revcmp(&seq[read_start..]);
        let seg_start = seg_end.saturating_sub(leading_seq.len());
        let seg_seq = &seg.sequence.as_ref().unwrap().as_bytes()[seg_start..seg_end];
        align_tip_inner(seg_seq, &leading_seq)
    };
    let head_clip = match chain.is_forward {
        true => tile.read_start - query.len(),
        false => seq.len() - tile.read_end - query.len(),
    };
    let len_seg = ops.iter().filter(|&&op| op != Op::Ins).count();
    (query, ops, len_seg, head_clip)
}

use crate::ALN_PARAMETER;
const BAND: usize = 20;
fn align_tip_inner(seg: &[u8], query: &[u8]) -> (Vec<u8>, Vec<Op>) {
    let mode = edlib_sys::AlignMode::Global;
    let task = edlib_sys::AlignTask::Alignment;
    let len = seg.len().min(query.len());
    if len == 0 {
        return (Vec::new(), Vec::new());
    }
    let aln = edlib_sys::align(query, seg, mode, task);
    let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
    let ops = kiley::bialignment::guided::global_guided(seg, query, &ops, BAND, ALN_PARAMETER).1;
    let mut ops =
        kiley::bialignment::guided::global_guided(seg, query, &ops, BAND, ALN_PARAMETER).1;
    ops.reverse();
    let mut start = 0;
    loop {
        match ops.last() {
            Some(&Op::Del) => assert!(ops.pop().is_some()),
            Some(&Op::Ins) => {
                start += 1;
                ops.pop();
            }
            _ => break,
        }
    }
    ops.reverse();
    (query[start..].to_vec(), ops)
}

fn align_tail_inner(seg: &[u8], query: &[u8]) -> (Vec<u8>, Vec<Op>) {
    let mode = edlib_sys::AlignMode::Global;
    let task = edlib_sys::AlignTask::Alignment;
    let len = seg.len().min(query.len());
    if len == 0 {
        return (Vec::new(), Vec::new());
    }
    let aln = edlib_sys::align(query, seg, mode, task);
    let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
    let ops = kiley::bialignment::guided::global_guided(seg, query, &ops, BAND, ALN_PARAMETER).1;
    let (_, mut ops) =
        kiley::bialignment::guided::global_guided(seg, query, &ops, BAND, ALN_PARAMETER);
    let mut poped = 0;
    loop {
        match ops.last() {
            Some(&Op::Del) => assert!(ops.pop().is_some()),
            Some(&Op::Ins) => {
                poped += 1;
                ops.pop();
            }
            _ => break,
        }
    }
    (query[..query.len() - poped].to_vec(), ops)
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

#[derive(Debug, Clone)]
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
    let mut tiles: Vec<Tile> = vec![];
    let (mut read_idx, mut seg_idx) = (chain.query_start_idx, chain.contig_start_idx);
    let mut prev_r_pos = None;
    for &op in chain.ops.iter() {
        match op {
            Op::Match => {
                let read_node = match chain.is_forward {
                    true => &read.nodes[read_idx],
                    false => &read.nodes[read.nodes.len() - read_idx - 1],
                };
                // Determine where the alignment occured.
                let node_start = read_node.position_from_start;
                let seg_node = &encs.tiles()[seg_idx];
                let (start, end) = convert_to_read_coordinate(read_node, seg_node);
                assert!(end <= read_node.seq().len());
                let (read_start, read_end) = (node_start + start, node_start + end);
                let (read_start, read_end) = match (chain.is_forward, prev_r_pos) {
                    (_, None) => (read_start, read_end),
                    (true, Some((_, e))) => (read_start.max(e), read_end.max(e)),
                    (false, Some((s, _))) => (read_start.min(s), read_end.min(s)),
                };
                if read_end - node_start > read_node.seq().len() {
                    error!("{}", chain.is_forward);
                    error!("{chain:?}");
                    for t in tiles.iter() {
                        let len = t.node.seq().len();
                        error!(
                            "{}-{},{len},{}-{}",
                            t.read_start, t.read_end, t.ctg_start, t.ctg_end
                        );
                    }
                    let (ctg_start, ctg_end) = seg_node.contig_range();
                    error!(
                        "{read_start}-{read_end},{start}-{end},{},{ctg_start}-{ctg_end}",
                        read_node.is_forward
                    );
                }
                assert!(read_end - node_start <= read_node.seq().len());
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
    assert!(
        r_end <= tile.node.seq().len(),
        "{},{},{}",
        r_end,
        tile.node.seq().len(),
        tile.node.is_forward,
    );
    let (start, end) = match tile.encoding.unit_range() {
        (true, start, end) => (start, end),
        (false, start, end) => (
            tile.encoding.unit_len() - end,
            tile.encoding.unit_len() - start,
        ),
    };
    assert!(start <= end, "{},{}", start, end);
    let alignment = crate::misc::ops_to_kiley(&tile.node.cigar);
    let r_len = alignment.iter().filter(|&&op| op != Op::Del).count();
    let u_len = alignment.iter().filter(|&&op| op != Op::Ins).count();
    assert_eq!(u_len, tile.encoding.unit_len());
    assert_eq!(r_len, tile.node.seq().len());
    let (mut r_pos, mut u_pos) = (0, 0);
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
    read_id: u64,
    contig: String,
    // 0-based.
    contig_start: usize,
    contig_end: usize,
    query: Vec<u8>,
    query_head_clip: usize,
    query_tail_clip: usize,
    ops: Vec<Op>,
    is_forward: bool,
}

impl Alignment {
    pub fn cigar(&self) -> String {
        let mut cigar = String::new();
        if self.query_head_clip != 0 {
            cigar.push_str(&format!("{}H", self.query_head_clip));
        }
        cigar.push_str(&format!("{}", crate::misc::kiley_op_to_ops(&self.ops)));
        if self.query_tail_clip != 0 {
            cigar.push_str(&format!("{}H", self.query_tail_clip));
        }
        cigar
    }
    pub fn new(
        read_id: u64,
        contig: String,
        (mut contig_start, mut contig_end): (usize, usize),
        (mut query_head_clip, mut query_tail_clip): (usize, usize),
        mut query: Vec<u8>,
        mut ops: Vec<Op>,
        is_forward: bool,
    ) -> Alignment {
        // Remove the last ins/dels
        loop {
            match ops.last() {
                Some(&Op::Del) => contig_end -= ops.pop().is_some() as usize,
                Some(&Op::Ins) => {
                    assert!(ops.pop().is_some());
                    assert!(query.pop().is_some());
                    query_tail_clip += 1;
                }
                _ => break,
            }
        }
        // Remove the leading ins/dels.
        ops.reverse();
        query.reverse();
        loop {
            match ops.last() {
                Some(&Op::Del) => contig_start += ops.pop().is_some() as usize,
                Some(&Op::Ins) => {
                    assert!(ops.pop().is_some());
                    assert!(query.pop().is_some());
                    query_head_clip += 1;
                }
                _ => break,
            }
        }
        ops.reverse();
        query.reverse();
        let q_len = ops.iter().filter(|&&op| op != Op::Del).count();
        let r_len = ops.iter().filter(|&&op| op != Op::Ins).count();
        assert_eq!(q_len, query.len());
        assert_eq!(r_len, contig_end - contig_start);
        Self {
            read_id,
            contig,
            contig_start,
            contig_end,
            query,
            query_head_clip,
            query_tail_clip,
            ops,
            is_forward,
        }
    }
}

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
        let ops = alignment(&query, &refr).0;
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
