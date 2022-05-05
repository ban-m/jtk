use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::ALN_PARAMETER;
const CONS_MIN_LENGTH: usize = 200;
// Directed edge between nodes.
type DEdge = ((u64, u64, bool), (u64, u64, bool));
// (unit,cluster,direction, if it is `from` part)
type DTip = (u64, u64, bool, bool);
#[derive(Debug, Clone)]
pub struct DenseEncodingConfig {
    pub len: usize,
    pub file: Option<String>,
}
impl DenseEncodingConfig {
    pub fn new(len: usize, file: Option<&str>) -> Self {
        Self {
            len,
            file: file.map(|f| f.to_string()),
        }
    }
}
pub trait DenseEncoding {
    fn dense_encoding_dev(&mut self, config: &DenseEncodingConfig);
}

impl DenseEncoding for DataSet {
    fn dense_encoding_dev(&mut self, config: &DenseEncodingConfig) {
        let original_assignments = log_original_assignments(self);
        // let msr = self.read_type.min_span_reads();
        // let min_lk = self.read_type.min_llr_value();
        // use crate::assemble::*;
        // let asm_config = AssembleConfig::new(1, 1000, false, true, msr, min_lk);
        // self.squish_small_contig(&asm_config, 3);
        let new_units = encode_polyploid_edges(self, config);
        for read in self.encoded_reads.iter_mut() {
            let orig = &original_assignments[&read.id];
            recover_original_assignments(read, orig);
        }
        crate::local_clustering::local_clustering_selected(self, &new_units);
    }
}

type TipAndUnit<'a> = HashMap<DTip, Vec<&'a [Unit]>>;
fn encode_polyploid_edges(ds: &mut DataSet, config: &DenseEncodingConfig) -> HashSet<u64> {
    let edge_units = enumerate_polyploid_edges(ds, config);
    let tip_units = {
        let mut tip_units: TipAndUnit = HashMap::new();
        for (&(from, to), units) in edge_units.iter() {
            let from = (from.0, from.1, from.2, true);
            let to = (to.0, to.1, to.2, false);
            tip_units.entry(from).or_default().push(units.as_slice());
            tip_units.entry(to).or_default().push(units.as_slice());
        }
        tip_units
    };
    let rawseq: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let readtype = ds.read_type;
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let rawseq = &rawseq[&read.id];
        edge_encode(read, rawseq, &edge_units, &tip_units, readtype, config)
    });
    let unit_ids: HashSet<_> = edge_units
        .values()
        .flat_map(|x| x.iter().map(|x| x.id))
        .collect();
    {
        let lens: HashMap<_, _> = edge_units
            .values()
            .flat_map(|x| x.iter().map(|x| (x.id, x.seq().len())))
            .collect();
        let mut counts: HashMap<_, u32> = lens.keys().map(|&id| (id, 0)).collect();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *counts.entry(node.unit).or_default() += 1;
        }
        for id in unit_ids.iter() {
            let (count, len) = (counts[id], lens[id]);
            debug!("DE\tCount\t{}\t{}\t{}", id, count, len);
        }
    }
    ds.selected_chunks
        .extend(edge_units.into_values().flatten());
    unit_ids
}

// TODO: we do not need raw seq `seq`. But I'm afraid read.rawseq() might be broken.
fn edge_encode(
    read: &mut EncodedRead,
    seq: &[u8],
    edges: &EdgeAndUnit,
    tips: &TipAndUnit,
    read_type: definitions::ReadType,
    _config: &DenseEncodingConfig,
) {
    let inserts = fill_edges_by_new_units(read, seq, edges, tips, &read_type);
    for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
        match idx + accum_inserts {
            pos if pos < read.nodes.len() => read.nodes.insert(idx + accum_inserts, node),
            _ => read.nodes.push(node),
        }
    }
    re_encode_read(read, seq);
}

fn re_encode_read(read: &mut EncodedRead, seq: &[u8]) {
    if !read.nodes.is_empty() {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use crate::encode::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| n.unit);
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
}

fn log_original_assignments(ds: &DataSet) -> HashMap<u64, Vec<(u64, u64)>> {
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
fn recover_original_assignments(read: &mut EncodedRead, log: &[(u64, u64)]) {
    let mut read = read.nodes.iter_mut();
    for &(unit, cluster) in log {
        for node in &mut read {
            // while let Some(node) = read.next() {
            if node.unit == unit {
                node.cluster = cluster;
                break;
            }
        }
    }
}

pub fn fill_edges_by_new_units(
    read: &EncodedRead,
    seq: &[u8],
    edges: &EdgeAndUnit,
    tips: &TipAndUnit,
    read_type: &definitions::ReadType,
) -> Vec<(usize, Node)> {
    const MARGIN: usize = 25;
    let len = seq.len();
    let mut inserts = vec![];
    // Head tip.
    if let Some(head) = read.nodes.first() {
        let (start, end) = (0, (head.position_from_start + MARGIN).min(len));
        let key = (head.unit, head.cluster, head.is_forward, false);
        if let Some(units) = tips.get(&key) {
            // --Tip--|Node[0]>|------
            // --Unit-|ToNode|-----
            for unit_info in units {
                let nodes = encode_edge(seq, start, end, true, unit_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (0, x)));
            }
        }
        // Here is a bug
        let key = (head.unit, head.cluster, !head.is_forward, true);
        if let Some(units) = tips.get(&key) {
            // |<Node[0]|-Tip--
            // |FromNode|-Unit-
            for unit_info in units {
                let nodes = encode_edge(seq, start, end, false, unit_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (0, x)));
            }
        }
    }
    for (idx, window) in read.nodes.windows(2).enumerate() {
        let start = window[0].position_from_start + window[0].seq().len();
        let end = window[1].position_from_start;
        let (start, end) = (start.max(MARGIN) - MARGIN, (end + MARGIN).min(len));
        let forward = get_forward_d_edge_from_window(window);
        let reverse = get_reverse_d_edge_from_window(window);
        let (unit_info, direction) = if edges.contains_key(&forward) {
            (&edges[&forward], true)
        } else if edges.contains_key(&reverse) {
            (&edges[&reverse], false)
        } else {
            continue;
        };
        if end <= start {
            let start = window[0].position_from_start;
            let end = window[1].position_from_start;
            let seqlen = window[0].seq().len();
            panic!("{}\t{}\t{}\t{}\t{}", read.id, start, end, seqlen, len);
        }
        for node in encode_edge(seq, start, end, direction, unit_info, read_type) {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    // Tail tip
    if let Some(tail) = read.nodes.last() {
        let idx = read.nodes.len();
        let (start, end) = ((tail.position_from_start + tail.seq().len()).min(len), len);
        let (start, end) = (start.saturating_sub(MARGIN), (end + MARGIN).min(len));
        let key = (tail.unit, tail.cluster, tail.is_forward, true);
        if let Some(units) = tips.get(&key) {
            // | Last>  |-Tip--
            // |FromNode|-Unit-
            for unit_info in units {
                let nodes = encode_edge(seq, start, end, true, unit_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (idx + 1, x)));
            }
        }
        let key = (tail.unit, tail.cluster, !tail.is_forward, false);
        if let Some(units) = tips.get(&key) {
            // --Tip-|<Node[0]|
            // -Unit-|ToNode|
            for unit_info in units {
                let nodes = encode_edge(seq, start, end, false, unit_info, read_type);
                inserts.extend(nodes.into_iter().map(|x| (idx + 1, x)));
            }
        }
    }
    inserts
}

// fn weak_resolve(read_type: definitions::ReadType) -> (usize, f64) {
//     match read_type {
//         ReadType::CCS => (3, 2f64),
//         ReadType::CLR => (4, 5f64),
//         ReadType::ONT => (4, 4f64),
//         ReadType::None => (4, 4f64),
//     }
// }

type EdgeAndUnit = HashMap<DEdge, Vec<Unit>>;
fn enumerate_polyploid_edges(ds: &DataSet, de_config: &DenseEncodingConfig) -> EdgeAndUnit {
    use crate::assemble::*;
    let msr = ds.read_type.min_span_reads();
    let min_lk = ds.read_type.min_llr_value();
    // let (min_span_reads, lk_ratio) = weak_resolve(ds.read_type);
    let config = AssembleConfig::new(1, 1000, false, true, msr, min_lk);
    let (records, summaries) = assemble(ds, &config);
    if let Some(file) = de_config.file.as_ref() {
        let header = gfa::Content::Header(gfa::Header::default());
        let header = gfa::Record::from_contents(header, vec![]);
        let records = std::iter::once(header).chain(records.clone()).collect();
        let gfa = gfa::GFA::from_records(records);
        if let Ok(mut wtr) = std::fs::File::create(file).map(std::io::BufWriter::new) {
            use std::io::Write;
            if let Err(why) = writeln!(&mut wtr, "{}", gfa) {
                eprintln!("{:?}", why);
            }
        }
    }
    let edges: HashMap<_, _> = summaries
        .iter()
        .filter(|summary| !summary.summary.is_empty())
        .filter(|summary| {
            let mut has_edges = [false, false];
            for rec in records.iter() {
                if let &gfa::Content::Edge(ref edge) = &rec.content {
                    if edge.sid1.id == summary.id {
                        has_edges[edge.end1.is_last as usize] |= true;
                    }
                    if edge.sid2.id == summary.id {
                        has_edges[edge.end2.is_last as usize] |= true;
                    }
                }
            }
            has_edges[0] && has_edges[1]
        })
        .flat_map(|summary| {
            let (total_cp, num) = summary
                .summary
                .iter()
                .filter_map(|x| x.copy_number)
                .fold((0, 0), |(cp, n), x| (cp + x, n + 1));
            match total_cp / num {
                0 | 1 => Vec::new(),
                copy_number => summary
                    .summary
                    .windows(2)
                    .map(|w| {
                        let from = (w[0].unit, w[0].cluster, w[0].strand);
                        let to = (w[1].unit, w[1].cluster, w[1].strand);
                        ((from, to), copy_number)
                    })
                    .collect(),
            }
        })
        .collect();
    debug!("DE\t{}\tEDGES", edges.len());
    let mut consensi_materials: HashMap<_, Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        assert_eq!(read.nodes.len(), read.edges.len() + 1);
        for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
            assert_eq!((edge.from, edge.to), (w[0].unit, w[1].unit));
            let forward = get_forward_d_edge_from_window(w);
            let reverse = get_forward_d_edge_from_window(w);
            if edges.contains_key(&forward) {
                let label = edge.label().to_vec();
                consensi_materials.entry(forward).or_default().push(label);
            } else if edges.contains_key(&reverse) {
                let label = bio_utils::revcmp(edge.label());
                consensi_materials.entry(reverse).or_default().push(label);
            }
        }
    }
    let cov_thr = ds.coverage.unwrap().ceil() as usize / 4;
    debug!("DE\tCand\t{}", consensi_materials.len());
    let mean_chunk_len = {
        let sum: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
        sum / ds.selected_chunks.len()
    };
    let mut consensi: Vec<_> = consensi_materials
        .into_par_iter()
        .filter_map(|(key, mut seqs)| {
            seqs.iter_mut()
                .for_each(|xs| xs.iter_mut().for_each(u8::make_ascii_uppercase));
            let mean = seqs.iter().map(|x| x.len()).sum::<usize>() / seqs.len();
            let len = seqs.len();
            let cons = consensus(seqs, cov_thr).map(|s| (key, s));
            if cons.is_none() {
                debug!("DE\tEDGE\t{}\t{}\tNG", len, mean);
            } else {
                debug!("DE\tEDGE\t{}\t{}\tOK", len, mean);
            }
            cons
        })
        .collect();
    consensi.sort_unstable_by_key(|x| x.0);
    let mut newly_defined_unit = HashMap::new();
    let mut max_unit_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    for (key, consensus) in consensi {
        let copy_num = edges[&key];
        let chunk_num = (consensus.len() as f64 / mean_chunk_len as f64).ceil();
        let chunk_len = (consensus.len() as f64 / chunk_num).ceil() as usize;
        let units: Vec<_> = consensus
            .chunks(chunk_len)
            .map(|seq| {
                max_unit_id += 1;
                Unit::new(max_unit_id, seq.to_vec(), copy_num)
            })
            .collect();
        let edge = format!("({},{})-({},{})", key.0 .0, key.0 .2, key.1 .0, key.1 .2);
        debug!(
            "DE\tNewUnit\t{}\t{}\t{}\t{}",
            consensus.len(),
            max_unit_id - units.len() as u64,
            max_unit_id,
            edge
        );
        newly_defined_unit.insert(key, units);
    }
    debug!("DE\tDefined\t{}", newly_defined_unit.len());
    newly_defined_unit
}

fn consensus(mut seqs: Vec<Vec<u8>>, cov_thr: usize) -> Option<Vec<u8>> {
    let pos = seqs.len() / 2;
    let median = seqs.select_nth_unstable_by_key(pos, |x| x.len()).1.len();
    let (upper, lower) = (2 * median, median.max(CONS_MIN_LENGTH) / 2);
    let idx = seqs.iter().position(|x| x.len() == median).unwrap();
    seqs.swap(0, idx);
    seqs.retain(|x| (lower..upper).contains(&x.len()));
    if seqs.len() <= cov_thr {
        return None;
    }
    let mean_len = seqs.iter().map(|x| x.len()).sum::<usize>() / seqs.len();
    let band_width = (mean_len / 20).max(50);
    let rough_contig = kiley::ternary_consensus_by_chunk(&seqs, band_width);
    let cons = kiley::bialignment::guided::polish_until_converge(&rough_contig, &seqs, band_width);
    (CONS_MIN_LENGTH < cons.len()).then(|| cons)
}

// w: windows of nodes with 2 length.
fn get_forward_d_edge_from_window(w: &[Node]) -> DEdge {
    let from = (w[0].unit, w[0].cluster, w[0].is_forward);
    let to = (w[1].unit, w[1].cluster, w[1].is_forward);
    (from, to)
}

// w: windows of nodes with 2 length.
fn get_reverse_d_edge_from_window(w: &[Node]) -> DEdge {
    let from = (w[1].unit, w[1].cluster, !w[1].is_forward);
    let to = (w[0].unit, w[0].cluster, !w[0].is_forward);
    (from, to)
}

// Note that the seq[start..end] can be much longer than the contig itself....
fn encode_edge(
    seq: &[u8],
    start: usize,
    end: usize,
    is_forward: bool,
    units: &[Unit],
    read_type: &definitions::ReadType,
) -> Vec<definitions::Node> {
    let contig: Vec<_> = units.iter().map(|x| x.seq()).fold(Vec::new(), |mut x, y| {
        x.extend(y);
        x
    });
    let ctg_orig_len = contig.len();
    // seq is equal to seq[start..end], revcmped if is_forward is false.
    let band = read_type.band_width(contig.len());
    let ((start, end, seq), (ctg_start, ctg_end, _), mut ops) =
        tune_position(start, end, seq, is_forward, &contig, band);
    let (break_points, _): (Vec<_>, _) =
        units
            .iter()
            .map(|x| x.seq().len())
            .fold((Vec::new(), 0), |(mut acc, x), y| {
                acc.push(x + y);
                (acc, x + y)
            });
    // Current position on the contg
    let mut xpos = ctg_start;
    // Nearest break;
    let mut target_idx = match break_points.iter().position(|&x| ctg_start < x) {
        Some(i) => i,
        None => return Vec::new(),
    };
    // Current edit operations inside the focal unit.
    let mut alignments = match target_idx {
        0 => vec![kiley::Op::Del; ctg_start],
        i => vec![kiley::Op::Del; ctg_start - break_points[i - 1]],
    };
    // Push deletion operations up to the last base.
    ops.extend(std::iter::repeat(kiley::Op::Del).take(ctg_orig_len - ctg_end));
    // Split the alignment into encoded nodes.
    // Current position of the query.
    let mut ypos = {
        let mut head_ins = 0;
        ops.reverse();
        while let Some(&kiley::Op::Ins) = ops.last() {
            ops.pop().unwrap();
            head_ins += 1;
        }
        ops.reverse();
        head_ins
    };
    // Encoded nodes.
    let mut nodes = vec![];
    let sim_thr = read_type.sim_thr();
    for op in ops {
        match op {
            kiley::Op::Match | kiley::Op::Mismatch => {
                xpos += 1;
                ypos += 1;
            }
            kiley::Op::Del => {
                xpos += 1;
            }
            kiley::Op::Ins => {
                ypos += 1;
            }
        }
        alignments.push(op);
        if xpos == break_points[target_idx] {
            // Reached the boundary.
            let unit = &units[target_idx];
            let (uid, _unitlen) = (unit.id, unit.seq().len());
            let ylen = alignments.iter().filter(|&&x| x != kiley::Op::Del).count();
            let cigar = crate::encode::compress_kiley_ops(&alignments);
            let percent_identity = {
                let (aln, mat) = alignments.iter().fold((0, 0), |(aln, mat), &op| match op {
                    kiley::Op::Match => (aln + 1, mat + 1),
                    _ => (aln + 1, mat),
                });
                mat as f64 / aln as f64
            };
            //if max_indel < gap_thr && 1f64 - sim_thr < percent_identity {
            if 1f64 - sim_thr < percent_identity {
                let position_from_start = match is_forward {
                    true => start + ypos - ylen,
                    false => end - ypos,
                };
                let seq = &seq[ypos - ylen..ypos];
                let cl = unit.cluster_num;
                let node = Node::new(uid, is_forward, seq, cigar, position_from_start, cl);
                nodes.push(node);
            }
            // Refresh.
            target_idx += 1;
            alignments.clear();
        }
        if target_idx == units.len() {
            // This is needed, as sometimes only insertions would be remain.
            break;
        }
    }
    if !is_forward {
        nodes.reverse()
    }
    nodes
}

// Note that the seq and the contig is sometimes like:
// seq:----------------------
// ctg:    --------------
// And sometimes like:
// seq:----------(Tip)
// ctg:   -------------------(Long edge)
// This.
// So, we need two-round infix alignment so that, in first alignment, the sequence would be truncted,
// and the second round the contig would be tructed.
type TunedPosition<'a> = (
    (usize, usize, Vec<u8>),
    (usize, usize, &'a [u8]),
    Vec<kiley::Op>,
);
fn tune_position<'a>(
    start: usize,
    end: usize,
    seq: &[u8],
    is_forward: bool,
    contig: &'a [u8],
    band: usize,
) -> TunedPosition<'a> {
    assert!(start < end);
    let mut seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    // First, let's truncate the `seq`
    let alignment = edlib_sys::edlib_align(contig, &seq, mode, task);
    let (seq_start, seq_end) = alignment.locations.unwrap()[0];
    let seq_end = seq_end + 1;
    // Modefy seq into truncate version.
    let pop_num = seq.len() - seq_end;
    (0..pop_num).map(|_| seq.pop().unwrap()).count();
    seq.reverse();
    (0..seq_start).map(|_| seq.pop().unwrap()).count();
    seq.reverse();
    // In the original coordinate.
    let (seq_start, seq_end) = match is_forward {
        true => (start + seq_start, start + seq_end),
        false => (end - seq_end, end - seq_start),
    };
    // Second, let's truncate the contig.
    let alignment = edlib_sys::edlib_align(&seq, contig, mode, task);
    let (ctg_start, ctg_end) = alignment.locations.unwrap()[0];
    let ctg_end = ctg_end + 1;
    let contig = &contig[ctg_start..ctg_end];
    let alignment = alignment.operations.unwrap();
    let edlib_to_op = {
        use kiley::Op::*;
        [Match, Ins, Del, Mismatch]
    };
    let ops: Vec<_> = alignment.iter().map(|&x| edlib_to_op[x as usize]).collect();
    let (_, ops) =
        kiley::bialignment::guided::global_guided(contig, &seq, &ops, band, ALN_PARAMETER);
    ((seq_start, seq_end, seq), (ctg_start, ctg_end, contig), ops)
}
