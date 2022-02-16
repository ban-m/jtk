use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
const CONS_MIN_LENGTH: usize = 200;
// Directed edge between nodes.
type DEdge = ((u64, u64, bool), (u64, u64, bool));
// (unit,cluster,direction, if it is `from` part)
type DTip = (u64, u64, bool, bool);
#[derive(Debug, Clone, Copy)]
pub struct DenseEncodingConfig {
    pub len: usize,
    // pub min_span_reads: usize,
}
impl DenseEncodingConfig {
    //pub fn new(len: usize, min_span_reads: usize) -> Self {
    pub fn new(len: usize) -> Self {
        Self {
            len,
            // min_span_reads,
        }
    }
}
pub trait DenseEncoding {
    fn dense_encoding_dev(&mut self, config: &DenseEncodingConfig);
    // fn dense_encoding(&mut self, config: &DenseEncodingConfig);
}

impl DenseEncoding for DataSet {
    fn dense_encoding_dev(&mut self, config: &DenseEncodingConfig) {
        let original_assignments = log_original_assignments(self);
        use crate::dirichlet_mixture::{ClusteringConfig, DirichletMixtureCorrection};
        let correction_config = ClusteringConfig::new(5, 10, 5);
        self.correct_clustering(&correction_config);
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
    let inserts = fill_edges_by_new_units_dev(read, seq, edges, tips, &read_type);
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

// fn get_average_unit_length(ds: &DataSet) -> usize {
//     let total: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
//     total / ds.selected_chunks.len()
// }

pub fn fill_edges_by_new_units(
    read: &EncodedRead,
    seq: &[u8],
    edge_encoding_patterns: &HashMap<DEdge, Vec<Unit>>,
    read_type: &definitions::ReadType,
) -> Vec<(usize, Node)> {
    let mut inserts = vec![];
    for (idx, window) in read.nodes.windows(2).enumerate() {
        // What is important is the direction, the label, and the start-end position. Thats all.
        let (start, end) = (
            window[0].position_from_start + window[0].seq().len(),
            window[1].position_from_start,
        );
        let (edge, direction) = edge_and_direction(window);
        let unit_info = match edge_encoding_patterns.get(&edge) {
            Some(units) => units,
            None => continue,
        };
        for node in encode_edge(seq, start, end, direction, unit_info, read_type) {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    inserts
}

pub fn fill_edges_by_new_units_dev(
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
            // if head.unit == 815 {
            //     debug!("DE\tTRACE\tHit\tTailNode");
            // }
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
        for node in encode_edge(seq, start, end, direction, unit_info, read_type) {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    // Tail tip
    if let Some(tail) = read.nodes.last() {
        let idx = read.nodes.len();
        let (start, end) = (tail.position_from_start + tail.seq().len(), seq.len());
        let (start, end) = (start.max(MARGIN) - MARGIN, (end + MARGIN).min(len));
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

// fn split_edges_into_units(
//     edges: &HashMap<DEdge, Vec<u8>>,
//     ave_unit_len: usize,
//     cluster_num: usize,
//     mut max_unit_id: u64,
// ) -> HashMap<DEdge, Vec<Unit>> {
//     let mut units_in_edges = HashMap::new();
//     for (key, consensus) in edges.iter() {
//         let mut units_in_edge = vec![];
//         let mut prev_bp = 0;
//         for break_pos in (1..)
//             .map(|i| i * ave_unit_len)
//             .take_while(|&break_pos| break_pos + ave_unit_len < consensus.len())
//             .chain(std::iter::once(consensus.len()))
//         {
//             max_unit_id += 1;
//             let seq = String::from_utf8_lossy(&consensus[prev_bp..break_pos]).to_string();
//             let unit = Unit::new(max_unit_id, seq, cluster_num);
//             prev_bp = break_pos;
//             units_in_edge.push(unit)
//         }
//         let count = units_in_edge.len();
//         log::debug!("NEWUNIT\t{:?}\t{}\t{}", key, count, consensus.len());
//         units_in_edges.insert(*key, units_in_edge);
//     }
//     units_in_edges
// }

// fn take_consensus_between_nodes<T: std::borrow::Borrow<EncodedRead>>(
//     reads: &[T],
//     nodes: &[(u64, u64)],
//     discard_thr: usize,
//     band_width: usize,
// ) -> HashMap<DEdge, Vec<u8>> {
//     let mut edges: HashMap<DEdge, Vec<Vec<u8>>> = HashMap::new();
//     for read in reads.iter().map(|r| r.borrow()) {
//         for (window, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
//             if !nodes.contains(&(window[0].unit, window[0].cluster))
//                 || !nodes.contains(&(window[1].unit, window[1].cluster))
//             {
//                 continue;
//             }
//             let (edge_entry, direction) = edge_and_direction(window);
//             let mut label = match direction {
//                 true => edge.label().to_vec(),
//                 false => bio_utils::revcmp(edge.label()),
//             };
//             // Register.
//             label.iter_mut().for_each(u8::make_ascii_uppercase);
//             edges.entry(edge_entry).or_insert_with(Vec::new).push(label);
//         }
//     }
//     edges.retain(|_, val| discard_thr < val.len());
//     debug!("EDGE\tDump edges.");
//     for (edge, ls) in edges.iter() {
//         log::debug!("EDGE\t{:?}\t{}", edge, ls.len());
//     }
//     // Create new units.
//     edges
//         .into_par_iter()
//         .filter_map(|(key, mut labels)| {
//             let cov_thr = labels.len() / COV_THR_FACTOR;
//             let mut lengths: Vec<_> = labels.iter().map(|x| x.len()).collect();
//             let (_, median, _) = lengths.select_nth_unstable(labels.len() / 2);
//             let (upper, lower) = (2 * *median, (*median / 2).max(CONS_MIN_LENGTH));
//             let med_idx = labels.len() / 2;
//             labels.select_nth_unstable_by_key(med_idx, |x| x.len());
//             labels.swap(0, med_idx);
//             let labels: Vec<&[u8]> = labels
//                 .iter()
//                 .filter(|ls| lower < ls.len() && ls.len() < upper)
//                 .map(|x| x.as_slice())
//                 .collect();
//             if labels.len() <= cov_thr {
//                 return None;
//             }
//             let rough_contig = kiley::ternary_consensus_by_chunk(&labels, band_width);
//             let rough_contig = kiley::bialignment::polish_until_converge_banded(
//                 &rough_contig,
//                 &labels,
//                 band_width,
//             );
//             match rough_contig.len() {
//                 0..=CONS_MIN_LENGTH => None,
//                 _ => Some((key, rough_contig)),
//             }
//         })
//         .collect()
// }

// Why I did this? The clustering is already squished if the LK gain is not sufficient.
// #[allow(dead_code)]
// fn squish_bad_clustering(ds: &mut DataSet, nodes: &HashSet<u64>, per_read_lk_gain: f64) {
//     let mut new_clustered: HashMap<_, Vec<&mut _>> = nodes.iter().map(|&n| (n, vec![])).collect();
//     for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
//         if let Some(bucket) = new_clustered.get_mut(&node.unit) {
//             bucket.push(&mut node.cluster);
//         }
//     }
//     for unit in ds.selected_chunks.iter_mut() {
//         if let Some(assignments) = new_clustered.get_mut(&unit.id) {
//             let threshold = assignments.len() as f64 * per_read_lk_gain;
//             if unit.score <= threshold {
//                 log::debug!(
//                     "Squishing\t{}\t{}\t{}",
//                     unit.id,
//                     unit.score,
//                     assignments.len()
//                 );
//                 unit.score = 0f64;
//                 assignments.iter_mut().for_each(|x| **x = 0);
//             }
//         }
//     }
// }

// fn enumerate_multitigs(
//     ds: &DataSet,
//     &DenseEncodingConfig { .. }: &DenseEncodingConfig,
// ) -> Vec<(usize, Vec<(u64, u64)>)> {
//     use crate::assemble::*;
//     let min_span_reads = ds.read_type.min_span_reads();
//     let config = AssembleConfig::new(1, 1000, false, true, min_span_reads);
//     let (_, summaries) = assemble(ds, &config);
//     summaries
//         .iter()
//         .filter(|summary| !summary.summary.is_empty())
//         .filter_map(|summary| {
//             let (total_cp, num) = summary
//                 .summary
//                 .iter()
//                 .filter_map(|x| x.copy_number)
//                 .fold((0, 0), |(cp, n), x| (cp + x, n + 1));
//             let copy_number = total_cp / num;
//             let nodes: Vec<_> = summary
//                 .summary
//                 .iter()
//                 .map(|x| (x.unit, x.cluster))
//                 .collect();
//             (1 < copy_number).then(|| (copy_number, nodes))
//         })
//         .collect()
// }

type EdgeAndUnit = HashMap<DEdge, Vec<Unit>>;
fn enumerate_polyploid_edges(ds: &DataSet, _config: &DenseEncodingConfig) -> EdgeAndUnit {
    use crate::assemble::*;
    let min_span_reads = ds.read_type.min_span_reads();
    let config = AssembleConfig::new(1, 1000, false, true, min_span_reads);
    let (_, summaries) = assemble(ds, &config);
    let edges: HashMap<_, _> = summaries
        .iter()
        .filter(|summary| !summary.summary.is_empty())
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
    let mut max_unit_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    debug!("DE\tCand\t{}", consensi_materials.len());
    let mean_chunk_len = {
        let sum: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
        sum / ds.selected_chunks.len()
    };
    let mut newly_defined_unit = HashMap::new();
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
    for (key, consensus) in consensi {
        let copy_num = edges[&key];
        let chunk_num = (consensus.len() as f64 / mean_chunk_len as f64).ceil();
        let chunk_len = (consensus.len() as f64 / chunk_num).ceil() as usize;
        let units: Vec<_> = consensus
            .chunks(chunk_len)
            .map(|seq| {
                let seq = String::from_utf8_lossy(seq).to_string();
                max_unit_id += 1;
                Unit::new(max_unit_id, seq, copy_num)
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
    match kiley::bialignment::polish_until_converge_banded(&rough_contig, &seqs, band_width) {
        seq if seq.len() < CONS_MIN_LENGTH => None,
        seq => Some(seq),
    }
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

// formatting the two nodes so that it is from small unit to unit with larger ID.
// If the direction is reversed, the 2nd argument would be false.
fn edge_and_direction(nodes: &[Node]) -> (DEdge, bool) {
    let (from, to) = match nodes {
        [from, to] => (from, to),
        _ => panic!(),
    };
    match from.unit.cmp(&to.unit) {
        std::cmp::Ordering::Less => {
            let from_elm = (from.unit, from.cluster, from.is_forward);
            let to_elm = (to.unit, to.cluster, !to.is_forward);
            let edge = (from_elm, to_elm);
            (edge, true)
        }
        std::cmp::Ordering::Equal if from.is_forward => {
            let from_elm = (from.unit, from.cluster, from.is_forward);
            let to_elm = (to.unit, to.cluster, !to.is_forward);
            ((from_elm, to_elm), true)
        }
        _ => {
            let from_elm = (to.unit, to.cluster, !to.is_forward);
            let to_elm = (from.unit, from.cluster, from.is_forward);
            ((from_elm, to_elm), false)
        }
    }
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
    // seq is equal to seq[start..end], revcmped if is_forward is false.
    let (start, end, seq, ctg_start, ctg_end) = tune_position(start, end, seq, is_forward, &contig);
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
    let ctg_orig_len = contig.len();
    let contig = &contig[ctg_start..ctg_end];
    let (_, mut ops) =
        kiley::bialignment::global_banded(contig, &seq, 2, -5, -6, -1, read_type.band_width());
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
            let (uid, unitlen) = (unit.id, unit.seq().len());
            let ylen = alignments.iter().filter(|&&x| x != kiley::Op::Del).count();
            let cigar = crate::encode::compress_kiley_ops(&alignments);
            let indel_mism = alignments
                .iter()
                .map(|&op| 1 - 2 * (op == kiley::Op::Match) as i32);
            let max_indel = crate::encode::max_region(indel_mism).max(0) as usize;
            let unitlen = unitlen as f64;
            let gap_thr = ((unitlen * crate::encode::INDEL_FRACTION).round() as usize)
                .max(crate::encode::MIN_INDEL_SIZE);
            let percent_identity = {
                let (aln, mat) = alignments.iter().fold((0, 0), |(aln, mat), &op| match op {
                    kiley::Op::Match => (aln + 1, mat + 1),
                    _ => (aln + 1, mat),
                });
                mat as f64 / aln as f64
            };
            if max_indel < gap_thr && 1f64 - sim_thr < percent_identity {
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
fn tune_position(
    start: usize,
    end: usize,
    seq: &[u8],
    is_forward: bool,
    contig: &[u8],
) -> (usize, usize, Vec<u8>, usize, usize) {
    let mut seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Location;
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
    // TODO:Maybe this is wrong...
    let (seq_start, seq_end) = match is_forward {
        true => (start + seq_start, start + seq_end),
        false => (end - seq_end, end - seq_start),
    };
    // Second, let's truncate the contig.
    let alignment = edlib_sys::edlib_align(&seq, contig, mode, task);
    let (ctg_start, ctg_end) = alignment.locations.unwrap()[0];
    let ctg_end = ctg_end + 1;
    (seq_start, seq_end, seq, ctg_start, ctg_end)
}
