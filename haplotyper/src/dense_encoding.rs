use definitions::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
const CONS_MIN_LENGTH: usize = 100;
const COV_THR_FACTOR: usize = 4;
// Directed edge between nodes.
type DEdge = ((u64, u64, bool), (u64, u64, bool));
#[derive(Debug, Clone, Copy)]
pub struct DenseEncodingConfig {
    pub len: usize,
}
impl DenseEncodingConfig {
    pub fn new(len: usize) -> Self {
        Self { len }
    }
}
pub trait DenseEncoding {
    fn dense_encoding(self, config: &DenseEncodingConfig) -> Self;
}

impl DenseEncoding for DataSet {
    fn dense_encoding(mut self, config: &DenseEncodingConfig) -> Self {
        let ave_unit_len = get_average_unit_length(&self);
        let original_assignments = log_original_assignments(&self);
        use crate::em_correction::ClusteringCorrection;
        let units = self
            .selected_chunks
            .iter()
            .filter_map(|u| (1.0 < u.score).then(|| u.id))
            .collect();
        self = self.correct_clustering_em_on_selected(10, 3, true, &units);
        // The maximum value of the previous unit.
        // If the unit id is greater than this, it is newly added one.
        let multi_tig = enumerate_diplotigs(&mut self, config.len);
        // Nodes to be clustered, or node newly added by filling edge.
        let mut to_clustering_nodes = HashSet::new();
        for (cluster_num, nodes) in multi_tig {
            log::debug!("CONTIG\t{:?}", nodes);
            let mut reads: Vec<_> = self
                .encoded_reads
                .iter_mut()
                .filter(|r| r.nodes.iter().any(|n| nodes.contains(&(n.unit, n.cluster))))
                .collect();
            // Create new nodes between nodes to make this region "dense."
            // (unit,cluster,is_tail).
            let edge_consensus = take_consensus_between_nodes(&reads, &nodes, 2);
            let max_unit_id: u64 = self.selected_chunks.iter().map(|x| x.id).max().unwrap();
            let edge_encoding_patterns =
                split_edges_into_units(&edge_consensus, ave_unit_len, max_unit_id);
            for (key, new_units) in edge_encoding_patterns.iter() {
                let contig = &edge_consensus[key];
                let mut prev_bp = 0;
                log::debug!("NEWUNIT\t{:?}\t{:?}\t{}", key, new_units, contig.len());
                for &(break_point, id) in new_units {
                    let seq = String::from_utf8_lossy(&contig[prev_bp..break_point]).to_string();
                    let unit = Unit::new(id, seq, cluster_num);
                    self.selected_chunks.push(unit);
                    prev_bp = break_point;
                }
            }
            let read_seq: HashMap<u64, &[u8]> =
                self.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
            // Encoding.
            for read in reads.iter_mut() {
                let seq = &read_seq[&read.id];
                let inserts =
                    fill_edges_by_new_units(&read, seq, &edge_encoding_patterns, &edge_consensus);
                // Encode.
                for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
                    read.nodes.insert(idx + accum_inserts, node);
                }
                // Formatting.
                re_encode_read(read, seq);
            }
            let filled_reads: HashSet<_> = reads.iter().map(|r| r.id).collect();
            crate::encode::deletion_fill::correct_unit_deletion_selected(&mut self, &filled_reads);
            to_clustering_nodes.extend(nodes.iter().map(|x| x.0));
            edge_encoding_patterns.values().for_each(|new_units| {
                to_clustering_nodes.extend(new_units.iter().map(|(_, id)| *id))
            });
        }
        for read in self.encoded_reads.iter_mut() {
            let orig = &original_assignments[&read.id];
            recover_original_assignments(read, &orig);
        }
        // Local clustering.
        crate::local_clustering::local_clustering_selected(&mut self, &to_clustering_nodes);
        squish_bad_clustering(&mut self, &to_clustering_nodes, 1f64);
        self
    }
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
        while let Some(node) = read.next() {
            if node.unit == unit {
                node.cluster = cluster;
                break;
            }
        }
    }
}

fn get_average_unit_length(ds: &DataSet) -> usize {
    let total: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
    total / ds.selected_chunks.len()
}

fn fill_edges_by_new_units(
    read: &EncodedRead,
    seq: &[u8],
    edge_encoding_patterns: &HashMap<DEdge, Vec<(usize, u64)>>,
    edge_consensus: &HashMap<DEdge, Vec<u8>>,
) -> Vec<(usize, Node)> {
    let mut inserts = vec![];
    for (idx, window) in read.nodes.windows(2).enumerate() {
        // What is important is the direction, the label, and the start-end position. Thats all.
        let (start, end) = (
            window[0].position_from_start + window[0].seq().len(),
            window[1].position_from_start,
        );
        let (edge, direction) = edge_and_direction(window);
        let contig = match edge_consensus.get(&edge) {
            Some(contig) => contig,
            None => continue,
        };
        let unit_info = &edge_encoding_patterns[&edge];
        log::debug!("TRY\t{:?}\t{}...{}bp({})", edge, start, end, end - start);
        for node in encode_edge(seq, start, end, direction, contig, unit_info) {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    inserts
}

fn split_edges_into_units(
    edges: &HashMap<DEdge, Vec<u8>>,
    ave_unit_len: usize,
    mut max_unit_id: u64,
) -> HashMap<DEdge, Vec<(usize, u64)>> {
    edges
        .iter()
        .map(|(&key, consensus)| {
            let break_positions: Vec<_> = (1..)
                .map(|i| i * ave_unit_len)
                .take_while(|&break_pos| break_pos + ave_unit_len < consensus.len())
                .chain(std::iter::once(consensus.len()))
                .map(|break_pos| {
                    max_unit_id += 1;
                    (break_pos, max_unit_id)
                })
                .collect();
            (key, break_positions)
        })
        .collect()
}

fn take_consensus_between_nodes<T: std::borrow::Borrow<EncodedRead>>(
    reads: &[T],
    nodes: &Vec<(u64, u64)>,
    discard_thr: usize,
) -> HashMap<DEdge, Vec<u8>> {
    let mut edges: HashMap<DEdge, Vec<Vec<u8>>> = HashMap::new();
    for read in reads.iter().map(|r| r.borrow()) {
        for (window, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
            if !nodes.contains(&(window[0].unit, window[0].cluster))
                || !nodes.contains(&(window[1].unit, window[1].cluster))
            {
                continue;
            }
            let (edge_entry, direction) = edge_and_direction(window);
            let mut label = match direction {
                true => edge.label().to_vec(),
                false => bio_utils::revcmp(edge.label()),
            };
            // Register.
            label.iter_mut().for_each(u8::make_ascii_uppercase);
            edges.entry(edge_entry).or_insert(vec![]).push(label);
        }
    }
    debug!("EDGE\tDump edges. Contains empty edges.");
    for (edge, ls) in edges.iter() {
        log::debug!("EDGE\t{:?}\t{}", edge, ls.len());
    }
    edges.retain(|_, val| discard_thr < val.len());
    // Create new units.
    edges
        .into_par_iter()
        .filter_map(|(key, labels)| {
            let cov_thr = labels.len() / COV_THR_FACTOR;
            let labels: Vec<&[u8]> = labels
                .iter()
                .filter(|ls| ls.len() > CONS_MIN_LENGTH)
                .map(|x| x.as_slice())
                .collect();
            if labels.len() < cov_thr {
                return None;
            }
            let rough_contig = kiley::ternary_consensus_by_chunk(&labels, 100);
            match rough_contig.len() {
                0..=CONS_MIN_LENGTH => None,
                _ => Some((key, rough_contig)),
            }
        })
        .collect()
}

fn squish_bad_clustering(ds: &mut DataSet, nodes: &HashSet<u64>, per_read_lk_gain: f64) {
    let mut new_clustered: HashMap<_, Vec<&mut _>> = nodes.iter().map(|&n| (n, vec![])).collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some(bucket) = new_clustered.get_mut(&node.unit) {
            bucket.push(&mut node.cluster);
        }
    }
    for unit in ds.selected_chunks.iter_mut() {
        if let Some(assignments) = new_clustered.get_mut(&unit.id) {
            // At least 1 LK for each element(CLRmode)
            let threshold = assignments.len() as f64 * per_read_lk_gain;
            if unit.score <= threshold {
                log::debug!(
                    "Squishing\t{}\t{}\t{}",
                    unit.id,
                    unit.score,
                    assignments.len()
                );
                unit.score = 0f64;
                assignments.iter_mut().for_each(|x| **x = 0);
            }
        }
    }
}

// Squish short contig, return the generated diplotigs.
fn enumerate_diplotigs(ds: &mut DataSet, len: usize) -> Vec<(usize, Vec<(u64, u64)>)> {
    use crate::assemble::*;
    let config = AssembleConfig::new(1, 1000, false, true);
    ds.squish_small_contig(&config, len);
    let (_, summaries) = assemble(ds, 0, &config);
    let mut multi_tig: Vec<_> = summaries
        .iter()
        .filter_map(|summary| {
            let (total_cp, num) = summary
                .summary
                .iter()
                .filter_map(|x| x.copy_number)
                .fold((0, 0), |(cp, n), x| (cp + x, n + 1));
            let copy_number = total_cp / num;
            let nodes: Vec<_> = summary
                .summary
                .iter()
                .map(|x| (x.unit, x.cluster))
                .collect();
            (1 < copy_number).then(|| (copy_number, nodes))
        })
        .collect();
    multi_tig.sort_by_key(|x| x.1.len());
    multi_tig
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

fn encode_edge(
    seq: &[u8],
    start: usize,
    end: usize,
    is_forward: bool,
    contig: &[u8],
    unit_info: &[(usize, u64)],
) -> Vec<definitions::Node> {
    let mut seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    use log::*;
    debug!("Encoding\t{}\t{}\t{:?}", seq.len(), contig.len(), unit_info);
    let band = contig.len() / 20;
    // TODO: These parameters shoule be changed depending on the sequencing tech.
    let (_, ops) = kiley::bialignment::global_banded(contig, &seq, 2, -2, -4, -2, band);
    // Split the alignment into encoded nodes.
    // Current position on the contg and the query.
    let (mut xpos, mut ypos) = (0, 0);
    // Current edit distance inside the focal unit.
    let mut edit_dist = 0;
    // Current edit operations inside the focal unit.
    let mut alignments = vec![];
    let mut target_idx = 0;
    // Encoded nodes.
    let mut nodes = vec![];
    for op in ops {
        match op {
            kiley::bialignment::Op::Mat => {
                edit_dist += (contig[xpos] != seq[ypos]) as u32;
                xpos += 1;
                ypos += 1;
            }
            kiley::bialignment::Op::Del => {
                edit_dist += 1;
                xpos += 1;
            }
            kiley::bialignment::Op::Ins => {
                edit_dist += 1;
                ypos += 1;
            }
        }
        alignments.push(op);
        if xpos == unit_info[target_idx].0 {
            // Reached the boundary.
            let (break_point, uid) = unit_info[target_idx];
            let unitlen = match target_idx {
                0 => break_point,
                _ => break_point - unit_info[target_idx - 1].0,
            };
            let ylen = alignments
                .iter()
                .filter(|&&x| x != kiley::bialignment::Op::Del)
                .count();
            let cigar = crate::encode::compress_kiley_ops(&alignments);
            let has_large_indel = cigar.iter().any(|op| match *op {
                Op::Ins(l) | Op::Del(l) => l > crate::encode::INDEL_THRESHOLD,
                _ => false,
            });
            use crate::encode::deletion_fill::ALIGN_LIMIT;
            let dist_thr = (unitlen as f64 * ALIGN_LIMIT).floor() as u32;
            if !has_large_indel && edit_dist < dist_thr {
                nodes.push(Node {
                    position_from_start: start + ypos - ylen,
                    unit: uid,
                    cluster: 0,
                    is_forward,
                    seq: String::from_utf8_lossy(&seq[ypos - ylen..ypos]).to_string(),
                    cigar,
                })
            }
            // Refresh.
            target_idx += 1;
            alignments.clear();
            edit_dist = 0;
        }
        if target_idx == unit_info.len() {
            // This is needed, as sometimes only insertions would be remain.
            break;
        }
    }
    if !is_forward {
        nodes.reverse()
    }
    nodes
}
