use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
const CONS_MIN_LENGTH: usize = 200;
const COV_THR_FACTOR: usize = 4;
// Directed edge between nodes.
type DEdge = ((u64, u64, bool), (u64, u64, bool));
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
    fn dense_encoding(&mut self, config: &DenseEncodingConfig);
}

impl DenseEncoding for DataSet {
    fn dense_encoding(&mut self, config: &DenseEncodingConfig) {
        let (band_width, sim_thr) = (self.read_type.band_width(), self.read_type.sim_thr());
        let read_type = self.read_type;
        let ave_unit_len = get_average_unit_length(self);
        let original_assignments = log_original_assignments(self);
        let units: HashSet<_> = self
            .selected_chunks
            .iter()
            .filter_map(|u| (1.0 < u.score).then(|| u.id))
            .collect();
        debug!("DE\t{}\tEMCorrection", units.len());
        {
            use crate::dirichlet_mixture::{ClusteringConfig, DirichletMixtureCorrection};
            let config = ClusteringConfig::new(5, 10, 3);
            self.correct_clustering_on_selected(&config, &units);
        }
        debug!("DE\t{:?}\tEnumDiplotig", config);
        // The maximum value of the previous unit.
        // If the unit id is greater than this, it is newly added one.
        let multi_tig = enumerate_multitigs(self, config);
        // Nodes to be clustered, or node newly added by filling edge.
        let mut to_clustering_nodes = HashSet::new();
        let cov = self
            .coverage
            .unwrap_or_else(|| panic!("do not have coverage"));
        let read_seq: HashMap<u64, &[u8]> =
            self.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
        for (cluster_num, nodes) in multi_tig {
            log::debug!("CONTIG\t{:?}", nodes);
            let mut reads: Vec<_> = self
                .encoded_reads
                .iter_mut()
                .filter(|r| r.nodes.iter().any(|n| nodes.contains(&(n.unit, n.cluster))))
                .collect();
            // Create new nodes between nodes to make this region "dense."
            // (unit,cluster,is_tail).
            let discard_thr = cov.floor() as usize / 2;
            let edge_consensus =
                take_consensus_between_nodes(&reads, &nodes, discard_thr, band_width);
            let max_unit_id: u64 = self.selected_chunks.iter().map(|x| x.id).max().unwrap();
            let edge_encoding_patterns =
                split_edges_into_units(&edge_consensus, ave_unit_len, cluster_num, max_unit_id);
            // Encoding.
            for read in reads.iter_mut() {
                let seq = &read_seq[&read.id];
                let inserts =
                    fill_edges_by_new_units(read, seq, &edge_encoding_patterns, &read_type);
                // Encode.
                for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
                    read.nodes.insert(idx + accum_inserts, node);
                }
                // Formatting.
                re_encode_read(read, seq);
            }
            to_clustering_nodes.extend(nodes.iter().map(|x| x.0));
            to_clustering_nodes.extend(edge_encoding_patterns.values().flatten().map(|x| x.id));
            self.selected_chunks.extend(
                edge_encoding_patterns
                    .into_iter()
                    .flat_map(|(_, units)| units),
            );
        }
        use crate::encode::deletion_fill::correct_unit_deletion;
        let filled_units = correct_unit_deletion(self, sim_thr);
        to_clustering_nodes.extend(filled_units);
        for read in self.encoded_reads.iter_mut() {
            let orig = &original_assignments[&read.id];
            recover_original_assignments(read, orig);
        }
        // Local clustering.
        debug!("LOCAL\tNEW\t{}", to_clustering_nodes.len());
        crate::local_clustering::local_clustering_selected(self, &to_clustering_nodes);
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
        for node in &mut read {
            // while let Some(node) = read.next() {
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
        // for node in encode_edge(seq, start, end, direction, contig, unit_info) {
        for node in encode_edge(seq, start, end, direction, unit_info, read_type) {
            // idx=0 -> Insert at the first edge. So, the index should be 1.
            inserts.push((idx + 1, node));
        }
    }
    inserts
}

fn split_edges_into_units(
    edges: &HashMap<DEdge, Vec<u8>>,
    ave_unit_len: usize,
    cluster_num: usize,
    mut max_unit_id: u64,
) -> HashMap<DEdge, Vec<Unit>> {
    let mut units_in_edges = HashMap::new();
    for (key, consensus) in edges.iter() {
        let mut units_in_edge = vec![];
        let mut prev_bp = 0;
        for break_pos in (1..)
            .map(|i| i * ave_unit_len)
            .take_while(|&break_pos| break_pos + ave_unit_len < consensus.len())
            .chain(std::iter::once(consensus.len()))
        {
            max_unit_id += 1;
            let seq = String::from_utf8_lossy(&consensus[prev_bp..break_pos]).to_string();
            let unit = Unit::new(max_unit_id, seq, cluster_num);
            prev_bp = break_pos;
            units_in_edge.push(unit)
        }
        let count = units_in_edge.len();
        log::debug!("NEWUNIT\t{:?}\t{}\t{}", key, count, consensus.len());
        units_in_edges.insert(*key, units_in_edge);
    }
    units_in_edges
}

fn take_consensus_between_nodes<T: std::borrow::Borrow<EncodedRead>>(
    reads: &[T],
    nodes: &[(u64, u64)],
    discard_thr: usize,
    band_width: usize,
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
            edges.entry(edge_entry).or_insert_with(Vec::new).push(label);
        }
    }
    edges.retain(|_, val| discard_thr < val.len());
    debug!("EDGE\tDump edges.");
    for (edge, ls) in edges.iter() {
        log::debug!("EDGE\t{:?}\t{}", edge, ls.len());
    }
    // Create new units.
    edges
        .into_par_iter()
        .filter_map(|(key, mut labels)| {
            let cov_thr = labels.len() / COV_THR_FACTOR;
            let mut lengths: Vec<_> = labels.iter().map(|x| x.len()).collect();
            let (_, median, _) = lengths.select_nth_unstable(labels.len() / 2);
            let (upper, lower) = (2 * *median, (*median / 2).max(CONS_MIN_LENGTH));
            let med_idx = labels.len() / 2;
            labels.select_nth_unstable_by_key(med_idx, |x| x.len());
            labels.swap(0, med_idx);
            let labels: Vec<&[u8]> = labels
                .iter()
                .filter(|ls| lower < ls.len() && ls.len() < upper)
                .map(|x| x.as_slice())
                .collect();
            if labels.len() <= cov_thr {
                return None;
            }
            let rough_contig = kiley::ternary_consensus_by_chunk(&labels, band_width);
            let rough_contig = kiley::bialignment::polish_until_converge_banded(
                &rough_contig,
                &labels,
                band_width,
            );
            match rough_contig.len() {
                0..=CONS_MIN_LENGTH => None,
                _ => Some((key, rough_contig)),
            }
        })
        .collect()
}

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

// Squish short contig, return the generated multi-tigs.
fn enumerate_multitigs(
    ds: &mut DataSet,
    &DenseEncodingConfig { len }: &DenseEncodingConfig,
) -> Vec<(usize, Vec<(u64, u64)>)> {
    use crate::assemble::*;
    let min_span_reads = match ds.read_type {
        ReadType::CCS => 1,
        ReadType::CLR => 4,
        ReadType::ONT => 1,
        ReadType::None => 3,
    };
    let config = AssembleConfig::new(1, 1000, false, true, min_span_reads);
    ds.squish_small_contig(&config, len);
    let (_, summaries) = assemble(ds, &config);
    let mut multi_tig: Vec<_> = summaries
        .iter()
        .filter(|summary| !summary.summary.is_empty())
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
    // TODO:Is this needed?
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
    units: &[Unit],
    read_type: &definitions::ReadType,
) -> Vec<definitions::Node> {
    let mut seq = if is_forward {
        seq[start..end].to_vec()
    } else {
        bio_utils::revcmp(&seq[start..end])
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let contig: Vec<_> = units.iter().map(|x| x.seq()).fold(Vec::new(), |mut x, y| {
        x.extend(y);
        x
    });
    let (break_points, _): (Vec<_>, _) =
        units
            .iter()
            .map(|x| x.seq().len())
            .fold((Vec::new(), 0), |(mut acc, x), y| {
                acc.push(x + y);
                (acc, x + y)
            });
    let (band, sim_thr) = (read_type.band_width(), read_type.sim_thr());
    let (_, ops) = kiley::bialignment::global_banded(&contig, &seq, 2, -5, -6, -1, band);
    // Split the alignment into encoded nodes.
    // Current position on the contg and the query.
    let (mut xpos, mut ypos) = (0, 0);
    // Current edit operations inside the focal unit.
    let mut alignments = vec![];
    let mut target_idx = 0;
    // Encoded nodes.
    let mut nodes = vec![];
    for op in ops {
        match op {
            kiley::bialignment::Op::Mat | kiley::bialignment::Op::Mism => {
                xpos += 1;
                ypos += 1;
            }
            kiley::bialignment::Op::Del => {
                xpos += 1;
            }
            kiley::bialignment::Op::Ins => {
                ypos += 1;
            }
        }
        alignments.push(op);
        if xpos == break_points[target_idx] {
            // Reached the boundary.
            let unit = &units[target_idx];
            let (uid, unitlen) = (unit.id, unit.seq().len());
            let ylen = alignments
                .iter()
                .filter(|&&x| x != kiley::bialignment::Op::Del)
                .count();
            // let edit_dist = alignments
            //     .iter()
            //     .filter(|&&x| x != kiley::bialignment::Op::Mat)
            //     .count();
            let cigar = crate::encode::compress_kiley_ops(&alignments);
            let indel_mism = alignments
                .iter()
                .map(|&op| 1 - 2 * (op == kiley::bialignment::Op::Mat) as i32);
            let max_indel = crate::encode::max_region(indel_mism).max(0) as usize;
            let unitlen = unitlen as f64;
            // let dist_thr = (unitlen * sim_thr).floor() as usize;
            let gap_thr = ((unitlen * crate::encode::INDEL_FRACTION).round() as usize)
                .max(crate::encode::MIN_INDEL_SIZE);
            let percent_identity = {
                let (aln, mat) = alignments.iter().fold((0, 0), |(aln, mat), &op| match op {
                    kiley::bialignment::Op::Mat => (aln + 1, mat + 1),
                    _ => (aln + 1, mat),
                });
                mat as f64 / aln as f64
            };
            //if max_indel < gap_thr && edit_dist < dist_thr {
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
