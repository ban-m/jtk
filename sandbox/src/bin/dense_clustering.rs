use rayon::prelude::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    use definitions::*;
    let mut ds: DataSet = std::fs::File::open(&args[1])
        .map(BufReader::new)
        .map(|x| serde_json::de::from_reader(x).unwrap())?;
    use haplotyper::assemble::*;
    let config = AssembleConfig::new(1, 1000, false, true);
    ds.squish_small_contig(&config);
    let (_records, summaries) = assemble(&ds, 0, &config);
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
    let (cp, nodes) = multi_tig.remove(0);
    let cluster_num = cp;
    log::debug!("CONTIG\t{:?}", nodes);
    let mut reads: Vec<_> = ds
        .encoded_reads
        .iter_mut()
        .filter(|r| r.nodes.iter().any(|n| nodes.contains(&(n.unit, n.cluster))))
        .collect();
    // Create new nodes between nodes to make this region "dense."
    // (unit,cluster,is_tail).
    type Edge = ((u64, u64, bool), (u64, u64, bool));
    use std::collections::HashMap;
    let mut edges: HashMap<Edge, Vec<Vec<u8>>> = HashMap::new();
    for read in reads.iter() {
        for (window, edge) in read.nodes.windows(2).zip(read.edges.iter()) {
            // Determine direction. From smaller unit, true strand to the other.
            let (from, to) = match window {
                [from, to] => (from, to),
                _ => panic!(),
            };
            if !(nodes.contains(&(from.unit, from.cluster))
                && nodes.contains(&(to.unit, to.cluster)))
            {
                continue;
            }
            if edge.offset < 0 {
                continue;
            }
            let label = edge.label();
            let (edge_entry, mut label) = match from.unit.cmp(&to.unit) {
                std::cmp::Ordering::Less => {
                    let from_elm = (from.unit, from.cluster, from.is_forward);
                    let to_elm = (to.unit, to.cluster, !to.is_forward);
                    ((from_elm, to_elm), label.to_vec())
                }
                std::cmp::Ordering::Equal if from.is_forward => {
                    let from_elm = (from.unit, from.cluster, from.is_forward);
                    let to_elm = (to.unit, to.cluster, !to.is_forward);
                    ((from_elm, to_elm), label.to_vec())
                }
                _ => {
                    let from_elm = (to.unit, to.cluster, !to.is_forward);
                    let to_elm = (from.unit, from.cluster, from.is_forward);
                    ((from_elm, to_elm), bio_utils::revcmp(label))
                }
            };
            // Register.
            label.iter_mut().for_each(u8::make_ascii_uppercase);
            edges.entry(edge_entry).or_insert(vec![]).push(label);
        }
    }
    for (edge, ls) in edges.iter() {
        log::debug!("EDGE\t{:?}\t{}", edge, ls.len());
    }
    edges.retain(|_, val| 2 < val.len());
    let ave_unit_len = {
        let total: usize = ds.selected_chunks.iter().map(|x| x.seq().len()).sum();
        total / ds.selected_chunks.len()
    };
    // Create new units.
    let mut max_unit_id: u64 = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    let edge_consensus: HashMap<Edge, Vec<u8>> = edges
        .into_par_iter()
        .filter_map(|(key, labels)| {
            // Take consensus(Rough, on labels)
            let rough_contig = kiley::ternary_consensus_by_chunk(&labels, 100);
            match rough_contig.len() {
                0..=100 => None,
                _ => Some((key, rough_contig)),
            }
        })
        .collect();
    let edge_encoding_patterns: HashMap<Edge, Vec<(usize, u64)>> = edge_consensus
        .iter()
        .map(|(&key, consensus)| {
            let break_positions: Vec<_> = (1..)
                .map(|i| i * ave_unit_len)
                .take_while(|&break_pos| break_pos + ave_unit_len < consensus.len())
                .chain(std::iter::once(consensus.len()))
                .map(|break_pos| {
                    max_unit_id += 1;
                    (break_pos, (max_unit_id))
                })
                .collect();
            (key, break_positions)
        })
        .collect();
    for (key, new_units) in edge_encoding_patterns.iter() {
        let contig = &edge_consensus[key];
        let mut prev_bp = 0;
        log::debug!("{:?}\t{:?}\t{}", key, new_units, contig.len());
        for &(break_point, id) in new_units {
            let seq = String::from_utf8_lossy(&contig[prev_bp..break_point]).to_string();
            let unit = Unit::new(id, seq, cluster_num);
            ds.selected_chunks.push(unit);
            prev_bp = break_point;
        }
    }
    let read_seq: HashMap<u64, &[u8]> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    use std::collections::HashSet;
    let filled_reads: HashSet<_> = reads.iter().map(|r| r.id).collect();
    // Encoding.
    for read in reads.iter_mut() {
        let seq = read_seq[&read.id];
        let mut inserts = vec![];
        for (idx, window) in read.nodes.windows(2).enumerate() {
            let (from, to) = match window {
                [from, to] => (from, to),
                _ => panic!(),
            };
            // What is important is the direction, the label, and the start-end position. Thats all.
            // let label = edge.label();
            let (start, end) = (
                from.position_from_start + from.seq().len(),
                to.position_from_start,
            );
            let (edge, direction) = match from.unit.cmp(&to.unit) {
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
            };
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
        let mut accum_inserts = 0;
        for (idx, node) in inserts {
            read.nodes.insert(idx + accum_inserts, node);
            accum_inserts += 1;
        }
        if !read.nodes.is_empty() {
            let mut nodes = vec![];
            nodes.append(&mut read.nodes);
            use haplotyper::encode::{nodes_to_encoded_read, remove_slippy_alignment};
            nodes.sort_by_key(|n| n.unit);
            nodes = remove_slippy_alignment(nodes);
            nodes.sort_by_key(|n| n.position_from_start);
            nodes = remove_slippy_alignment(nodes);
            **read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
        }
    }
    // Filling up deletions.
    haplotyper::encode::deletion_fill::correct_unit_deletion_selected(&mut ds, &filled_reads);
    // nodes to be clustered.
    let nodes: HashSet<_> = {
        let mut nodes: HashSet<_> = nodes.iter().map(|x| x.0).collect();
        edge_encoding_patterns
            .values()
            .for_each(|new_units| nodes.extend(new_units.iter().map(|(_, id)| *id)));
        nodes
    };
    log::debug!("Cluster\t{:?}", nodes);
    // Local clustering.
    haplotyper::local_clustering::local_clustering_selected(&mut ds, &nodes);
    // Squish bad mods.
    let mut new_clustered: HashMap<_, Vec<&mut _>> = nodes.iter().map(|&n| (n, vec![])).collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        if let Some(bucket) = new_clustered.get_mut(&node.unit) {
            bucket.push(&mut node.cluster);
        }
    }
    for unit in ds.selected_chunks.iter_mut() {
        if let Some(assignments) = new_clustered.get_mut(&unit.id) {
            // At least 1 LK for each element(CLRmode)
            let len = assignments.len();
            if unit.score <= assignments.len() as f64 {
                log::debug!("Squishing\t{}\t{}\t{}", unit.id, unit.score, len,);
                unit.score = 0f64;
                assignments.iter_mut().for_each(|x| **x = 0);
            }
        }
    }
    // EM clustering.
    use haplotyper::em_correction::ClusteringCorrection;
    ds = ds.correct_clustering_em_on_selected(10, 3, false, &nodes);
    ds = ds.correct_clustering_em(10, 3, false);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
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
            let cigar = haplotyper::encode::compress_kiley_ops(&alignments);
            use definitions::Op;
            let has_large_indel = cigar.iter().any(|op| match *op {
                Op::Ins(l) | Op::Del(l) => l > haplotyper::encode::INDEL_THRESHOLD,
                _ => false,
            });
            let dist_thr = (unitlen as f64 * 0.35).floor() as u32;
            if !has_large_indel && edit_dist < dist_thr {
                use definitions::Node;
                // Encoding.
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
