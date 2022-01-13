use definitions::DataSet;
use definitions::{Edge, EncodedRead, Node, Op, RawRead, Unit};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::*;
pub mod deletion_fill;
/// Expected upper bound of the disparsity between polished segment vs raw reads.
/// TODO: should be changed based on the read types.
pub const CLR_CTG_SIM: f64 = 0.25;
pub const CLR_CLR_SIM: f64 = 0.35;
/// This is a parameter only valid for last program.
/// TODO: As the Licencing issue, maybe we should remove these parameters as well as dependencies for last.
pub const MARGIN: usize = 50;
// Any alignment having deletion longer than ALLOWED_END_GAP would be discarded.
// Increasing this value would be more "abundant" encoding,
// but it would be problem in local clustering, which requiring almost all the
// unit is correctly globally aligned.
const ALLOWED_END_GAP: usize = 50;
// /// Any alignment having Insertion/Deletion longer than INDEL_THRESHOLD would be discarded.
// pub const INDEL_THRESHOLD: usize = 50;
/// Any alignment having Insertion/Deletion longer than unit.seq() * INDEL_FRACTION would be discarded.
/// for 2Kbp length, the threshold is 30bp.
pub const INDEL_FRACTION: f64 = 1f64 / 40f64;
// pub const INDEL_FRACTION: f64 = 0.015;
pub const MIN_INDEL_SIZE: usize = 10;
pub trait Encode {
    fn encode(&mut self, threads: usize, sim_thr: f64);
}

impl Encode for definitions::DataSet {
    fn encode(&mut self, threads: usize, sim_thr: f64) {
        encode_by_mm2(self, threads, sim_thr).unwrap();
        deletion_fill::correct_unit_deletion(self, sim_thr);
        debug!("Encoded {} reads.", self.encoded_reads.len());
        assert!(self.encoded_reads.iter().all(is_uppercase));
    }
}

pub fn encode_by_mm2(ds: &mut definitions::DataSet, p: usize, sim_thr: f64) -> std::io::Result<()> {
    let mm2 = mm2_alignment(ds, p)?;
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u)).collect();
    let alignments: Vec<_> = String::from_utf8_lossy(&mm2)
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(|a| a.tstart < ALLOWED_END_GAP && a.tlen - a.tend < ALLOWED_END_GAP)
        .filter_map(|aln| {
            // Check the max-indel.
            use bio_utils::sam;
            let tname = aln.tname.parse::<u64>().unwrap();
            let chunk_len = chunks.get(&tname)?.seq().len();
            let dist_thr = (chunk_len as f64 * sim_thr).floor() as usize;
            let gap_thr =
                ((chunk_len as f64 * INDEL_FRACTION).round() as usize).max(MIN_INDEL_SIZE);
            let cigar = sam::parse_cigar_string(aln.get_tag("cg")?.1);
            let indel_iter = cigar.iter().map(|op| match *op {
                sam::Op::Mismatch(l) | sam::Op::Deletion(l) | sam::Op::Insertion(l) => l as i32,
                sam::Op::Align(l) | sam::Op::Match(l) => -(l as i32),
                _ => 0,
            });
            let max_indel = max_region(indel_iter).max(0) as usize;
            // let no_large_indel = max_indel < INDEL_THRESHOLD;
            let no_large_indel = max_indel < gap_thr;
            let dist: usize = cigar
                .iter()
                .map(|op| match *op {
                    sam::Op::Deletion(l) | sam::Op::Insertion(l) | sam::Op::Mismatch(l) => l,
                    _ => 0,
                })
                .sum();
            (no_large_indel && dist < dist_thr).then(|| aln)
        })
        .collect();
    encode_by(ds, &alignments);
    Ok(())
}

pub fn encode_by(ds: &mut DataSet, alignments: &[bio_utils::paf::PAF]) {
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u)).collect();
    let mut bucket: HashMap<_, Vec<_>> = HashMap::new();
    for aln in alignments.iter() {
        bucket.entry(aln.qname.clone()).or_default().push(aln);
    }
    if log_enabled!(log::Level::Debug) {
        let no_aln_reads: Vec<_> = ds
            .raw_reads
            .iter()
            .filter(|read| bucket.get(&read.name).is_none())
            .collect();
        let lensum: usize = no_aln_reads.iter().map(|x| x.seq().len()).sum();
        debug!("ENCODE\tNotEncoded\t{}\t{}", lensum, no_aln_reads.len());
    }
    ds.encoded_reads = ds
        .raw_reads
        .par_iter()
        .filter_map(|read| {
            bucket
                .get(&read.name)
                .and_then(|alns| encode_read_by_paf(read, alns, &chunks))
        })
        .collect::<Vec<_>>();
}

fn encode_read_by_paf(
    read: &RawRead,
    alns: &[&bio_utils::paf::PAF],
    units: &HashMap<u64, &Unit>,
) -> Option<EncodedRead> {
    let mut seq: Vec<_> = read.seq().to_vec();
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let nodes = encode_read_to_nodes_by_paf(&seq, alns, units)?;
    nodes_to_encoded_read(read.id, nodes, &seq)
}

pub fn nodes_to_encoded_read(id: u64, nodes: Vec<Node>, seq: &[u8]) -> Option<EncodedRead> {
    let leading_gap = {
        let start_pos = nodes.first()?.position_from_start;
        seq[..start_pos].to_vec()
    };
    let trailing_gap = {
        let last_node = nodes.last()?;
        let end_pos = last_node.position_from_start + last_node.query_length();
        seq[end_pos..].to_vec()
    };
    let edges: Vec<_> = nodes.windows(2).map(|w| Edge::from_nodes(w, seq)).collect();
    Some(EncodedRead {
        id,
        original_length: seq.len(),
        leading_gap,
        trailing_gap,
        nodes,
        edges,
    })
}

fn encode_read_to_nodes_by_paf(
    seq: &[u8],
    alns: &[&bio_utils::paf::PAF],
    units: &HashMap<u64, &Unit>,
) -> Option<Vec<Node>> {
    let mut nodes: Vec<_> = alns
        .iter()
        .filter_map(|aln| {
            let tname = aln.tname.parse::<u64>().ok()?;
            encode_paf(seq, aln, units.get(&tname)?)
        })
        .collect();
    (!nodes.is_empty()).then(|| {
        nodes.sort_by_key(|n| n.unit);
        let mut nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        remove_slippy_alignment(nodes)
    })
}

fn encode_paf(seq: &[u8], aln: &bio_utils::paf::PAF, unit: &Unit) -> Option<Node> {
    use bio_utils::sam;
    let cigar = sam::parse_cigar_string(aln.get_tag("cg")?.1);
    // It is padded sequence. Should be adjusted.
    let (mut leading, aligned, mut trailing) = if aln.relstrand {
        let aligned = seq[aln.qstart..aln.qend].to_vec();
        let start = aln.qstart.saturating_sub(2 * aln.tstart);
        let end = (aln.qend + 2 * (aln.tlen - aln.tend)).min(seq.len());
        let leading = seq[start..aln.qstart].to_vec();
        let trailing = seq[aln.qend..end].to_vec();
        (leading, aligned, trailing)
    } else {
        let aligned = bio_utils::revcmp(&seq[aln.qstart..aln.qend]);
        let start = aln.qstart.saturating_sub(2 * (aln.tlen - aln.tend));
        let end = (aln.qend + 2 * aln.tstart).min(seq.len());
        let leading = bio_utils::revcmp(&seq[aln.qend..end]);
        let trailing = bio_utils::revcmp(&seq[start..aln.qstart]);
        (leading, aligned, trailing)
    };
    let mut ops = vec![];
    if aln.tstart > 0 {
        let (_, mut lops) =
            kiley::bialignment::edit_dist_slow_ops_semiglobal(&unit.seq()[..aln.tstart], &leading);
        // Leading Operations
        // Remove leading insertions.
        while lops[0] == kiley::bialignment::Op::Ins {
            leading.remove(0);
            lops.remove(0);
        }
        // We should have some element, as there are `aln.tstart` reference sequences, and this value is positive.
        assert!(!lops.is_empty());
        // Compress the resulting elements.
        ops.extend(compress_kiley_ops(&lops));
    } else {
        // If the alignment is propery starts from unit begining,
        // we do not need leading sequence.
        assert!(leading.is_empty());
    }
    // Main alignment.
    for op in cigar {
        let op = match op {
            sam::Op::Align(l) | sam::Op::Match(l) | sam::Op::Mismatch(l) => Op::Match(l),
            sam::Op::Insertion(l) => Op::Ins(l),
            sam::Op::Deletion(l) => Op::Del(l),
            _ => panic!("{:?}", op),
        };
        ops.push(op);
    }
    if aln.tlen != aln.tend {
        // Trailing operations
        let (_, mut tops) =
            kiley::bialignment::edit_dist_slow_ops_semiglobal(&unit.seq()[aln.tend..], &trailing);
        // Remove while trailing insertions.
        while tops.last() == Some(&kiley::bialignment::Op::Ins) {
            tops.pop();
            trailing.pop();
        }
        // We shoulud have remaining alignment.
        assert!(!tops.is_empty());
        ops.extend(compress_kiley_ops(&tops));
    } else {
        // Likewise, if the alignment propery stopped at the end of the units,
        // we do not need trailing sequence.
        assert!(trailing.is_empty());
    }
    let position_from_start = if aln.relstrand {
        aln.qstart - leading.len()
    } else {
        aln.qstart - trailing.len()
    };
    // Combine leading, aligned, and trailing sequences into one.
    leading.extend(aligned);
    leading.extend(trailing);
    let query_length = ops
        .iter()
        .map(|op| match op {
            Op::Ins(l) | Op::Match(l) => *l,
            _ => 0,
        })
        .sum::<usize>();
    assert_eq!(query_length, leading.len());
    let cl = unit.cluster_num;
    let node = Node::new(
        unit.id,
        aln.relstrand,
        &leading,
        ops,
        position_from_start,
        cl,
    );
    Some(node)
}

pub fn compress_kiley_ops(k_ops: &[kiley::bialignment::Op]) -> Vec<Op> {
    use kiley::bialignment;
    fn is_the_same(op1: bialignment::Op, op2: bialignment::Op) -> bool {
        use bialignment::Op::*;
        match (op1, op2) {
            (Mat, Mism) | (Mism, Mat) => true,
            (x, y) if x == y => true,
            _ => false,
        }
    }
    assert!(!k_ops.is_empty());
    let (mut current_op, mut len) = (k_ops[0], 1);
    let mut ops = vec![];
    for &op in k_ops.iter().skip(1) {
        if is_the_same(op, current_op) {
            len += 1;
        } else {
            match current_op {
                bialignment::Op::Del => ops.push(Op::Del(len)),
                bialignment::Op::Ins => ops.push(Op::Ins(len)),
                bialignment::Op::Mat => ops.push(Op::Match(len)),
                bialignment::Op::Mism => ops.push(Op::Match(len)),
            }
            current_op = op;
            len = 1;
        }
    }
    match current_op {
        bialignment::Op::Del => ops.push(Op::Del(len)),
        bialignment::Op::Ins => ops.push(Op::Ins(len)),
        bialignment::Op::Mat => ops.push(Op::Match(len)),
        bialignment::Op::Mism => ops.push(Op::Match(len)),
    }
    ops
}

pub fn remove_slippy_alignment(nodes: Vec<Node>) -> Vec<Node> {
    fn score(n: &Node) -> i32 {
        n.cigar
            .iter()
            .map(|op| match *op {
                Op::Match(x) => x as i32,
                Op::Del(x) | Op::Ins(x) => -(x as i32),
            })
            .sum::<i32>()
    }
    // Remove overlapping same units.
    let mut deduped_nodes = vec![];
    let mut nodes = nodes.into_iter();
    let mut prev = nodes.next().unwrap();
    for node in nodes {
        let is_disjoint = prev.position_from_start + prev.query_length() < node.position_from_start;
        if prev.unit != node.unit || prev.is_forward != node.is_forward || is_disjoint {
            deduped_nodes.push(prev);
            prev = node;
        } else if score(&prev) < score(&node) {
            prev = node;
        }
    }
    deduped_nodes.push(prev);
    deduped_nodes
}

pub fn mm2_alignment(ds: &definitions::DataSet, p: usize) -> std::io::Result<Vec<u8>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut c_dir = std::env::current_dir()?;
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir)?;
    // Create reference and reads.
    let (reference, reads) = {
        let mut reference = c_dir.clone();
        reference.push("units.fa");
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        for unit in ds.selected_chunks.iter() {
            writeln!(&mut wtr, ">{}\n{}", unit.id, &unit.seq)?;
        }
        let mut reads = c_dir.clone();
        reads.push("reads.fa");
        let mut wtr = std::fs::File::create(&reads).map(BufWriter::new)?;
        for read in ds.raw_reads.iter() {
            writeln!(&mut wtr, ">{}\n{}", read.name, &read.seq)?;
        }
        let reference = reference.into_os_string().into_string().unwrap();
        let reads = reads.into_os_string().into_string().unwrap();
        (reference, reads)
    };
    use crate::minimap2;
    let threads = format!("{}", p);
    use definitions::ReadType;
    let mut args = vec!["-t", &threads, "-c", "--eqx", "-P"];
    match ds.read_type {
        ReadType::CCS => args.extend(vec!["-H"]),
        ReadType::CLR => args.extend(vec!["-H", "-k", "15"]),
        ReadType::ONT => args.extend(vec!["-k", "15"]),
        _ => {}
    };
    let mm2 = minimap2::minimap2_args(&reference, &reads, &args);
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    Ok(mm2)
}

fn is_uppercase(read: &definitions::EncodedRead) -> bool {
    let nodes = read
        .nodes
        .iter()
        .all(|node| node.seq.chars().all(|c| c.is_ascii_uppercase()));
    let edges = read
        .edges
        .iter()
        .all(|edge| edge.label.chars().all(|c| c.is_ascii_uppercase()));
    nodes && edges
}

// The maximum value of sum of a range in xs,
// If the sequence is empty, return i32::MIN
pub fn max_region<T: std::iter::Iterator<Item = i32>>(xs: T) -> i32 {
    // The max value of the sum of the range ending at i.
    let mut right = i32::MIN;
    // The max value of the sum of the range ending before i.
    let mut left = i32::MIN;
    for x in xs {
        left = right.max(left);
        right = match right.is_negative() {
            true => x,
            false => right + x,
        };
    }
    right.max(left)
}

// pub fn encode(read: &RawRead, alignments: &[&LastTAB], units: &[Unit]) -> Option<EncodedRead> {
//     let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
//     for &aln in alignments.iter().filter(|aln| aln.seq1_matchlen() > MARGIN) {
//         let r_name = aln.seq1_name().to_string();
//         let q_direction = aln.seq2_direction().is_forward();
//         buckets.entry((r_name, q_direction)).or_default().push(aln);
//     }
//     let seq: Vec<_> = read.seq().iter().map(|c| c.to_ascii_uppercase()).collect();
//     let mut nodes: Vec<Node> = vec![];
//     for (_, alns) in buckets {
//         if let Some(mut ns) = encode_alignment(alns, units, &seq) {
//             nodes.append(&mut ns);
//         }
//     }
//     nodes.sort_by_key(|e| e.position_from_start);
//     let edges: Vec<_> = nodes
//         .windows(2)
//         .map(|w| Edge::from_nodes(w, &seq))
//         .collect();
//     let leading_gap: Vec<_> = {
//         let start = nodes.first()?.position_from_start;
//         seq[..start].to_vec()
//     };
//     let trailing_gap: Vec<_> = {
//         let end = nodes.last()?;
//         let end = end.position_from_start + end.seq.len();
//         seq[end..].to_vec()
//     };
//     let len = nodes.iter().map(|n| n.seq.len()).sum::<usize>() as i64;
//     let edge_len = edges.iter().map(|n| n.offset).sum::<i64>();
//     let chunked_len = (len + edge_len) as usize;
//     assert_eq!(
//         read.seq().len(),
//         chunked_len + leading_gap.len() + trailing_gap.len()
//     );
//     Some(EncodedRead {
//         original_length: read.seq().len(),
//         id: read.id,
//         edges,
//         nodes,
//         leading_gap,
//         trailing_gap,
//     })
// }

// fn encode_alignment(mut alns: Vec<&LastTAB>, units: &[Unit], seq: &[u8]) -> Option<Vec<Node>> {
//     let (is_forward, unit_id): (_, u64) = {
//         let aln = alns.get(0)?;
//         let is_forward = aln.seq2_direction().is_forward();
//         let unit_id = aln.seq1_name().parse().ok()?;
//         (is_forward, unit_id)
//     };
//     let seq = if is_forward {
//         seq.to_vec()
//     } else {
//         bio_utils::revcmp(seq)
//     };
//     let unit = units.iter().find(|u| u.id == unit_id)?;
//     let mut result = vec![];
//     while let Some((position_from_start, s, cigar)) =
//         encode_alignment_by_chaining(&mut alns, unit, &seq)
//     {
//         let node = if is_forward {
//             Node {
//                 position_from_start,
//                 unit: unit_id,
//                 cluster: 0,
//                 seq: s,
//                 is_forward,
//                 cigar,
//             }
//         } else {
//             let position_from_start = seq.len() - position_from_start - s.len();
//             Node {
//                 position_from_start,
//                 unit: unit_id,
//                 cluster: 0,
//                 seq: s,
//                 is_forward,
//                 cigar,
//             }
//         };
//         result.push(node);
//     }
//     Some(result)
// }

// #[derive(Debug, Clone)]
// enum DAGNode<'a> {
//     Aln(&'a LastTAB),
//     Start,
//     End,
// }

// Query start, Query sequence, Cigar.
// type Alignment = (usize, String, Vec<Op>);
// fn encode_alignment_by_chaining(
//     alns: &mut Vec<&LastTAB>,
//     unit: &Unit,
//     read: &[u8],
// ) -> Option<Alignment> {
//     if alns.is_empty() {
//         return None;
//     };
//     let unitseq: Vec<_> = unit.seq().iter().map(|x| x.to_ascii_uppercase()).collect();
//     let mut dag_nodes = vec![DAGNode::Start, DAGNode::End];
//     dag_nodes.extend(alns.iter().map(|&a| DAGNode::Aln(a)));
//     let chain = chaining(&dag_nodes, unit.seq().len());
//     let (query_start, query_end) = {
//         let chain = chain.iter().filter_map(|n| match n {
//             DAGNode::Aln(aln) => Some(aln),
//             _ => None,
//         });
//         // First aln.
//         let start = chain.clone().next().unwrap().seq2_start();
//         // Last aln.
//         let end_aln = chain.clone().last().unwrap();
//         let end = end_aln.seq2_start() + end_aln.seq2_matchlen();
//         (start, end)
//     };
//     // Check if this chain covers sufficient fraction of the unit.
//     let cover_length = chain
//         .iter()
//         .map(|chain| match chain {
//             DAGNode::Aln(x) => x.seq1_matchlen(),
//             _ => 0,
//         })
//         .sum::<usize>();
//     if cover_length < unit.seq().len() * 9 / 10 {
//         return None;
//     }
//     let (mut q_pos, mut r_pos) = (query_start, 0);
//     let mut cigar = vec![];
//     assert!(matches!(*chain[0], DAGNode::Start));
//     for node in chain.iter().skip(1) {
//         let (q_target, r_target) = match node {
//             &&DAGNode::End => (query_end, unit.seq().len()),
//             DAGNode::Aln(x) => (x.seq2_start(), x.seq1_start()),
//             _ => panic!(),
//         };
//         assert!(q_pos <= q_target);
//         if r_pos <= r_target {
//             let query = &read[q_pos..q_target];
//             let refr = &unitseq[r_pos..r_target];
//             cigar.extend(alignment(query, refr, 1, -1, -1));
//         } else {
//             // Step back.
//             let query_pop_len = pop_cigar_by(&mut cigar, r_pos - r_target);
//             cigar.push(Op::Ins(q_target - q_pos + query_pop_len));
//         }
//         if let DAGNode::Aln(x) = node {
//             cigar.extend(convert_aln_to_cigar(x));
//             r_pos = x.seq1_start() + x.seq1_matchlen();
//             q_pos = x.seq2_start() + x.seq2_matchlen();
//         }
//     }
//     cigar.reverse();
//     let mut query_start = query_start;
//     while let Some(&Op::Ins(l)) = cigar.last() {
//         query_start += l;
//         cigar.pop();
//     }
//     cigar.reverse();
//     let mut query_end = query_end;
//     while let Some(&Op::Ins(l)) = cigar.last() {
//         query_end -= l;
//         cigar.pop();
//     }
//     assert_eq!(consumed_reference_length(&cigar), unit.seq().len());
//     let query = String::from_utf8_lossy(&read[query_start..query_end]).to_string();
//     // Remove used alignment.
//     let chain: Vec<&LastTAB> = chain
//         .iter()
//         .filter_map(|n| match n {
//             DAGNode::Aln(a) => Some(*a),
//             _ => None,
//         })
//         .collect();
//     alns.retain(|&aln| {
//         chain
//             .iter()
//             .all(|&a| a as *const LastTAB != aln as *const LastTAB)
//     });
//     Some((query_start, query, cigar))
// }

// /// Public interface.
// pub fn join_alignments(
//     alns: &mut Vec<&LastTAB>,
//     refr: &[u8],
//     read: &[u8],
// ) -> (usize, usize, Vec<Op>) {
//     assert!(!alns.is_empty());
//     let mut dag_nodes = vec![DAGNode::Start, DAGNode::End];
//     dag_nodes.extend(alns.iter().map(|&a| DAGNode::Aln(a)));
//     let chain = chaining(&dag_nodes, refr.len());
//     let (query_start, refr_start) = {
//         let first_chain = chain
//             .iter()
//             .find_map(|n| match n {
//                 DAGNode::Aln(aln) => Some(aln),
//                 _ => None,
//             })
//             .unwrap();
//         // First aln.
//         let query_start = first_chain.seq2_start();
//         let refr_start = first_chain.seq1_start();
//         (query_start, refr_start)
//     };
//     let (mut q_pos, mut r_pos) = (query_start, refr_start);
//     let mut cigar = vec![];
//     assert!(matches!(*chain[0], DAGNode::Start));
//     for node in chain.iter().skip(1) {
//         let (q_target, r_target) = match node {
//             &&DAGNode::End => break,
//             DAGNode::Aln(x) => (x.seq2_start(), x.seq1_start()),
//             _ => panic!(),
//         };
//         assert!(q_pos <= q_target);
//         if r_pos <= r_target {
//             // Usual chaining
//             let (query, refr) = (&read[q_pos..q_target], &refr[r_pos..r_target]);
//             cigar.extend(alignment(query, refr, 1, -1, -1));
//         } else {
//             // Step back a little bit.
//             let query_pop_len = pop_cigar_by(&mut cigar, r_pos - r_target);
//             cigar.push(Op::Ins(q_target - q_pos + query_pop_len));
//         }
//         if let DAGNode::Aln(x) = node {
//             cigar.extend(convert_aln_to_cigar(x));
//             r_pos = x.seq1_start() + x.seq1_matchlen();
//             q_pos = x.seq2_start() + x.seq2_matchlen();
//         }
//     }
//     let chain: Vec<&LastTAB> = chain
//         .iter()
//         .filter_map(|node| match node {
//             DAGNode::Aln(x) => Some(*x),
//             _ => None,
//         })
//         .collect();
//     alns.retain(|&aln| {
//         chain
//             .iter()
//             .all(|&a| a as *const LastTAB != aln as *const LastTAB)
//     });
//     (query_start, refr_start, cigar)
// }

// fn pop_cigar_by(cigar: &mut Vec<Op>, ref_len: usize) -> usize {
//     assert!(ref_len > 0);
//     let mut query_pop_len = 0;
//     let mut refr_pop_len = 0;
//     let mut op = None;
//     while refr_pop_len < ref_len {
//         op = cigar.pop();
//         match op {
//             Some(Op::Del(l)) => refr_pop_len += l,
//             Some(Op::Ins(l)) => query_pop_len += l,
//             Some(Op::Match(l)) => {
//                 refr_pop_len += l;
//                 query_pop_len += l;
//             }
//             None => panic!("{}\t{}\t{}", ref_len, refr_pop_len, line!()),
//         }
//     }
//     let overflow = refr_pop_len - ref_len;
//     if overflow > 0 {
//         match op {
//             Some(Op::Del(_)) => cigar.push(Op::Del(overflow)),
//             Some(Op::Match(_)) => {
//                 assert!(query_pop_len >= overflow);
//                 query_pop_len -= overflow;
//                 cigar.push(Op::Match(overflow))
//             }
//             _ => panic!("{}", line!()),
//         }
//     }
//     query_pop_len
// }

// pub fn recover(query: &[u8], refr: &[u8], ops: &[Op]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
//     let (mut q, mut al, mut r) = (vec![], vec![], vec![]);
//     let (mut q_pos, mut r_pos) = (0, 0);
//     for op in ops {
//         match *op {
//             Op::Match(l) => {
//                 al.extend(
//                     query[q_pos..q_pos + l]
//                         .iter()
//                         .zip(&refr[r_pos..r_pos + l])
//                         .map(|(x, y)| if x == y { b'|' } else { b'X' }),
//                 );
//                 q.extend(query[q_pos..q_pos + l].iter().copied());
//                 r.extend(refr[r_pos..r_pos + l].iter().copied());
//                 q_pos += l;
//                 r_pos += l;
//             }
//             Op::Del(l) => {
//                 al.extend(vec![b' '; l]);
//                 q.extend(vec![b' '; l]);
//                 r.extend(refr[r_pos..r_pos + l].iter().copied());
//                 r_pos += l;
//             }
//             Op::Ins(l) => {
//                 al.extend(vec![b' '; l]);
//                 q.extend(query[q_pos..q_pos + l].iter().copied());
//                 r.extend(vec![b' '; l]);
//                 q_pos += l;
//             }
//         }
//     }
//     (q, al, r)
// }

// #[allow(dead_code)]
// fn get_read_range(alns: &[&LastTAB], unit: &Unit, read: &[u8]) -> Option<(usize, usize)> {
//     let query_start = {
//         let first = alns.iter().min_by_key(|aln| aln.seq1_start())?;
//         let remaining = first.seq1_start();
//         let query_position = first.seq2_start();
//         query_position.max(2 * remaining) - 2 * remaining
//     };
//     let query_end = {
//         let last = alns
//             .iter()
//             .max_by_key(|aln| aln.seq1_start() + aln.seq1_matchlen())?;
//         let remaining = unit.seq().len() - (last.seq1_start() + last.seq1_matchlen());
//         let query_position = last.seq2_start() + last.seq2_matchlen();
//         (query_position + 2 * remaining).min(read.len())
//     };
//     Some((query_start, query_end))
// }

// fn chaining<'a, 'b>(nodes: &'b [DAGNode<'a>], unit_len: usize) -> Vec<&'b DAGNode<'a>> {
//     let edges: Vec<Vec<(usize, i64)>> = nodes
//         .iter()
//         .map(|u| {
//             nodes
//                 .iter()
//                 .enumerate()
//                 .filter_map(|(idx, v)| compute_edge(u, v, unit_len).map(|w| (idx, w)))
//                 .collect()
//         })
//         .collect();
//     let order = topological_sort(&edges);
//     // By initializing by 0, we assume that we can start anywhere...
//     let (mut dp, mut parent) = (vec![-1; nodes.len()], vec![0; nodes.len()]);
//     for i in order {
//         for &(j, score) in edges[i].iter() {
//             let aln_j = match nodes[j] {
//                 DAGNode::Aln(aln) => aln.score() as i64,
//                 _ => 0,
//             };
//             let from_i_to_j = dp[i] + score + aln_j;
//             if dp[j] <= from_i_to_j {
//                 dp[j] = from_i_to_j;
//                 parent[j] = i;
//             }
//         }
//     }
//     // i <= 1, tracing back from the end node.
//     assert!(matches!(nodes[1], DAGNode::End));
//     let mut path = vec![];
//     let mut current_node = 1;
//     while current_node != 0 {
//         path.push(&nodes[current_node]);
//         current_node = parent[current_node];
//     }
//     path.push(&nodes[current_node]);
//     path.reverse();
//     path
// }

// // Compute edge score between from u to v. If no edge possible, return None.
// // Caution: It accepts `slippy` alignments, in other words,
// // sometimes an edge would be drawn even if aln1's end position is larger than
// // aln2's start position.
// fn compute_edge<'a>(u: &DAGNode<'a>, v: &DAGNode<'a>, _unit_len: usize) -> Option<i64> {
//     match u {
//         DAGNode::End => None,
//         DAGNode::Start => Some(0),
//         DAGNode::Aln(aln1) => match v {
//             DAGNode::End => Some(0),
//             DAGNode::Start => None,
//             DAGNode::Aln(aln2) => {
//                 let u_refr_end = (aln1.seq1_start() + aln1.seq1_matchlen()).max(MARGIN) - MARGIN;
//                 let u_query_end = aln1.seq2_start() + aln1.seq2_matchlen();
//                 let u_refr_end_radius = u_refr_end + MARGIN;
//                 let u_query_end_radius = u_query_end + MARGIN;
//                 let v_refr_start = aln2.seq1_start();
//                 let v_query_start = aln2.seq2_start();
//                 if u_refr_end <= v_refr_start
//                     && v_refr_start < u_refr_end_radius
//                     && u_query_end <= v_query_start
//                     && v_query_start < u_query_end_radius
//                 {
//                     Some(-((v_refr_start - u_refr_end + v_query_start - u_query_end) as i64))
//                 } else {
//                     None
//                 }
//             }
//         },
//     }
// }

// // Return nodes in topological order. Note that the graph IS connected.
// fn topological_sort(edges: &[Vec<(usize, i64)>]) -> Vec<usize> {
//     let len = edges.len();
//     let mut arrived = vec![false; len];
//     // 0 is the start node.
//     let mut stack = vec![0];
//     let mut order = vec![];
//     'dfs: while !stack.is_empty() {
//         let node = *stack.last().unwrap();
//         if !arrived[node] {
//             arrived[node] = true;
//         }
//         for &(to, _) in &edges[node] {
//             if !arrived[to] {
//                 stack.push(to);
//                 continue 'dfs;
//             }
//         }
//         let last = stack.pop().unwrap();
//         order.push(last);
//     }
//     order.reverse();
//     order
// }

// fn alignment(query: &[u8], refr: &[u8], mat: i64, mism: i64, gap: i64) -> Vec<Op> {
//     assert!(query.iter().all(|c| c.is_ascii_uppercase()));
//     assert!(refr.iter().all(|c| c.is_ascii_uppercase()));
//     let mut dp = vec![vec![0; query.len() + 1]; refr.len() + 1];
//     for j in 0..=query.len() {
//         dp[0][j] = j as i64 * gap;
//     }
//     for (i, row) in dp.iter_mut().enumerate() {
//         row[0] = i as i64 * gap;
//     }
//     for (i, r) in refr.iter().enumerate() {
//         for (j, q) in query.iter().enumerate() {
//             let match_score = if r == q { mat } else { mism };
//             let max = (dp[i][j] + match_score)
//                 .max(dp[i][j + 1] + gap)
//                 .max(dp[i + 1][j] + gap);
//             dp[i + 1][j + 1] = max;
//         }
//     }
//     // Traceback.
//     let (mut q_pos, mut r_pos) = (query.len(), refr.len());
//     let mut ops = vec![];
//     while 0 < q_pos && 0 < r_pos {
//         let current = dp[r_pos][q_pos];
//         if current == dp[r_pos - 1][q_pos] + gap {
//             ops.push(2);
//             r_pos -= 1;
//         } else if current == dp[r_pos][q_pos - 1] + gap {
//             ops.push(1);
//             q_pos -= 1;
//         } else {
//             let match_score = if query[q_pos - 1] == refr[r_pos - 1] {
//                 mat
//             } else {
//                 mism
//             };
//             assert_eq!(current, dp[r_pos - 1][q_pos - 1] + match_score);
//             ops.push(0);
//             q_pos -= 1;
//             r_pos -= 1;
//         }
//     }
//     ops.extend(std::iter::repeat(1).take(q_pos));
//     ops.extend(std::iter::repeat(2).take(r_pos));
//     // for _ in 0..q_pos {
//     //     ops.push(1);
//     // }
//     // for _ in 0..r_pos {
//     //     ops.push(2);
//     // }
//     compress(ops)
// }

// fn compress(mut ops: Vec<u8>) -> Vec<Op> {
//     let mut cigar = vec![];
//     while !ops.is_empty() {
//         let last = ops.pop().unwrap();
//         let mut count = 1;
//         while let Some(&res) = ops.last() {
//             if res == last {
//                 count += 1;
//                 ops.pop();
//             } else {
//                 break;
//             }
//         }
//         match last {
//             0 => cigar.push(Op::Match(count)),
//             1 => cigar.push(Op::Ins(count)),
//             2 => cigar.push(Op::Del(count)),
//             _ => panic!(),
//         }
//     }
//     cigar
// }

// fn convert_aln_to_cigar(aln: &lasttab::LastTAB) -> Vec<Op> {
//     aln.alignment()
//         .iter()
//         .map(|op| match op {
//             lasttab::Op::Seq1In(l) => Op::Ins(*l),
//             lasttab::Op::Seq2In(l) => Op::Del(*l),
//             lasttab::Op::Match(l) => Op::Match(*l),
//         })
//         .collect()
// }

// fn consumed_reference_length(cigar: &[Op]) -> usize {
//     cigar
//         .iter()
//         .map(|op| match op {
//             Op::Match(l) | Op::Del(l) => *l,
//             Op::Ins(_) => 0,
//         })
//         .sum::<usize>()
// }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn works() {}
    #[test]
    fn max_range_operation_test() {
        let ops = [
            Op::Match(10),
            Op::Del(5),
            Op::Ins(2),
            Op::Match(1),
            Op::Del(5),
            Op::Match(10),
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Ins(l) => *l as i32 * -1,
            Op::Del(l) => *l as i32 * 2,
        });
        let max_del = max_region(iter);
        assert_eq!(max_del, 10 + 10 - 3);
        let ops = [
            Op::Ins(10),
            Op::Del(5),
            Op::Ins(2),
            Op::Match(1),
            Op::Del(5),
            Op::Match(10),
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Del(l) => *l as i32 * -1,
            Op::Ins(l) => *l as i32 * 2,
        });
        let max_in = max_region(iter);
        assert_eq!(max_in, 20);
        let ops = [
            Op::Ins(10),
            Op::Del(5),
            Op::Ins(2),    // 19
            Op::Match(1),  // 18
            Op::Del(5),    // 13
            Op::Match(10), // 3
            Op::Ins(100),  // 203
        ];
        let iter = ops.iter().map(|x| match x {
            Op::Match(l) | Op::Del(l) => *l as i32 * -1,
            Op::Ins(l) => *l as i32 * 2,
        });
        let max_in = max_region(iter);
        assert_eq!(max_in, 203);
    }
    // #[test]
    // fn alignment_check() {
    //     let query = b"AAAAA";
    //     let reference = b"AAAAA";
    //     let res = alignment(query, reference, 1, -1, -1);
    //     assert_eq!(res, vec![Op::Match(query.len())]);
    //     let query = b"AAAAA";
    //     let reference = b"AACAA";
    //     let res = alignment(query, reference, 1, -1, -1);
    //     assert_eq!(res, vec![Op::Match(query.len())]);
    //     let query = b"AACCAAAA";
    //     let refer = b"AAGCAA";
    //     let res = alignment(query, refer, 1, -1, -1);
    //     assert_eq!(res, vec![Op::Match(refer.len()), Op::Ins(2)]);
    //     let query = b"ACGCGCGCAA";
    //     let refer = b"GCGCGC";
    //     let res = alignment(query, refer, 1, -1, -1);
    //     assert_eq!(res, vec![Op::Ins(2), Op::Match(6), Op::Ins(2)]);
    //     let query = b"GCGCGC";
    //     let refer = b"ACGCGCGCAA";
    //     let res = alignment(query, refer, 1, -1, -1);
    //     assert_eq!(res, vec![Op::Del(2), Op::Match(6), Op::Del(2)]);
    //     let query = b"CGCTGCGCAAAAA";
    //     let refer = b"AAAAAGCGCGCT";
    //     let res = alignment(query, refer, 1, -1, -1);
    //     let ans = vec![
    //         Op::Match(1),
    //         Op::Del(4),
    //         Op::Match(2),
    //         Op::Ins(1),
    //         Op::Match(5),
    //         Op::Ins(4),
    //     ];
    //     assert_eq!(res, ans)
    // }
}
