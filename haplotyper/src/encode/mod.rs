use definitions::DataSet;
use definitions::{Chunk, Edge, EncodedRead, Node, Op, RawRead};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::*;
pub mod deletion_fill;
pub const MARGIN: usize = 50;
use log::*;

// Any alignment having deletion longer than ALLOWED_END_GAP would be discarded.
// Increasing this value would be more "abundant" encoding,
// but it would be problem in local clustering, which requiring almost all the
// chunk is correctly globally aligned.
const ALLOWED_END_GAP: usize = 25;
pub trait Encode {
    fn encode(&mut self, threads: usize, sim_thr: f64, sd_of_error: f64);
}

impl Encode for definitions::DataSet {
    fn encode(&mut self, threads: usize, sim_thr: f64, sd_of_error: f64) {
        debug!("ENCODE\tErrorRate\t{sim_thr}");
        encode_by_mm2(self, threads, sim_thr).unwrap();
        let config =
            deletion_fill::CorrectDeletionConfig::new(false, Some(sim_thr), Some(sd_of_error));
        deletion_fill::correct_chunk_deletion(self, &config);
        // debug!("Encoded {} reads.", self.encoded_reads.len());
        assert!(self.encoded_reads.iter().all(is_uppercase));
        if log_enabled!(log::Level::Debug) {
            let encoded: HashSet<_> = self.encoded_reads.iter().map(|r| r.id).collect();
            let no_aln_reads: Vec<_> = self
                .raw_reads
                .iter()
                .filter(|&read| (!encoded.contains(&read.id))).map(|read| read.seq().len())
                .collect();
            let lensum: usize = no_aln_reads.iter().sum();
            debug!("ENCODE\tNotEncoded\t{}\t{}", lensum, no_aln_reads.len());
        }
    }
}

pub fn encode_by_mm2(ds: &mut definitions::DataSet, p: usize, sim_thr: f64) -> std::io::Result<()> {
    let mm2 = mm2_alignment(ds, p)?;
    let alignments: Vec<_> = String::from_utf8_lossy(&mm2)
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(|a| a.tstart < ALLOWED_END_GAP && a.tlen - a.tend < ALLOWED_END_GAP)
        .filter_map(|aln| {
            use bio_utils::sam;
            let cigar = sam::parse_cigar_string(aln.get_tag("cg")?.1);
            let (aln_len, mat_num) = cigar.iter().fold((0, 0), |(aln, mat), &op| match op {
                sam::Op::Align(x)
                | sam::Op::Insertion(x)
                | sam::Op::Deletion(x)
                | sam::Op::Mismatch(x) => (aln + x, mat),
                sam::Op::Match(x) => (aln + x, mat + x),
                _ => unreachable!(),
            });
            let percent_identity = mat_num as f64 / aln_len as f64;
            (1f64 - sim_thr < percent_identity).then_some(aln)
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
    chunks: &HashMap<u64, &Chunk>,
) -> Option<EncodedRead> {
    let mut seq: Vec<_> = read.seq().to_vec();
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let nodes = encode_read_to_nodes_by_paf(&seq, alns, chunks)?;
    nodes_to_encoded_read(read.id, nodes, &seq)
}

pub fn nodes_to_encoded_read(id: u64, nodes: Vec<Node>, seq: &[u8]) -> Option<EncodedRead> {
    let leading_gap: Vec<_> = {
        let start_pos = nodes.first()?.position_from_start;
        seq.iter()
            .take(start_pos)
            .map(u8::to_ascii_uppercase)
            .collect()
    };
    let trailing_gap: Vec<_> = {
        let last_node = nodes.last()?;
        let end_pos = last_node.position_from_start + last_node.query_length();
        seq.iter()
            .skip(end_pos)
            .map(u8::to_ascii_uppercase)
            .collect()
    };
    let edges: Vec<_> = nodes.windows(2).map(|w| Edge::from_nodes(w, seq)).collect();
    Some(EncodedRead {
        id,
        original_length: seq.len(),
        leading_gap: leading_gap.into(),
        trailing_gap: trailing_gap.into(),
        nodes,
        edges,
    })
}

fn encode_read_to_nodes_by_paf(
    seq: &[u8],
    alns: &[&bio_utils::paf::PAF],
    chunks: &HashMap<u64, &Chunk>,
) -> Option<Vec<Node>> {
    let mut nodes: Vec<_> = alns
        .iter()
        .filter_map(|aln| {
            let tname = aln.tname.parse::<u64>().ok()?;
            encode_paf(seq, aln, chunks.get(&tname)?)
        })
        .collect();
    (!nodes.is_empty()).then(|| {
        nodes.sort_by_key(|n| n.chunk);
        let mut nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        remove_slippy_alignment(nodes)
    })
}

fn split_query(seq: &[u8], aln: &bio_utils::paf::PAF) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    if aln.relstrand {
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
    }
}

fn leading_alignment(refr: &[u8], mut leading: Vec<u8>) -> (Vec<Op>, Vec<u8>) {
    let mut lops = semiglobal(refr, &leading);
    lops.reverse();
    leading.reverse();
    while lops.last() == Some(&kiley::Op::Ins) {
        assert!(leading.pop().is_some());
        assert!(lops.pop().is_some());
    }
    lops.reverse();
    leading.reverse();
    (crate::misc::kiley_op_to_ops(&lops).0, leading)
}

fn trailing_alignment(refr: &[u8], mut trailing: Vec<u8>) -> (Vec<Op>, Vec<u8>) {
    let mut tops = semiglobal(refr, &trailing);
    while tops.last() == Some(&kiley::Op::Ins) {
        assert!(tops.pop().is_some());
        assert!(trailing.pop().is_some());
    }
    (crate::misc::kiley_op_to_ops(&tops).0, trailing)
}

fn encode_paf(seq: &[u8], aln: &bio_utils::paf::PAF, chunk: &Chunk) -> Option<Node> {
    use bio_utils::sam;
    let cigar = sam::parse_cigar_string(aln.get_tag("cg")?.1);
    let (leading, aligned, trailing) = split_query(seq, aln);
    let mut ops = vec![];
    assert!(0 < aln.tstart || leading.is_empty());
    let (leading_aln, leading) = leading_alignment(&chunk.seq()[..aln.tstart], leading);
    ops.extend(leading_aln);
    let cigar = cigar.iter().map(|&op| match op {
        sam::Op::Align(l) | sam::Op::Match(l) | sam::Op::Mismatch(l) => Op::Match(l),
        sam::Op::Insertion(l) => Op::Ins(l),
        sam::Op::Deletion(l) => Op::Del(l),
        _ => panic!("{:?}", op),
    });
    ops.extend(cigar);
    assert!(aln.tlen != aln.tend || trailing.is_empty());
    let (trailing_aln, trailing) = trailing_alignment(&chunk.seq()[aln.tend..], trailing);
    ops.extend(trailing_aln);
    let position_from_start = match aln.relstrand {
        true => aln.qstart - leading.len(),
        false => aln.qstart - trailing.len(),
    };
    let seq = [leading, aligned, trailing].concat();
    check_length(&ops, seq.len(), chunk.seq().len());
    let cl = chunk.cluster_num;
    let node = Node::new(chunk.id, aln.relstrand, seq, ops, position_from_start, cl);
    Some(node)
}

fn check_length(ops: &[Op], query_len: usize, chunk_len: usize) {
    let (mut query_length, mut chunk_length) = (0, 0);
    for op in ops.iter() {
        match op {
            Op::Match(x) => {
                query_length += x;
                chunk_length += x;
            }
            Op::Del(x) => chunk_length += x,
            Op::Ins(x) => query_length += x,
        }
    }
    assert_eq!(query_length, query_len);
    assert_eq!(chunk_length, chunk_len);
}

// Usually, the query is *longer* than the reference.
fn semiglobal(refr: &[u8], query: &[u8]) -> Vec<kiley::Op> {
    if refr.is_empty() {
        return vec![kiley::Op::Ins; query.len()];
    } else if query.is_empty() {
        return vec![kiley::Op::Del; refr.len()];
    }
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    // This is *reverse*. So we should fix it later.
    let aln = edlib_sys::align(refr, query, mode, task);
    let (start, end) = aln.location().unwrap();
    let end = end + 1;
    let leading = std::iter::repeat(kiley::Op::Ins).take(start);
    let aln = aln.operations().unwrap().iter();
    use kiley::Op;
    const OPS: [kiley::Op; 4] = [Op::Match, Op::Del, Op::Ins, Op::Mismatch];
    let aln = aln.map(|&op| OPS[op as usize]);
    let trailing = std::iter::repeat(kiley::Op::Ins).take(query.len() - end);
    leading.chain(aln).chain(trailing).collect()
}

pub fn remove_overlapping_encoding(mut nodes: Vec<Node>) -> Vec<Node> {
    fn start_end_iden(node: &definitions::Node) -> (usize, usize, f64) {
        let start = node.position_from_start;
        let (total, indel): (usize, usize) =
            node.cigar
                .0
                .iter()
                .fold((0, 0), |(total, indel), op| match op {
                    Op::Match(l) => (total + l, indel),
                    Op::Del(l) => (total + l, indel + l),
                    Op::Ins(l) => (total + l, indel + l),
                });
        let identity = 1f64 - indel as f64 / total as f64;
        (start, start + node.seq().len(), identity)
    }
    loop {
        let should_be_removed = nodes.windows(2).enumerate().find_map(|(i, w)| {
            let (former_start, former_end, former_identity) = start_end_iden(&w[0]);
            let (latter_start, latter_end, latter_identity) = start_end_iden(&w[1]);
            if (former_start <= latter_start && latter_end < former_end)
                || (latter_start <= former_start && former_end < latter_end)
            {
                if former_identity < latter_identity {
                    Some(i)
                } else {
                    Some(i + 1)
                }
            } else {
                None
            }
        });
        if let Some(idx) = should_be_removed {
            let _removed = nodes.remove(idx);
        } else {
            break;
        }
    }
    nodes
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
    // Remove overlapping same chunks.
    let mut deduped_nodes = vec![];
    let mut nodes = nodes.into_iter();
    let mut prev = nodes.next().unwrap();
    for node in nodes {
        let is_disjoint = prev.position_from_start + prev.query_length() < node.position_from_start;
        if prev.chunk != node.chunk || prev.is_forward != node.is_forward || is_disjoint {
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
        reference.push("chunks.fa");
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        for chunk in ds.selected_chunks.iter() {
            writeln!(wtr, ">{}\n{}", chunk.id, &chunk.seq)?;
        }
        let mut reads = c_dir.clone();
        reads.push("reads.fa");
        let mut wtr = std::fs::File::create(&reads).map(BufWriter::new)?;
        for read in ds.raw_reads.iter() {
            writeln!(wtr, ">{}\n{}", read.name, &read.seq)?;
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
        ReadType::CCS => args.extend(vec!["-H", "-k", "18"]),
        ReadType::CLR => args.extend(vec!["-H", "-k", "15"]),
        ReadType::ONT => args.extend(vec!["-k", "17"]),
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
        .all(|node| node.seq.iter().all(|c| c.is_ascii_uppercase()));
    let edges = read
        .edges
        .iter()
        .all(|edge| edge.label.iter().all(|c| c.is_ascii_uppercase()));
    nodes && edges
}
