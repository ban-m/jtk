use bio_utils::lasttab;
use bio_utils::lasttab::LastTAB;
use std::collections::HashMap;
pub trait Encode {
    fn encode(self, alignmnets: &[LastTAB]) -> Self;
}

impl Encode for definitions::DataSet {
    fn encode(mut self, alignments: &[LastTAB]) -> Self {
        let alignments_each_reads: HashMap<String, Vec<&LastTAB>> = distribute(alignments);
        let encoded_reads: Vec<_> = self
            .raw_reads
            .iter()
            .filter_map(|read| {
                let alns = alignments_each_reads.get(&read.name)?;
                encode(read, alns, &self.selected_chunks)
            })
            .collect();
        debug!("Encoding {} reads.", encoded_reads.len());
        self.encoded_reads = encoded_reads;
        self
    }
}

fn distribute<'a>(alignments: &'a [LastTAB]) -> HashMap<String, Vec<&'a LastTAB>> {
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for alignment in alignments {
        let q_name = alignment.seq2_name().to_string();
        buckets.entry(q_name).or_default().push(alignment);
    }
    buckets
}

use definitions::{Edge, EncodedRead, Node, Op, RawRead, Unit};
fn encode(read: &RawRead, alignments: &[&LastTAB], units: &[Unit]) -> Option<EncodedRead> {
    let mut nodes: Vec<_> = alignments
        .iter()
        .filter(|aln| aln.seq1_matchlen() > aln.seq1_len() * 98 / 100)
        .filter_map(|aln| encode_alignment(aln, units, read))
        .collect();
    nodes.sort_by_key(|e| e.position_from_start);
    let edges: Vec<_> = nodes
        .windows(2)
        .map(|w| Edge::from_nodes(w, read.seq()))
        .collect();
    let leading_gap = nodes.first()?.position_from_start;
    let trailing_gap = {
        let last = nodes.last()?;
        read.seq().len() - last.position_from_start + consumed_reference_length(&last.cigar)
    };
    Some(EncodedRead {
        original_length: read.seq().len(),
        id: read.id,
        edges,
        nodes,
        leading_gap,
        trailing_gap,
    })
}

fn encode_alignment(aln: &LastTAB, _units: &[Unit], read: &RawRead) -> Option<Node> {
    let position = aln.seq2_start_from_forward();
    let unit: u64 = aln.seq1_name().parse().ok()?;
    let is_forward = aln.seq2_direction().is_forward();
    let seq = {
        let start = aln.seq2_start();
        let end = start + aln.seq2_matchlen();
        let seq = if is_forward {
            read.seq().to_vec()
        } else {
            bio_utils::revcmp(read.seq())
        };
        seq[start..end].to_vec()
    };
    let cigar = convert_aln_to_cigar(aln);
    Some(Node {
        position_from_start: position,
        unit,
        cluster: 0,
        seq: String::from_utf8_lossy(&seq).to_string(),
        is_forward,
        cigar,
    })
}

fn convert_aln_to_cigar(aln: &lasttab::LastTAB) -> Vec<Op> {
    let mut cigar = if aln.seq1_start_from_forward() != 0 {
        vec![Op::Del(aln.seq1_start_from_forward())]
    } else {
        vec![]
    };
    cigar.extend(aln.alignment().into_iter().map(|op| match op {
        lasttab::Op::Seq1In(l) => Op::Ins(l),
        lasttab::Op::Seq2In(l) => Op::Del(l),
        lasttab::Op::Match(l) => Op::Match(l),
    }));
    let reflen = consumed_reference_length(&cigar);
    assert!(reflen <= aln.seq1_len(), "{} > {}", reflen, aln.seq1_len());
    if aln.seq1_len() > reflen {
        cigar.push(Op::Del(aln.seq1_len() - reflen))
    }
    cigar
}

fn consumed_reference_length(cigar: &[Op]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            Op::Match(l) | Op::Del(l) => *l,
            Op::Ins(_) => 0,
        })
        .sum::<usize>()
}
