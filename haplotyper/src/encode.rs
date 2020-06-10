use bio_utils::sam::Sam;
use std::collections::HashMap;
pub trait Encode {
    fn encode(self, alignmnets: &[Sam]) -> Self;
}

impl Encode for definitions::DataSet {
    fn encode(mut self, alignments: &[Sam]) -> Self {
        let alignments_each_reads: HashMap<String, Vec<&Sam>> = distribute(alignments);
        let encoded_reads: Vec<_> = self
            .raw_reads
            .iter()
            .filter_map(|read| {
                let alns = alignments_each_reads.get(&read.name)?;
                encode(read, alns, &self.selected_chunks)
            })
            .collect();
        self.encoded_reads = encoded_reads;
        self
    }
}

fn distribute<'a>(alignments: &'a [Sam]) -> HashMap<String, Vec<&'a Sam>> {
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for alignment in alignments {
        let q_name = alignment.q_name().to_string();
        buckets.entry(q_name).or_default().push(alignment);
    }
    buckets
}

use definitions::{Edge, EncodedRead, Node, Op, RawRead, Unit};
fn encode(read: &RawRead, alignments: &[&Sam], units: &[Unit]) -> Option<EncodedRead> {
    let mut nodes: Vec<_> = alignments
        .iter()
        .filter_map(|aln| encode_alignment(aln, units, read))
        .collect();
    nodes.sort_by_key(|e| e.position_from_start);
    let edges: Vec<_> = nodes.windows(2).map(Edge::from_nodes).collect();
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

fn encode_alignment(aln: &Sam, units: &[Unit], read: &RawRead) -> Option<Node> {
    None
}

fn consumed_reference_length(cigar: &[Op]) -> usize {
    cigar
        .iter()
        .map(|op| match op {
            Op::Match(l) | Op::MisMatch(l) | Op::Del(l) => *l,
            Op::Ins(l) => 0,
        })
        .sum::<usize>()
}
