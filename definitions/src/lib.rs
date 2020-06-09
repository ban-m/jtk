//! Definitions -- A tiny interface for HLA-typing problem.
//! Roughly speaking, we incorporate with other programs, pass messages, or interact with other CLI via JSON object format. Specifically, the message is encoded only one, possibly large, structure named [DataSet](DataSet)

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DataSet {
    pub raw_reads: Vec<RawRead>,
    pub hic_pairs: Vec<HiCPair>,
    pub selected_chunks: Vec<Unit>,
    pub encoded_reads: Vec<EncodedRead>,
    pub hic_edges: Vec<HiCEdge>,
}

impl DataSet {
    pub fn with_param(
        raw_reads: Vec<RawRead>,
        hic_pairs: Vec<HiCPair>,
        selected_chunks: Vec<Unit>,
        encoded_reads: Vec<EncodedRead>,
        hic_edges: Vec<HiCEdge>,
    ) -> Self {
        Self {
            raw_reads,
            hic_pairs,
            selected_chunks,
            encoded_reads,
            hic_edges,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawRead {
    pub name: String,
    pub desc: String,
    pub id: u64,
    pub seq: String,
}

impl RawRead {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HiCPair {
    pub pair1: u64,
    pub pair2: u64,
    pub pair_id: u64,
    pub seq1: String,
    pub seq2: String,
}

impl HiCPair {
    pub fn seq1(&self) -> &[u8] {
        self.seq1.as_bytes()
    }
    pub fn seq2(&self) -> &[u8] {
        self.seq2.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Unit {
    pub id: u64,
    pub seq: String,
    pub cluster: u32,
}
impl Unit {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EncodedRead {
    pub original_length: usize,
    pub leading_gap: usize,
    pub trailing_gap: usize,
    pub id: u64,
    pub edges: Vec<Edge>,
    pub nodes: Vec<Node>,
}

impl EncodedRead {
    pub fn is_gappy(&self) -> bool {
        self.nodes.is_empty()
    }
    pub fn encoded_rate(&self) -> f64 {
        let encoded_length = self.encoded_length();
        encoded_length as f64 / self.original_length as f64
    }
    pub fn encoded_length(&self) -> usize {
        let sum = self.nodes.iter().map(|n| n.query_length()).sum::<usize>();
        let offset = self
            .edges
            .iter()
            .map(|e| e.offset)
            .filter(|&e| e < 0)
            .sum::<i64>();
        let length = sum as i64 + offset;
        if length < 0 {
            0
        } else {
            length as usize
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge {
    pub from: u64,
    pub to: u64,
    pub offset: i64,
    pub label: String,
}

impl Edge {
    pub fn label(&self) -> &[u8] {
        self.label.as_bytes()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Node {
    pub unit: u64,
    pub seq: String,
    pub is_forward: bool,
    pub cigar: Vec<Op>,
}
impl Node {
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
    pub fn query_length(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                Op::MisMatch(l) | Op::Match(l) | Op::Ins(l) => *l,
                Op::Del(_) => 0,
            })
            .sum::<usize>()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Copy)]
pub enum Op {
    Match(usize),
    MisMatch(usize),
    /// Deletion with respect to the reference.
    Del(usize),
    /// Insertion with respect to the reference.
    Ins(usize),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HiCEdge {
    pub pair_id: u64,
    pub pair1: u64,
    pub pair2: u64,
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
