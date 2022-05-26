//! Definitions -- A tiny interface for HLA-typing problem.
//! Roughly speaking, we incorporate with other programs, pass messages, or interact with other CLI via JSON object format. Specifically, the message is encoded only one, possibly large, structure named [DataSet](DataSet)

use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DataSet {
    /// The path to the input file.
    pub input_file: String,
    /// Masking information
    pub masked_kmers: MaskInfo,
    /// If Some(x), it is the estimated coverage per haplotype.
    /// Note that this value is a estimation, so please do not relay heavily on this value.
    pub coverage: Option<f64>,
    /// The FASTA-like record for input. If this value is empy, please reload the original
    /// sequence from either `input_file` or `encoded_reads`.
    pub raw_reads: Vec<RawRead>,
    /// The HiC reads.
    pub hic_pairs: Vec<HiCPair>,
    /// The chunks selected.
    pub selected_chunks: Vec<Unit>,
    /// The reads encoded by selected chunks.
    pub encoded_reads: Vec<EncodedRead>,
    /// The edge of HiC.
    pub hic_edges: Vec<HiCEdge>,
    /// Depricated: This value is previously assigned by `global_clustering` method.
    /// As the research has proceed, we realize that `phasing` reads is unnessesary, or even harmful to
    /// the assembly. So, in the future refactoring, this value and `global_clustering` method would be
    /// removed.
    pub assignments: Vec<Assignment>,
    /// The type of the reads.
    pub read_type: ReadType,
    // /// Estimated parameter for a statistical model (used to phase local region)
    // pub model_params: Option<ModelParameters>,
    /// Estimated error rate.
    pub error_rate: ErrorRate,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MaskInfo {
    pub k: usize,
    /// Occurent threshold(k-mers occurring more than this value would be masked)
    pub thr: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Copy)]
pub enum ReadType {
    CCS,
    CLR,
    ONT,
    None,
}

pub const CLR_BAND_WIDTH: usize = 100;
pub const HIFI_BAND_WIDTH: usize = 80;
pub const ONT_BAND_WIDTH: usize = 100;

pub const CLR_CTG_SIM: f64 = 0.15;
pub const CLR_CLR_SIM: f64 = 0.20;
pub const HIFI_SIM_THR: f64 = 0.05;
pub const ONT_SIM_THR: f64 = 0.15;

pub const CLR_BAND_FRAC: f64 = 0.03;
pub const ONT_BAND_FRAC: f64 = 0.03;
pub const HIFI_BAND_FRAC: f64 = 0.01;

impl ReadType {
    pub fn sim_thr(&self) -> f64 {
        match *self {
            ReadType::CCS => HIFI_SIM_THR,
            ReadType::None | ReadType::CLR => CLR_CLR_SIM,
            ReadType::ONT => ONT_SIM_THR,
        }
    }
    pub fn band_width(&self, len: usize) -> usize {
        let len = len as f64
            * match *self {
                ReadType::CCS => HIFI_BAND_FRAC,
                ReadType::CLR => CLR_BAND_FRAC,
                ReadType::ONT => ONT_BAND_FRAC,
                ReadType::None => CLR_BAND_FRAC,
            };
        len.ceil() as usize
    }
    pub fn min_span_reads(&self) -> usize {
        match *self {
            ReadType::CCS => 1,
            ReadType::CLR => 3,
            ReadType::ONT => 2,
            ReadType::None => 3,
        }
    }
    pub fn min_llr_value(&self) -> f64 {
        match *self {
            ReadType::CCS => 0.1f64,
            ReadType::CLR => 1f64,
            ReadType::ONT => 1f64,
            ReadType::None => 1f64,
        }
    }
}

impl std::default::Default for DataSet {
    fn default() -> Self {
        Self {
            input_file: String::new(),
            coverage: None,
            raw_reads: vec![],
            hic_pairs: vec![],
            selected_chunks: vec![],
            encoded_reads: vec![],
            hic_edges: vec![],
            assignments: vec![],
            read_type: ReadType::None,
            masked_kmers: MaskInfo::default(),
            // model_params: None,
            error_rate: ErrorRate::default(),
        }
    }
}

impl DataSet {
    /// Return an empty dataset.
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_minimum_data(
        input_file: &str,
        raw_reads: Vec<RawRead>,
        read_type: ReadType,
    ) -> Self {
        let assignments: Vec<_> = raw_reads
            .iter()
            .map(|r| Assignment {
                id: r.id,
                cluster: 0,
            })
            .collect();
        Self {
            input_file: input_file.to_string(),
            coverage: None,
            raw_reads,
            hic_pairs: vec![],
            selected_chunks: vec![],
            encoded_reads: vec![],
            hic_edges: vec![],
            assignments,
            read_type,
            masked_kmers: MaskInfo::default(),
            // model_params: None,
            error_rate: ErrorRate::guess(read_type),
        }
    }
    /// Sanity check function. Call it to ensure that some properties indeed holds.
    /// Currently, the following properties are checked.
    /// 1: The input file exists.
    /// 2: Every encoded read has its original read.
    /// 3: Every encoded read can be correctly recovered into the original sequence.
    /// Of course, these should be hold at any steps in the pipeline,
    /// So, this is just a checking function.
    pub fn sanity_check(&self) {
        let units: HashSet<_> = self.selected_chunks.iter().map(|c| c.id).collect();
        self.encoded_reads
            .iter()
            .flat_map(|r| r.nodes.iter())
            .for_each(|n| assert!(units.contains(&n.unit)));
        // assert!(std::path::Path::new(&self.input_file).exists());
        self.encoded_reads_can_be_recovered();
        use std::collections::HashSet;
        let mut units = HashSet::new();
        for unit in self.selected_chunks.iter() {
            assert!(!units.contains(&unit.id));
            units.insert(unit.id);
        }
        for chunk in self.selected_chunks.iter() {
            assert!(chunk.cluster_num <= chunk.copy_num);
        }
        use std::collections::HashMap;
        let max_cl_num: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        for node in self.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            assert!(node.cluster <= max_cl_num[&node.unit] as u64);
        }
    }
    fn encoded_reads_can_be_recovered(&self) {
        use std::collections::HashMap;
        let seq: HashMap<_, _> = self.raw_reads.iter().map(|x| (x.id, x.seq())).collect();
        self.encoded_reads.iter().for_each(|read| {
            let orig: Vec<_> = match seq.get(&read.id) {
                Some(res) => res.iter().map(u8::to_ascii_uppercase).collect(),
                None => panic!(),
            };
            let recover: Vec<_> = read
                .recover_raw_read()
                .iter()
                .map(u8::to_ascii_uppercase)
                .collect();
            assert_eq!(orig.len(), recover.len(), "{}", read.nodes.len());
            let orig_len = read.original_length;
            assert_eq!(orig.len(), orig_len);
            let not_matched = orig
                .iter()
                .zip(recover.iter())
                .enumerate()
                .filter(|(_, (x, y))| x != y);
            for (idx, _) in not_matched {
                eprintln!("{}", idx);
            }
            if orig != recover {
                eprintln!("{}", read.leading_gap.len());
                eprintln!("{}", read.trailing_gap.len());
                panic!();
            }
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawRead {
    /// Name of the read. It is the `id` in the fasta/fastq file.
    pub name: String,
    pub desc: String,
    /// The id of the read. It is automatically given by jtk program.
    pub id: u64,
    /// Sequence. It is a string on an alphabet of A,C,G,T,a,c,g,t.
    /// (i.e., lowercase included)
    pub seq: DNASeq,
}

impl RawRead {
    pub fn seq(&self) -> &[u8] {
        &self.seq.0
    }
}

impl std::fmt::Display for RawRead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} {} {}\n{}", self.name, self.desc, self.id, self.seq)
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
    /// Unit sequence.
    pub seq: DNASeq,
    /// Current estimation of the cluster number.
    pub cluster_num: usize,
    /// Current estimation of the copy number. This is an upper bound of the cluster number.
    /// (For me: cluster_num is not always the copy number. There can be a homologous chunk.)
    pub copy_num: usize,
    /// Local clustering score. If not clustered, zero.
    pub score: f64,
}

use serde_with::DeserializeFromStr;
use serde_with::SerializeDisplay;
#[derive(Debug, Clone, SerializeDisplay, DeserializeFromStr, Default)]
pub struct DNASeq(Vec<u8>);

impl std::fmt::Display for DNASeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", std::str::from_utf8(&self.0).unwrap())
    }
}

impl DNASeq {
    pub fn pop(&mut self) -> Option<u8> {
        self.0.pop()
    }
    pub fn extend<I: std::iter::IntoIterator<Item = u8>>(&mut self, iter: I) {
        self.0.extend(iter);
    }
    pub fn iter(&self) -> std::slice::Iter<'_, u8> {
        self.0.iter()
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn as_mut(&mut self) -> &mut [u8] {
        self.0.as_mut()
    }
    pub fn as_slice(&self) -> &[u8] {
        self.0.as_slice()
    }
}

impl std::str::FromStr for DNASeq {
    type Err = u64;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(DNASeq(s.as_bytes().to_vec()))
    }
}

impl std::convert::From<Vec<u8>> for DNASeq {
    fn from(seq: Vec<u8>) -> Self {
        Self(seq)
    }
}

impl std::convert::Into<Vec<u8>> for DNASeq {
    fn into(self) -> Vec<u8> {
        self.0
    }
}

impl Unit {
    pub fn new(id: u64, seq: Vec<u8>, copy_num: usize) -> Self {
        Self {
            id,
            seq: seq.into(),
            copy_num,
            cluster_num: 1,
            score: 0f64,
        }
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq.0
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct EncodedRead {
    pub id: u64,
    pub original_length: usize,
    pub leading_gap: DNASeq,
    pub trailing_gap: DNASeq,
    pub edges: Vec<Edge>,
    pub nodes: Vec<Node>,
}

impl std::fmt::Display for EncodedRead {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "{}({}bp)", self.id, self.original_length)?;
        write!(f, "{}bp gap|", self.leading_gap.len())?;
        for e in self.edges.iter() {
            write!(f, "{} ", e)?;
        }
        write!(f, "|{} bp gap", self.trailing_gap.len())
    }
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
    /// Remove the i-th node.
    pub fn remove(&mut self, i: usize) {
        assert!(i < self.nodes.len());
        assert_eq!(self.nodes.len(), self.edges.len() + 1);
        let len = self.nodes.len();
        let removed_node = self.nodes.remove(i);
        if self.nodes.is_empty() {
            assert!(self.edges.is_empty());
            self.leading_gap.extend(removed_node.original_seq());
            return;
        }
        if i + 1 == len {
            let removed_edge = self.edges.remove(i - 1);
            let node_seq = removed_node.original_seq();
            let skip = match removed_edge.offset {
                x if x < 0 => -x as usize,
                _ => 0,
            };
            let trailing_seq: Vec<_> = removed_edge
                .label()
                .iter()
                .chain(node_seq.iter())
                .chain(self.trailing_gap.iter())
                .skip(skip)
                .copied()
                .collect();
            self.trailing_gap = trailing_seq.into();
        } else if i == 0 {
            let removed_edge = self.edges.remove(i);
            self.leading_gap.extend(removed_node.original_seq());
            self.leading_gap
                .extend(removed_edge.label().iter().copied());
            for _ in 0..(-removed_edge.offset) {
                self.leading_gap.pop();
            }
        } else {
            let removed_edge = self.edges.remove(i);
            // First add the removed sequence and edges into the current edges.
            let mut edge: Vec<_> = self.edges[i - 1]
                .label
                .iter()
                .chain(removed_node.original_seq().iter())
                .chain(removed_edge.label.iter())
                .copied()
                .collect();
            // Then, fix the offset by the offset.
            if self.edges[i - 1].offset < 0 {
                edge.reverse();
                for _ in 0..-self.edges[i - 1].offset {
                    edge.pop();
                }
                edge.reverse();
            }
            if removed_edge.offset < 0 {
                for _ in 0..-removed_edge.offset {
                    edge.pop();
                }
            }
            self.edges[i - 1].to = removed_edge.to;
            //self.edges[i - 1].label = edge.iter().map(|&x| x as char).collect();
            self.edges[i - 1].label = edge.to_vec().into();
            self.edges[i - 1].offset += removed_node.seq().len() as i64 + removed_edge.offset;
        }
        assert_eq!(self.nodes.len(), self.edges.len() + 1);
    }
    pub fn recover_raw_read(&self) -> Vec<u8> {
        let mut original_seq: Vec<_> = self.leading_gap.clone().into();
        for (n, e) in self.nodes.iter().zip(self.edges.iter()) {
            assert_eq!(n.unit, e.from);
            original_seq.extend(n.original_seq());
            for _ in 0..(-e.offset).max(0) {
                original_seq.pop();
            }
            original_seq.extend(e.label());
        }
        if let Some(n) = self.nodes.last() {
            original_seq.extend(n.original_seq());
        }
        original_seq.extend(self.trailing_gap.iter());
        original_seq
    }
    /// Return true if this read contains (unit,cluster)-node. Linear time.
    pub fn contains(&self, (unit, cluster): (u64, u64)) -> bool {
        self.nodes
            .iter()
            .any(|n| n.unit == unit && n.cluster == cluster)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Edge {
    pub from: u64,
    pub to: u64,
    pub offset: i64,
    /// This is a string on an alphabet of A,C,G,T. There might be some lowercase characters.
    pub label: DNASeq,
}

impl std::fmt::Display for Edge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}({}){}", self.from, self.offset, self.to)
    }
}

impl Edge {
    pub fn label(&self) -> &[u8] {
        &self.label.0
    }
    pub fn from_nodes(ns: &[Node], seq: &[u8]) -> Self {
        let (from, to) = match *ns {
            [ref from, ref to] => (from, to),
            _ => unreachable!(),
        };
        let end = from.position_from_start + from.query_length();
        let start = to.position_from_start;
        let label = if start <= end {
            Vec::new()
        } else {
            seq.iter()
                .take(start)
                .skip(end)
                .map(u8::to_ascii_uppercase)
                .collect()
        };
        Edge {
            from: from.unit,
            to: to.unit,
            offset: start as i64 - end as i64,
            label: label.into(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Node {
    /// 0-index.
    pub position_from_start: usize,
    pub unit: u64,
    pub cluster: u64,
    /// Sequence. No lowercase included. Already rev-comped.
    pub seq: DNASeq,
    pub is_forward: bool,
    pub cigar: Ops,
    // TODO: This member should have a different name, such as
    // `likelihood gain` or something.
    pub posterior: Vec<f64>,
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}-{}({}bp,{},{})",
            self.unit,
            self.cluster,
            self.seq.len(),
            self.is_forward as u8,
            self.position_from_start,
        )
    }
}

impl Node {
    /// Create a new instance of node.
    /// As usually it is created before the local clustering,
    /// It only accepts the minimum requirements, i.e, the position in the read, the units, the direction, the sequence, and the cigar string.
    pub fn new(
        unit: u64,
        is_forward: bool,
        seq: &[u8],
        cigar: Vec<Op>,
        position_from_start: usize,
        cluster_num: usize,
    ) -> Self {
        // This is cluster number. Because each node would be assigned to a cluster, not a copy.
        // Actually, we can not determine each copy of in the genome, right?
        let post_prob = (cluster_num.max(1) as f64).recip().ln();
        Self {
            position_from_start,
            unit,
            cluster: 0,
            seq: seq.to_vec().into(),
            is_forward,
            cigar: cigar.into(),
            posterior: vec![post_prob; cluster_num],
        }
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq.0
    }
    pub fn original_seq(&self) -> Vec<u8> {
        if self.is_forward {
            self.seq.clone().into()
        } else {
            self.seq
                .iter()
                .rev()
                .map(|x| match x.to_ascii_uppercase() {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'T' => b'A',
                    _ => panic!(),
                })
                .collect()
        }
    }
    pub fn query_length(&self) -> usize {
        self.cigar
            .iter()
            .map(|op| match op {
                Op::Match(l) | Op::Ins(l) => *l,
                Op::Del(_) => 0,
            })
            .sum::<usize>()
    }
    /// Return (node path, alignment, unit path)
    pub fn recover(&self, unit: &Unit) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (read, unit) = (self.seq(), unit.seq());
        let (mut q, mut al, mut r) = (vec![], vec![], vec![]);
        let (mut q_pos, mut r_pos) = (0, 0);
        let match_char = |(x, y): (&u8, &u8)| {
            if x.to_ascii_uppercase() == y.to_ascii_uppercase() {
                b'|'
            } else {
                b'X'
            }
        };
        for op in self.cigar.iter() {
            match *op {
                Op::Match(l) => {
                    al.extend(
                        read[q_pos..q_pos + l]
                            .iter()
                            .zip(&unit[r_pos..r_pos + l])
                            .map(match_char),
                    );
                    q.extend(read[q_pos..q_pos + l].iter().copied());
                    r.extend(unit[r_pos..r_pos + l].iter().copied());
                    q_pos += l;
                    r_pos += l;
                }
                Op::Del(l) => {
                    al.extend(vec![b' '; l]);
                    q.extend(vec![b' '; l]);
                    r.extend(unit[r_pos..r_pos + l].iter().copied());
                    r_pos += l;
                }
                Op::Ins(l) => {
                    al.extend(vec![b' '; l]);
                    q.extend(read[q_pos..q_pos + l].iter().copied());
                    r.extend(vec![b' '; l]);
                    q_pos += l;
                }
            }
        }
        (q, al, r)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Copy, Eq, PartialEq)]
pub enum Op {
    Match(usize),
    /// Deletion with respect to the reference.
    Del(usize),
    /// Insertion with respect to the reference.
    Ins(usize),
}

#[derive(Debug, Clone, SerializeDisplay, DeserializeFromStr, Default)]
pub struct Ops(Vec<Op>);
impl std::fmt::Display for Ops {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for &x in self.0.iter() {
            match x {
                Op::Match(x) => write!(f, "{x}M")?,
                Op::Del(x) => write!(f, "{x}D")?,
                Op::Ins(x) => write!(f, "{x}I")?,
            }
        }
        Ok(())
    }
}

impl Ops {
    pub fn iter(&self) -> std::slice::Iter<'_, Op> {
        self.0.iter()
    }
}
impl std::convert::From<Vec<Op>> for Ops {
    fn from(xs: Vec<Op>) -> Self {
        Ops(xs)
    }
}
impl std::convert::Into<Vec<Op>> for Ops {
    fn into(self) -> Vec<Op> {
        self.0
    }
}

impl std::str::FromStr for Ops {
    type Err = u64;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ops = vec![];
        let mut num = 0;
        for x in s.bytes() {
            if x.is_ascii_digit() {
                num = 10 * num + (x - b'0') as usize;
            } else {
                match x {
                    b'M' => ops.push(Op::Match(num)),
                    b'D' => ops.push(Op::Del(num)),
                    b'I' => ops.push(Op::Ins(num)),
                    _ => panic!("{}", x),
                }
                num = 0;
            }
        }
        Ok(Ops(ops))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HiCEdge {
    pub pair_id: u64,
    pub pair1: u64,
    pub pair2: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Assignment {
    pub id: u64,
    pub cluster: usize,
}

impl Assignment {
    pub fn new(id: u64, cluster: usize) -> Self {
        Self { id, cluster }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ModelParameters {
    pub del_open: f64,
    pub del_ext: f64,
    pub ins_open: f64,
    pub ins_ext: f64,
    pub match_prob: f64,
}

/// The error rate with respect to the length of the alignment.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ErrorRate {
    pub del: f64,
    pub del_sd: f64,
    pub ins: f64,
    pub ins_sd: f64,
    pub mismatch: f64,
    pub mism_sd: f64,
    pub total: f64,
    pub total_sd: f64,
}

impl std::fmt::Display for ErrorRate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Del:{:.3}/{:.3}", self.del, self.del_sd)?;
        writeln!(f, "Ins:{:.3}/{:.3}", self.ins, self.ins_sd)?;
        writeln!(f, "Mism:{:.3}/{:.3}", self.mismatch, self.mism_sd)?;
        write!(f, "Total:{:.3}/{:.3}", self.total, self.total_sd)
    }
}

const CCS_ERROR_RATE: ErrorRate = ErrorRate {
    del: 0.005,
    del_sd: 0.001,
    ins: 0.005,
    ins_sd: 0.001,
    mismatch: 0.005,
    mism_sd: 0.001,
    total: 0.01,
    total_sd: 0.005,
};

const CLR_ERROR_RATE: ErrorRate = ErrorRate {
    del: 0.07,
    del_sd: 0.02,
    ins: 0.06,
    ins_sd: 0.02,
    mismatch: 0.02,
    mism_sd: 0.01,
    total: 0.15,
    total_sd: 0.03,
};

const ONT_ERROR_RATE: ErrorRate = ErrorRate {
    del: 0.01,
    del_sd: 0.005,
    ins: 0.01,
    ins_sd: 0.005,
    mismatch: 0.01,
    mism_sd: 0.005,
    total: 0.03,
    total_sd: 0.008,
};

impl ErrorRate {
    pub fn new(
        (del, del_sd): (f64, f64),
        (ins, ins_sd): (f64, f64),
        (mism, mism_sd): (f64, f64),
        (total, total_sd): (f64, f64),
    ) -> Self {
        Self {
            del,
            del_sd,
            ins,
            ins_sd,
            mismatch: mism,
            mism_sd,
            total,
            total_sd,
        }
    }
    pub fn guess(readtype: ReadType) -> Self {
        match readtype {
            ReadType::CCS => CCS_ERROR_RATE,
            ReadType::CLR => CLR_ERROR_RATE,
            ReadType::ONT => ONT_ERROR_RATE,
            ReadType::None => CLR_ERROR_RATE,
        }
    }
    pub fn sum(&self) -> f64 {
        self.del + self.ins + self.mismatch
    }
}

impl std::ops::Add for ErrorRate {
    type Output = ErrorRate;
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            del: self.del + rhs.del,
            del_sd: self.del_sd + rhs.del_sd,
            ins: self.ins + rhs.ins,
            ins_sd: self.ins_sd + rhs.ins_sd,
            mismatch: self.mismatch + rhs.mismatch,
            mism_sd: self.mism_sd + rhs.mism_sd,
            total: self.total + rhs.total,
            total_sd: self.total_sd + rhs.total_sd,
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
