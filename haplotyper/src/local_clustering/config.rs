#![allow(dead_code)]
use definitions::ReadType;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
// use std::collections::HashMap;
const STABLE_LIMIT: u32 = 6;
const VARIANT_NUMBER: usize = 2;
const P_VALUE: f64 = 0.01;
const RETRY_LIMIT: u64 = 4;
#[derive(Debug, Clone)]
pub struct ClusteringConfig<F: Fn(u8, u8) -> i32> {
    pub cluster_num: usize,
    pub subchunk_length: usize,
    pub limit: u64,
    pub alnparam: AlignmentParameters<F>,
    pub id: u64,
    pub stable_limit: u32,
    pub variant_num: usize,
    pub read_type: ReadType,
    pub p_value: f64,
    pub retry_limit: u64,
    pub retain_current_clustering: bool,
}

impl ClusteringConfig<fn(u8, u8) -> i32> {
    pub fn with_default(
        dataset: &definitions::DataSet,
        cluster_num: usize,
        subchunk_length: usize,
    ) -> Self {
        let id: u64 = thread_rng().gen::<u64>() % 100_000;
        Self {
            cluster_num,
            subchunk_length,
            limit: 600,
            alnparam: DEFAULT_ALN,
            id,
            stable_limit: STABLE_LIMIT,
            variant_num: VARIANT_NUMBER,
            read_type: dataset.read_type,
            p_value: P_VALUE,
            retry_limit: 10,
            retain_current_clustering: false,
        }
    }
    pub fn default() -> Self {
        let id: u64 = thread_rng().gen::<u64>() % 100_000;
        Self {
            cluster_num: 3,
            subchunk_length: 100,
            limit: 2000,
            alnparam: DEFAULT_ALN,
            id,
            stable_limit: STABLE_LIMIT,
            variant_num: VARIANT_NUMBER,
            read_type: ReadType::CCS,
            p_value: P_VALUE,
            retry_limit: RETRY_LIMIT,
            retain_current_clustering: false,
        }
    }
    pub fn ccs(dataset: &definitions::DataSet, cluster_num: usize, subchunk_length: usize) -> Self {
        Self::with_default(dataset, cluster_num, subchunk_length)
    }
    pub fn clr(dataset: &definitions::DataSet, cluster_num: usize, subchunk_length: usize) -> Self {
        let mut c = Self::with_default(dataset, cluster_num, subchunk_length);
        c.read_type = ReadType::CLR;
        c
    }
    pub fn ont(dataset: &definitions::DataSet, cluster_num: usize, subchunk_length: usize) -> Self {
        let mut c = Self::clr(dataset, cluster_num, subchunk_length);
        c.read_type = ReadType::ONT;
        c
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentParameters<F>
where
    F: Fn(u8, u8) -> i32,
{
    pub ins: i32,
    pub del: i32,
    pub score: F,
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        3
    } else {
        -7
    }
}

pub const DEFAULT_ALN: AlignmentParameters<fn(u8, u8) -> i32> = AlignmentParameters {
    ins: -4,
    del: -4,
    score,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Op {
    Match,
    Mism,
    Del,
    In,
}

fn to_ops(node: &definitions::Node, unit: &definitions::Unit) -> Vec<Op> {
    assert!(node.unit == unit.id);
    let (mut r, mut q) = (0, 0);
    let refr = unit.seq();
    let query = node.seq();
    node.cigar
        .iter()
        .flat_map(|op| match op {
            definitions::Op::Match(l) => {
                let ops: Vec<_> = refr[r..r + l]
                    .iter()
                    .zip(query[q..q + l].iter())
                    .map(|(a, b)| if a == b { Op::Match } else { Op::Mism })
                    .collect();
                r += l;
                q += l;
                ops
            }
            definitions::Op::Ins(l) => {
                q += l;
                vec![Op::In; *l]
            }
            definitions::Op::Del(l) => {
                r += l;
                vec![Op::Del; *l]
            }
        })
        .collect()
}

fn base_freq(rs: &[definitions::RawRead]) -> [f64; 4] {
    let tot = rs.iter().map(|read| read.seq().len()).sum::<usize>() as f64;
    let mut base_count = [0.; 4];
    for base in rs.iter().flat_map(|e| e.seq().iter()) {
        match base {
            b'A' => base_count[0] += 1.,
            b'C' => base_count[1] += 1.,
            b'G' => base_count[2] += 1.,
            b'T' => base_count[3] += 1.,
            _ => {}
        }
    }
    base_count.iter_mut().for_each(|e| *e /= tot);
    base_count
}
