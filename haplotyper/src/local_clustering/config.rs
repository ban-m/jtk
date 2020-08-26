use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
const SAMPLE_RATE: f64 = 0.01;
const INITIAL_BETA: f64 = 0.0005;
const BETA_INCREASE: f64 = 1.03;
const REPEAT_NUM: usize = 1;
const MAX_BETA: f64 = 0.8;
const GIBBS_PRIOR: f64 = 0.01;
const STABLE_LIMIT: u32 = 6;
const VARIANT_NUMBER: usize = 2;
const DEFAULT_SAMPLE_NUM: usize = 50;
#[derive(Debug, Clone)]
pub struct ClusteringConfig<F: Fn(u8, u8) -> i32> {
    pub cluster_num: usize,
    pub subchunk_length: usize,
    pub limit: u64,
    pub sample_rate: f64,
    pub alnparam: AlignmentParameters<F>,
    pub poa_config: poa_hmm::Config,
    pub seed: u64,
    pub id: u64,
    pub beta_increase: f64,
    pub initial_beta: f64,
    pub max_beta: f64,
    pub repeat_num: usize,
    pub gibbs_prior: f64,
    pub stable_limit: u32,
    pub variant_num: usize,
    pub read_type: ReadType,
    pub sample_num: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadType {
    CCS,
    CLR,
    ONT,
}

impl ClusteringConfig<fn(u8, u8) -> i32> {
    pub fn with_default(
        dataset: &definitions::DataSet,
        cluster_num: usize,
        subchunk_length: usize,
        limit: u64,
    ) -> Self {
        let seed = ((subchunk_length << cluster_num) / 17) as u64 * limit;
        let id: u64 = thread_rng().gen::<u64>() % 100_000;
        let bf = base_freq(&dataset.raw_reads);
        let units: HashMap<u64, &definitions::Unit> =
            dataset.selected_chunks.iter().map(|u| (u.id, u)).collect();
        let opss: Vec<_> = dataset
            .encoded_reads
            .iter()
            .flat_map(|rs| rs.nodes.iter())
            .filter_map(|node| Some(to_ops(node, units.get(&node.unit)?)))
            .collect();
        let config = summarize_operations(opss, bf);
        // let config = poa_hmm::PACBIO_CONFIG;
        debug!("Config:{}", config);
        Self {
            cluster_num,
            subchunk_length,
            limit,
            alnparam: DEFAULT_ALN,
            poa_config: config,
            sample_rate: SAMPLE_RATE,
            seed,
            id,
            beta_increase: BETA_INCREASE,
            initial_beta: INITIAL_BETA,
            max_beta: MAX_BETA,
            repeat_num: REPEAT_NUM,
            gibbs_prior: GIBBS_PRIOR,
            stable_limit: STABLE_LIMIT,
            variant_num: VARIANT_NUMBER,
            read_type: ReadType::CCS,
            sample_num: DEFAULT_SAMPLE_NUM,
        }
    }
    pub fn default() -> Self {
        let id: u64 = thread_rng().gen::<u64>() % 100_000;
        let config = poa_hmm::PACBIO_CONFIG;
        // debug!("Config:{}", config);
        Self {
            cluster_num: 3,
            subchunk_length: 100,
            limit: 2000,
            alnparam: DEFAULT_ALN,
            poa_config: config,
            sample_rate: SAMPLE_RATE,
            seed: 10000,
            id,
            beta_increase: BETA_INCREASE,
            initial_beta: INITIAL_BETA,
            max_beta: MAX_BETA,
            repeat_num: REPEAT_NUM,
            gibbs_prior: GIBBS_PRIOR,
            stable_limit: STABLE_LIMIT,
            variant_num: VARIANT_NUMBER,
            read_type: ReadType::CCS,
            sample_num: DEFAULT_SAMPLE_NUM,
        }
    }
    pub fn ccs(
        dataset: &definitions::DataSet,
        cluster_num: usize,
        subchunk_length: usize,
        limit: u64,
    ) -> Self {
        Self::with_default(dataset, cluster_num, subchunk_length, limit)
    }
    pub fn clr(
        dataset: &definitions::DataSet,
        cluster_num: usize,
        subchunk_length: usize,
        limit: u64,
    ) -> Self {
        let mut c = Self::with_default(dataset, cluster_num, subchunk_length, limit);
        c.initial_beta = 0.0001;
        c.max_beta = 0.3;
        c.beta_increase = 1.01;
        c.repeat_num = 1;
        c.stable_limit = 8;
        c.read_type = ReadType::CLR;
        c
    }
    pub fn ont(
        dataset: &definitions::DataSet,
        cluster_num: usize,
        subchunk_length: usize,
        limit: u64,
    ) -> Self {
        let mut c = Self::clr(dataset, cluster_num, subchunk_length, limit);
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

fn summarize_operations(opss: Vec<Vec<Op>>, base_freq: [f64; 4]) -> poa_hmm::Config {
    // match + mismatch.
    let mut matchmis = 0;
    let mut num_mis = 0;
    let mut num_seq = 0;
    let mut num_del = 0;
    let mut num_in = 0;
    let mut mm_after_mm = 0;
    let mut in_after_mm = 0;
    let mut in_after_del = 0;
    let mut in_after_in = 0;
    let mut del_after_mm = 0;
    let mut del_after_del = 0;
    // let mut del_after_in = 0;
    use Op::*;
    for ops in opss {
        num_seq += 1;
        matchmis += ops
            .iter()
            .filter(|&e| match e {
                Match | Mism => true,
                _ => false,
            })
            .count();
        num_mis += ops
            .iter()
            .filter(|&e| match e {
                Mism => true,
                _ => false,
            })
            .count();
        num_del += ops
            .iter()
            .filter(|&e| match e {
                Del => true,
                _ => false,
            })
            .count();
        num_in += ops
            .iter()
            .filter(|&e| match e {
                In => true,
                _ => false,
            })
            .count();
        for before_after in ops.windows(2) {
            let b = before_after[0];
            let a = before_after[1];
            match (b, a) {
                (Match, Mism) | (Match, Match) | (Mism, Match) | (Mism, Mism) => mm_after_mm += 1,
                (Mism, Del) | (Match, Del) => del_after_mm += 1,
                (Del, Del) => del_after_del += 1,
                //(In, Del) => del_after_in += 1,
                (Mism, In) | (Match, In) => in_after_mm += 1,
                (In, In) => in_after_in += 1,
                (Del, In) => in_after_del += 1,
                _ => {}
            }
        }
    }
    let div = |x, y| x as f64 / y as f64;
    let p_mismatch = div(num_mis, matchmis);
    matchmis -= num_seq;
    //eprintln!("{}\t{}\t{}", matchmis, num_in,num_del);
    let p_match = div(mm_after_mm, matchmis);
    let p_start_in = div(in_after_mm, matchmis);
    let p_start_del = div(del_after_mm, matchmis);
    let p_ext_in = div(in_after_in, num_in);
    let p_ext_del = div(del_after_del, num_del);
    let p_del_to_in = div(in_after_del, num_del);
    poa_hmm::Config {
        mismatch: p_mismatch,
        base_freq,
        p_match,
        p_ins: p_start_in,
        p_del: p_start_del,
        p_extend_ins: p_ext_in,
        p_extend_del: p_ext_del,
        p_del_to_ins: p_del_to_in,
    }
}
