use definitions::ReadType;
use definitions::*;
use poa_hmm::POA;
use rand::{seq::SliceRandom, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct PolishUnitConfig {
    read_type: ReadType,
    consensus_size: usize,
    filter_size: usize,
    rep_num: usize,
    seed: u64,
}

impl PolishUnitConfig {
    pub fn new(read_type: ReadType, consensus_size: usize, iteration: usize) -> Self {
        Self {
            rep_num: iteration,
            seed: 342309,
            read_type,
            filter_size: consensus_size,
            consensus_size,
        }
    }
}
/// Polishing units by partial order alignment graph.
/// Note that after calling this function,
/// all the encoded reads would be removed.
/// This removal is to force users to encode reads by aligning the
/// newly poilished units again.
pub trait PolishUnit {
    fn polish_unit(self, c: &PolishUnitConfig) -> Self;
}

impl PolishUnit for DataSet {
    fn polish_unit(mut self, c: &PolishUnitConfig) -> Self {
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, vec![]))
            .collect();
        let subchunk_len: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, unit.seq().len() / 100))
            .collect();
        for read in self.encoded_reads.iter() {
            for node in read.nodes.iter() {
                if let Some(res) = pileups.get_mut(&node.unit) {
                    res.push(node);
                }
            }
        }
        let result: HashMap<_, _> = pileups
            .par_iter()
            .filter_map(|(id, pileup)| {
                let len = subchunk_len[id];
                if pileup.len() > c.filter_size {
                    consensus(pileup, len, c).map(|c| (id, c))
                } else {
                    None
                }
            })
            .collect();
        debug!("{}=>{}", pileups.len(), result.len());
        self.selected_chunks
            .retain(|unit| result.contains_key(&unit.id));
        self.selected_chunks.iter_mut().for_each(|unit| {
            unit.seq = result[&unit.id].clone();
        });
        self.encoded_reads.clear();
        self
    }
}

use definitions::Node;
fn consensus(pileup: &[&Node], len: usize, c: &PolishUnitConfig) -> Option<String> {
    let subchunks: Vec<_> = pileup
        .iter()
        .map(|n| super::local_clustering::node_to_subchunks(n, 100))
        .collect();
    let mut chunks = vec![vec![]; len];
    for sc in subchunks.iter() {
        for c in sc.iter().filter(|c| 90 < c.seq.len() && c.seq.len() < 120) {
            if c.pos < len {
                chunks[c.pos].push(c.seq.as_slice());
            }
        }
    }
    if chunks.iter().any(|cs| cs.len() < c.consensus_size) {
        None
    } else {
        let consensus: Vec<_> = chunks
            .into_iter()
            .flat_map(|cs| consensus_chunk(&cs, c))
            .collect();
        let consensus: String = consensus.into_iter().map(|n| n as char).collect();
        Some(consensus)
    }
}

fn consensus_chunk(pileup: &[&[u8]], c: &PolishUnitConfig) -> Vec<u8> {
    assert!(pileup.len() >= c.consensus_size);
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(c.seed);
    let subseq: Vec<_> = (0..c.rep_num)
        .map(|_| {
            let subchunk: Vec<_> = pileup
                .choose_multiple(&mut rng, c.consensus_size)
                .copied()
                .collect();
            POA::from_slice(&subchunk, &vec![1.; subchunk.len()], (-2, -2, &score)).consensus()
        })
        .collect();
    let subseq: Vec<_> = subseq.iter().map(|e| e.as_slice()).collect();
    POA::from_slice(&subseq, &vec![1.; subseq.len()], (-2, -2, &score)).consensus()
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        1
    } else {
        -1
    }
}
