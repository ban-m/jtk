use super::local_clustering::ReadType;
use definitions::*;
use poa_hmm::POA;
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
    pub fn new(readtype: &str, consensus_size: usize, filter_size: usize) -> Self {
        let read_type = match readtype {
            "ONT" => ReadType::ONT,
            "CCS" => ReadType::CCS,
            "CLR" => ReadType::CLR,
            _ => unreachable!(),
        };
        let (rep_num, seed) = (5, 349309);
        Self {
            rep_num,
            seed,
            read_type,
            filter_size,
            consensus_size,
        }
    }
}
/// Polishing units by partial order alignment graph.
/// Note that after calling this function,
/// all the encoded reads would be removed.
/// This removal is to force users to encode reads by aligning the
/// newly poilished units again.
pub trait PolishUnits {
    fn polish_units(self, c: &PolishUnitConfig) -> Self;
}

impl PolishUnits for DataSet {
    fn polish_units(mut self, c: &PolishUnitConfig) -> Self {
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, vec![]))
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
            .filter(|(_, pileup)| pileup.len() > c.filter_size)
            .filter_map(|(id, pileup)| consensus(pileup, c).map(|c| (id, c)))
            .collect();
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
fn consensus(pileup: &[&Node], c: &PolishUnitConfig) -> Option<String> {
    let subchunks: Vec<_> = pileup
        .iter()
        .map(|n| super::local_clustering::node_to_subchunks(n, 100))
        .collect();
    let max_pos = subchunks
        .iter()
        .filter_map(|n| n.iter().map(|n| n.pos).max())
        .max()
        .unwrap();
    let mut chunks = vec![vec![]; max_pos + 1];
    for sc in subchunks.iter() {
        for c in sc.iter() {
            chunks[c.pos].push(c.seq.as_slice());
        }
    }
    if chunks.iter().any(|cs| cs.len() < c.consensus_size) {
        None
    } else {
        let consensus: String = chunks
            .into_iter()
            .flat_map(|cs| consensus_chunk(&cs, c))
            .map(|n| n as char)
            .collect();
        Some(consensus)
    }
}

fn consensus_chunk(pileup: &[&[u8]], c: &PolishUnitConfig) -> Vec<u8> {
    assert!(pileup.len() >= c.consensus_size);
    use rand::{seq::SliceRandom, SeedableRng};
    use rand_xoshiro::Xoshiro256PlusPlus;
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
