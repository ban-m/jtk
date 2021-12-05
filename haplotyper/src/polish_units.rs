use definitions::ReadType;
use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct PolishUnitConfig {
    filter_size: usize,
    read_type: ReadType,
    consensus_size: usize,
}

// TODO: Remove Readtype from the argument.
impl PolishUnitConfig {
    pub fn new(read_type: ReadType, filter_size: usize, consensus_size: usize) -> Self {
        Self {
            read_type,
            filter_size,
            consensus_size,
        }
    }
    pub fn read_type(&self) -> ReadType {
        self.read_type
    }
}
/// Polishing units or Taking consensus.
/// Note that after calling this function,
/// all the encoded reads would be removed.
/// This removal is to force users to encode reads by aligning the
/// newly poilished units again.
pub trait PolishUnit {
    fn polish_unit(&mut self, c: &PolishUnitConfig);
    fn consensus_unit(&mut self, c: &PolishUnitConfig);
}

impl PolishUnit for DataSet {
    fn polish_unit(&mut self, c: &PolishUnitConfig) {
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, vec![]))
            .collect();
        for read in self.encoded_reads.iter() {
            for node in read.nodes.iter() {
                if let Some(res) = pileups.get_mut(&node.unit) {
                    res.push(node.seq());
                }
            }
        }
        self.selected_chunks
            .retain(|n| matches!(pileups.get(&n.id),Some(xs) if c.filter_size <= xs.len()));
        self.selected_chunks.par_iter_mut().for_each(|unit| {
            if let Some(pileup) = pileups.get(&unit.id) {
                let seqs = &pileup[..c.consensus_size.min(pileup.len())];
                let cons = kiley::bialignment::polish_until_converge_banded(unit.seq(), seqs, 100);
                unit.seq = String::from_utf8(cons).unwrap();
            }
        });
        // TODO: WHY? We do not need to remove these alignemnt! Just update them!
        // TODO: 1. Remove unused chunks. Change IDs in encoded read.
        // 1. Change alignment in the encoded read.
        // 2. Filtering out errorneous nodes by consider the distance and the indel patterns.
        self.encoded_reads.clear();
    }
    fn consensus_unit(&mut self, c: &PolishUnitConfig) {
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, vec![]))
            .collect();
        for read in self.encoded_reads.iter() {
            for node in read.nodes.iter() {
                if let Some(res) = pileups.get_mut(&node.unit) {
                    res.push(node.seq());
                }
            }
        }
        let result: HashMap<_, _> = pileups
            .par_iter_mut()
            .filter_map(|(id, pileup)| {
                if pileup.len() > c.filter_size {
                    let med_idx = pileup.len() / 2;
                    pileup.select_nth_unstable_by_key(med_idx, |x| x.len());
                    pileup.swap(0, med_idx);
                    let seqs = &pileup[..c.consensus_size.min(pileup.len())];
                    let cons = kiley::ternary_consensus_by_chunk(seqs, 100);
                    let cons = kiley::bialignment::polish_until_converge_banded(&cons, seqs, 100);
                    Some((id, String::from_utf8(cons).unwrap()))
                } else {
                    None
                }
            })
            .collect();
        self.selected_chunks
            .retain(|unit| result.contains_key(&unit.id));
        self.selected_chunks.iter_mut().for_each(|unit| {
            unit.seq = result[&unit.id].clone();
        });
        // TODO: WHY? We do not need to remove these alignemnt! Just update them!
        // TODO: 1. Remove unused chunks. Change IDs in encoded read.
        // 1. Change alignment in the encoded read.
        // 2. Filtering out errorneous nodes by consider the distance and the indel patterns.
        self.encoded_reads.clear();
    }
}

// use definitions::Node;
// #[allow(dead_code)]
// fn consensus(pileup: &[&Node], len: usize, c: &PolishUnitConfig) -> Option<String> {
//     let subchunks: Vec<_> = pileup
//         .iter()
//         .map(|n| super::local_clustering::node_to_subchunks(n, 100))
//         .collect();
//     let mut chunks = vec![vec![]; len];
//     for sc in subchunks.iter() {
//         for c in sc.iter().filter(|c| 90 < c.seq.len() && c.seq.len() < 120) {
//             if c.pos < len {
//                 chunks[c.pos].push(c.seq.as_slice());
//             }
//         }
//     }
//     if chunks.iter().any(|cs| cs.len() < c.consensus_size) {
//         None
//     } else {
//         let consensus: String = chunks
//             .into_iter()
//             .flat_map(|cs| consensus_chunk(&cs, c))
//             .map(|n| n as char)
//             .collect();
//         Some(consensus)
//     }
// }

// fn consensus_chunk(pileup: &[&[u8]], c: &PolishUnitConfig) -> Vec<u8> {
//     assert!(pileup.len() >= c.consensus_size);
//     let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(c.seed);
//     let subseq: Vec<_> = (0..c.rep_num)
//         .map(|_| {
//             let subchunk: Vec<_> = pileup
//                 .choose_multiple(&mut rng, c.consensus_size)
//                 .copied()
//                 .collect();
//             POA::from_slice(&subchunk, &vec![1.; subchunk.len()], (-2, -2, &score)).consensus()
//         })
//         .collect();
//     let subseq: Vec<_> = subseq.iter().map(|e| e.as_slice()).collect();
//     POA::from_slice(&subseq, &vec![1.; subseq.len()], (-2, -2, &score)).consensus()
// }

// fn score(x: u8, y: u8) -> i32 {
//     if x == y {
//         1
//     } else {
//         -1
//     }
// }
