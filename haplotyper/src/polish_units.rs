use super::local_clustering::ReadType;
use definitions::*;
use poa_hmm::POA;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct PolishUnitConfig {
    read_type: ReadType,
    consensus_size: usize,
    rep_num: usize,
    seed: u64,
}

impl PolishUnitConfig {
    pub fn new(readtype: &str, consensus_size: usize) -> Self {
        let read_type = match readtype {
            "ONT" => ReadType::ONT,
            "CCS" => ReadType::CCS,
            "CLR" => ReadType::CLR,
            _ => unreachable!(),
        };
        let (rep_num, seed) = (3, 349309);
        Self {
            rep_num,
            seed,
            read_type,
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
                    res.push(node.seq.as_bytes());
                }
            }
        }
        let result: HashMap<_, _> = pileups
            .par_iter()
            .map(|(id, pileup)| {
                let consensus = consensus(pileup, c);
                (id, consensus)
            })
            .collect();
        self.selected_chunks.iter_mut().for_each(|unit| {
            if let Some(res) = result.get(&unit.id) {
                unit.seq = res.clone();
            }
        });
        self.encoded_reads.clear();
        self
    }
}

fn consensus(pileup: &[&[u8]], c: &PolishUnitConfig) -> String {
    if pileup.len() <= c.consensus_size {
        String::from_utf8_lossy(&POA::from_slice_default(&pileup).consensus()).to_string()
    } else {
        use rand::{seq::SliceRandom, SeedableRng};
        use rand_xoshiro::Xoshiro256PlusPlus;
        let subchunks = c.rep_num * pileup.len() / c.consensus_size;
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(c.seed);
        let subseq: Vec<_> = (0..subchunks)
            .map(|_| {
                let subchunk: Vec<_> = pileup
                    .choose_multiple(&mut rng, c.consensus_size)
                    .copied()
                    .collect();
                POA::from_slice_default(&subchunk).consensus()
            })
            .collect();
        String::from_utf8_lossy(&POA::from_vec_default(&subseq).consensus()).to_string()
    }
}
