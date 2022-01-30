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
                let len_sum: usize = pileup.iter().map(|x| x.len()).sum();
                let mean_len: usize = len_sum / pileup.len();
                let radius = c.read_type.band_width().min(mean_len / 20);
                let seqs = &pileup[..c.consensus_size.min(pileup.len())];
                let cons =
                    kiley::bialignment::polish_until_converge_banded(unit.seq(), seqs, radius);
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
                    let median_len: usize = pileup[0].len();
                    let radius = match c.read_type {
                        ReadType::CCS => (median_len / 50).max(40),
                        ReadType::None | ReadType::CLR => (median_len / 20).max(100),
                        ReadType::ONT => (median_len / 30).max(50),
                    };
                    let cons = kiley::ternary_consensus_by_chunk(seqs, radius);
                    let cons =
                        kiley::bialignment::polish_until_converge_banded(&cons, seqs, radius);
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
