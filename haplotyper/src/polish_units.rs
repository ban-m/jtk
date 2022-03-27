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
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if let Some(res) = pileups.get_mut(&node.unit) {
                    res.push(node);
                }
            }
        }
        let unit_seqs: HashMap<_, _> = self.selected_chunks.iter().map(|u| (u.id, u)).collect();
        let mut polished_nodes: HashMap<_, _> = pileups
            .into_par_iter()
            .filter(|x| c.filter_size < x.1.len())
            .map(|(id, mut pileup)| {
                let unit = unit_seqs.get(&id).unwrap();
                let radius = c.read_type.band_width(unit.seq().len());
                pileup.sort_by_cached_key(|node| {
                    let (_, aln, _) = node.recover(unit);
                    aln.iter().filter(|&&x| x != b'|').count()
                });
                let (seqs, mut ops): (Vec<_>, Vec<_>) = pileup
                    .iter()
                    .map(|n| (n.seq(), crate::local_clustering::ops_to_kiley_ops(&n.cigar)))
                    .unzip();
                use kiley::bialignment::guided::polish_until_converge_with;
                let cons = polish_until_converge_with(unit.seq(), &seqs, &mut ops, radius);
                pileup
                    .iter_mut()
                    .zip(ops)
                    .for_each(|(n, ops)| n.cigar = crate::encode::compress_kiley_ops(&ops));
                (id, cons)
            })
            .collect();
        self.selected_chunks
            .retain(|n| polished_nodes.contains_key(&n.id));
        self.selected_chunks.iter_mut().for_each(|n| {
            let seq = polished_nodes.remove(&n.id).unwrap();
            n.seq = String::from_utf8(seq).unwrap();
        });
        use std::collections::HashSet;
        let is_in: HashSet<_> = self.selected_chunks.iter().map(|c| c.id).collect();
        self.encoded_reads.iter_mut().for_each(|read| {
            let mut idx = 0;
            loop {
                match read.nodes.get(idx) {
                    Some(n) if is_in.contains(&n.unit) => idx += 1,
                    Some(_) => read.remove(idx),
                    None => return,
                }
            }
        });
        // Remove node from reads.
        // self.selected_chunks
        //     .retain(|n| c.filter_size < pileups[&n.id].len());
        // self.selected_chunks.par_iter_mut().for_each(|unit| {
        //     // Never panic.
        //     let pileup = &pileups[&unit.id];
        //     let len_sum: usize = pileup.iter().map(|x| x.seq().len()).sum();
        //     let mean_len: usize = len_sum / pileup.len();
        //     let radius = c.read_type.band_width(mean_len);
        //     let mut seqs_with_diff: Vec<_> = pileup
        //         .iter()
        //         .map(|node| {
        //             let (_, aln, _) = node.recover(unit);
        //             let diff: usize = aln.iter().filter(|&&x| x != b'|').count();
        //             (diff, node.seq())
        //         })
        //         .collect();
        //     seqs_with_diff.sort_by_key(|x| x.0);
        //     let seqs: Vec<_> = seqs_with_diff
        //         .into_iter()
        //         .take(c.consensus_size)
        //         .map(|x| x.1)
        //         .collect();
        //     let cons = kiley::bialignment::guided::polish_until_converge(unit.seq(), &seqs, radius);
        //     unit.seq = String::from_utf8(cons).unwrap();
        // });
        // // TODO: WHY? We do not need to remove these alignemnt! Just update them!
        // // TODO: 1. Remove unused chunks. Change IDs in encoded read.
        // // 1. Change alignment in the encoded read.
        // // 2. Filtering out errorneous nodes by consider the distance and the indel patterns.
        // self.encoded_reads.clear();
    }
    fn consensus_unit(&mut self, c: &PolishUnitConfig) {
        let ed_ops = [
            kiley::Op::Match,
            kiley::Op::Ins,
            kiley::Op::Del,
            kiley::Op::Mismatch,
        ];
        // First, take consensus.
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|unit| (unit.id, vec![]))
            .collect();
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if let Some(res) = pileups.get_mut(&node.unit) {
                    res.push(node);
                }
            }
        }
        pileups.retain(|_, pileup| c.filter_size < pileup.len());
        let result: HashMap<_, _> = pileups
            .par_iter_mut()
            .map(|(id, pileup)| {
                let med_idx = pileup.len() / 2;
                pileup.select_nth_unstable_by_key(med_idx, |x| x.seq().len());
                pileup.swap(0, med_idx);
                let draft = {
                    let seqs: Vec<_> = pileup
                        .iter()
                        .map(|x| x.seq())
                        .take(c.consensus_size)
                        .collect();
                    let median_len: usize = seqs[0].len();
                    let radius = c.read_type.band_width(median_len);
                    kiley::ternary_consensus_by_chunk(&seqs, radius)
                };
                pileup.iter_mut().for_each(|node| {
                    let mode = edlib_sys::AlignMode::Global;
                    let task = edlib_sys::AlignTask::Alignment;
                    let aln = edlib_sys::edlib_align(node.seq(), &draft, mode, task);
                    let aln = aln.operations.unwrap();
                    let k_ops: Vec<_> = aln.iter().map(|&op| ed_ops[op as usize]).collect();
                    node.cigar = crate::encode::compress_kiley_ops(&k_ops);
                });
                (id, String::from_utf8(draft).unwrap())
            })
            .collect();
        self.selected_chunks
            .retain(|unit| result.contains_key(&unit.id));
        self.selected_chunks.iter_mut().for_each(|unit| {
            unit.seq = result[&unit.id].clone();
        });
        self.polish_unit(c);
    }
}
