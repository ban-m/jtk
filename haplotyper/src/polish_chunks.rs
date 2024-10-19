use definitions::ReadType;
use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct PolishChunkConfig {
    filter_size: usize,
    read_type: ReadType,
    #[allow(dead_code)]
    consensus_size: usize,
}

impl PolishChunkConfig {
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
/// Polishing chunks or Taking consensus.
/// Note that after calling this function,
/// all the encoded reads would be removed.
/// This removal is to force users to encode reads by aligning the
/// newly poilished chunks again.
pub trait PolishChunk {
    fn polish_chunk(&mut self, c: &PolishChunkConfig);
    fn consensus_chunk(&mut self, c: &PolishChunkConfig);
}

impl PolishChunk for DataSet {
    fn polish_chunk(&mut self, c: &PolishChunkConfig) {
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|chunk| (chunk.id, vec![]))
            .collect();
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if let Some(res) = pileups.get_mut(&node.chunk) {
                    res.push(node);
                }
            }
        }
        let chunk_seqs: HashMap<_, _> = self.selected_chunks.iter().map(|u| (u.id, u)).collect();
        let mut polished_nodes: HashMap<_, _> = pileups
            .into_par_iter()
            .filter(|x| c.filter_size < x.1.len())
            .map(|(id, mut pileup)| {
                let chunk = chunk_seqs.get(&id).unwrap();
                let radius = c.read_type().band_width(chunk.seq().len());
                pileup.sort_by_cached_key(|node| {
                    let (_, aln, _) = node.recover(chunk);
                    aln.iter().filter(|&&x| x != b'|').count()
                });
                let (seqs, mut ops): (Vec<_>, Vec<_>) = pileup
                    .iter()
                    .map(|n| (n.seq(), crate::misc::ops_to_kiley(&n.cigar)))
                    .unzip();
                use kiley::bialignment::guided::polish_until_converge_with;
                let cons = polish_until_converge_with(chunk.seq(), &seqs, &mut ops, radius);
                pileup
                    .iter_mut()
                    .zip(ops)
                    .for_each(|(n, ops)| n.cigar = crate::misc::kiley_op_to_ops(&ops));
                (id, cons)
            })
            .collect();
        self.selected_chunks
            .retain(|n| polished_nodes.contains_key(&n.id));
        self.selected_chunks.iter_mut().for_each(|n| {
            n.seq = polished_nodes.remove(&n.id).unwrap().into();
        });
        use std::collections::HashSet;
        let is_in: HashSet<_> = self.selected_chunks.iter().map(|c| c.id).collect();
        self.encoded_reads.iter_mut().for_each(|read| {
            let mut idx = 0;
            loop {
                match read.nodes.get(idx) {
                    Some(n) if is_in.contains(&n.chunk) => idx += 1,
                    Some(_) => read.remove(idx),
                    None => return,
                }
            }
        });
    }
    fn consensus_chunk(&mut self, c: &PolishChunkConfig) {
        // First, take consensus.
        let mut pileups: HashMap<_, Vec<_>> = self
            .selected_chunks
            .iter()
            .map(|chunk| (chunk.id, vec![]))
            .collect();
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if let Some(res) = pileups.get_mut(&node.chunk) {
                    res.push(node);
                }
            }
        }
        pileups.retain(|_, pileup| c.filter_size < pileup.len());
        self.selected_chunks
            .retain(|chunk| pileups.contains_key(&chunk.id));
        self.selected_chunks.iter_mut().for_each(|chunk| {
            let pileup = pileups.get_mut(&chunk.id).unwrap();
            let med_idx = pileup.len() / 2;
            pileup.select_nth_unstable_by_key(med_idx, |x| x.seq().len());
            pileup.swap(0, med_idx);
            chunk.seq = pileup[0].seq().to_vec().into();
            pileup.iter_mut().for_each(|node| {
                let mode = edlib_sys::AlignMode::Global;
                let task = edlib_sys::AlignTask::Alignment;
                let aln = edlib_sys::align(node.seq(), chunk.seq(), mode, task);
                let ops = crate::misc::edlib_to_kiley(aln.operations().unwrap());
                node.cigar = crate::misc::kiley_op_to_ops(&ops);
            });
        });
        self.polish_chunk(c);
    }
}
