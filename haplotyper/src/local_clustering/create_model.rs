use super::ChunkedUnit;
use super::ClusteringConfig;
use poa_hmm::POA;
use rand::seq::SliceRandom;
use rand::Rng;
use rayon::prelude::*;
pub fn get_models<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
    use_position: &[bool],
    picked: Option<usize>,
) -> Vec<Vec<POA>> {
    let mut chunks: Vec<Vec<Vec<&[u8]>>> = vec![vec![vec![]; chain_len]; c.cluster_num];
    for (_, read) in data.iter().enumerate().filter(|&(i, _)| Some(i) != picked) {
        for chunk in read.chunks.iter().filter(|chunk| use_position[chunk.pos]) {
            if read.cluster < c.cluster_num {
                chunks[read.cluster][chunk.pos].push(chunk.seq.as_slice());
            } else {
                // Wild card
                chunks
                    .iter_mut()
                    .for_each(|cs| cs[chunk.pos].push(chunk.seq.as_slice()));
            }
        }
    }
    // TODO:Modify so that we need not to allocate unnessesary arraies.
    let chunks: Vec<Vec<Vec<Vec<u8>>>> = chunks
        .into_iter()
        .map(|cluster| {
            cluster
                .into_iter()
                .map(|mut cs| {
                    cs.shuffle(rng);
                    cs.iter().map(|e| e.to_vec()).collect()
                })
                .collect()
        })
        .collect();
    chunks
        .into_par_iter()
        .map(|cluster| {
            cluster
                .into_iter()
                .zip(use_position.iter())
                .map(|(cs, &b)| {
                    if b && !cs.is_empty() {
                        construct_poa(&cs, c)
                    } else {
                        POA::default()
                    }
                })
                .collect()
        })
        .collect()
}

fn construct_poa<F>(cs: &[Vec<u8>], c: &ClusteringConfig<F>) -> POA
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let &super::AlignmentParameters {
        ins,
        del,
        ref score,
    } = &c.alnparam;
    let param = (ins, del, score);
    // let cs = &cs[..cs.len().min(c.sample_num)];
    POA::from_vec(cs, &vec![1.; cs.len()], param)
    // use super::config::ReadType;
    // match c.read_type {
    //     ReadType::CCS => POA::from_vec(cs, &vec![1.; cs.len()], param),
    //     _ => {
    //         let max_len = cs.iter().map(|s| s.len()).max().unwrap_or(0);
    //         let node_num_thr = (max_len as f64 * 1.5).floor() as usize;
    //         cs.iter()
    //             .fold(POA::default(), |x, y| {
    //                 let res = if x.nodes().len() > node_num_thr {
    //                     x.add(y, 1., param).remove_node(0.5)
    //                 } else {
    //                     x.add(y, 1., param)
    //                 };
    //                 res
    //             })
    //             .remove_node(0.7)
    //             .finalize()
    //     }
    // }
}
