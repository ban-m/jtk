use super::ChunkedUnit;
use super::ClusteringConfig;
use poa_hmm::POA;
use rand::distributions::Standard;
use rand::seq::SliceRandom;
use rand::Rng;
//use rayon::prelude::*;

fn select<R: Rng>(choises: &[usize], rng: &mut R, cl: usize, pick: f64) -> usize {
    *choises
        .choose_weighted(rng, |&k| if k == cl { 1. + pick } else { pick })
        .unwrap()
}

pub fn get_models<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
    use_position: &[bool],
    update_data: &[bool],
) -> Vec<Vec<POA>> {
    let mut chunks: Vec<_> = vec![vec![vec![]; chain_len]; c.cluster_num];
    let choises: Vec<usize> = (0..c.cluster_num).collect();
    for (read, _) in data.iter().zip(update_data).filter(|&(_, b)| !b) {
        let chosen = select(&choises, rng, read.cluster, c.gibbs_prior);
        for chunk in read.chunks.iter().filter(|chunk| use_position[chunk.pos]) {
            chunks[chosen][chunk.pos].push(chunk.seq.as_slice());
        }
    }
    let seeds: Vec<_> = rng.sample_iter(Standard).take(chain_len).collect();
    let &super::AlignmentParameters {
        ins,
        del,
        ref score,
    } = &c.alnparam;
    let param = (ins, del, score);
    chunks
        .iter()
        .map(|cluster| {
            cluster
                .iter()
                .zip(seeds.iter())
                .zip(use_position.iter())
                .map(|((cs, &s), &b)| {
                    if b {
                        POA::generate_banded(cs, param, 10, s)
                    } else {
                        POA::default()
                    }
                })
                .collect()
        })
        .collect()
}
