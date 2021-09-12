use super::ChunkedUnit;
use super::ClusteringConfig;
use poa_hmm::POA;
use rand::seq::SliceRandom;
use rand::Rng;
use rayon::prelude::*;
const CONSENSUS_NUMBER: usize = 15;
const MAX_COVERAGE: usize = 40;
pub fn get_models<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    (cluster_num, chain_len): (usize, usize),
    rng: &mut R,
    c: &ClusteringConfig<F>,
    use_position: &[bool],
    picked: Option<usize>,
) -> Vec<Vec<POA>> {
    let mut reference = vec![vec![]; chain_len];
    let mut chunks: Vec<Vec<Vec<&[u8]>>> = vec![vec![vec![]; chain_len]; cluster_num];
    for (_, read) in data.iter().enumerate().filter(|&(i, _)| Some(i) != picked) {
        for chunk in read.chunks.iter().filter(|chunk| use_position[chunk.pos]) {
            chunks[read.cluster][chunk.pos].push(chunk.seq.as_slice());
            reference[chunk.pos].push(chunk.seq.as_slice());
        }
    }
    chunks.iter_mut().for_each(|cluster| {
        cluster.iter_mut().for_each(|cs| {
            cs.shuffle(rng);
        });
    });
    let reference: Vec<_> = reference
        .into_par_iter()
        .zip(use_position.par_iter())
        .map(|(cs, &b)| {
            if b && !cs.is_empty() {
                Some(construct_poa(&None, &cs, c))
            } else {
                None
            }
        })
        .collect();
    chunks
        .into_par_iter()
        .map(|cluster| {
            cluster
                .into_iter()
                .zip(use_position.iter())
                .zip(reference.iter())
                .map(|((cs, &b), refr)| {
                    if b && !cs.is_empty() {
                        construct_poa(refr, &cs, c)
                    } else {
                        POA::default()
                    }
                })
                .collect()
        })
        .collect()
}

pub fn construct_poa<F>(base: &Option<POA>, cs: &[&[u8]], c: &ClusteringConfig<F>) -> POA
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let &super::AlignmentParameters {
        ins,
        del,
        ref score,
    } = &c.alnparam;
    let param = (ins, del, score);
    use definitions::ReadType;
    let base = base
        .as_ref()
        .map(|e| POA::new(&e.consensus(), 1.))
        .unwrap_or_default();
    match c.read_type {
        ReadType::CCS => base.update(cs, &vec![1.; cs.len()], param),
        _ => {
            let mut cs: Vec<_> = cs.iter().copied().collect();
            let mut rng: rand_xoshiro::Xoroshiro128StarStar = rand::SeedableRng::seed_from_u64(924);
            let cs: Vec<_> = (0..CONSENSUS_NUMBER)
                .map(|_| {
                    cs.shuffle(&mut rng);
                    let max_len = cs.iter().map(|s| s.len()).max().unwrap_or(0);
                    let node_num_thr = (max_len as f64 * 1.5).floor() as usize;
                    cs.iter()
                        .take(MAX_COVERAGE)
                        .fold(base.clone(), |x, y| {
                            let res = if x.nodes().len() > node_num_thr {
                                x.add(y, 1., param).remove_node(0.4)
                            } else {
                                x.add(y, 1., param)
                            };
                            res
                        })
                        .remove_node(0.4)
                        .finalize()
                        .consensus()
                })
                .collect();
            let cs: Vec<_> = cs.iter().map(|cs| cs.as_slice()).collect();
            POA::default().update_thr(&cs, &vec![1.; cs.len()], param, 0.8, 1.5)
        }
    }
}
