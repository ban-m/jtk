// use log::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    // How many draft-consensus we generate to get a consensus.
    // pub consensus_size: usize,
    // How many reads we use to get a draft-consensus
    // pub consensus_number: usize,
    // Number of cluster. We run assembly-based clustering with k.
    pub cluster_num: usize,
    // How many assembly-based clustering we exec.
    pub repeat_number: usize,
    pub band_size: usize,
    pub seed: u64,
}

impl ClusteringConfig {
    pub fn new(cluster_num: usize, repeat_number: usize, band_size: usize, seed: u64) -> Self {
        Self {
            cluster_num,
            repeat_number,
            band_size,
            seed,
        }
    }
}

use definitions::Node;
pub fn clustering_by_assemble(data: &mut [Node], config: &ClusteringConfig) -> usize {
    let seqs: Vec<_> = data.iter().map(|x| x.seq().to_vec()).collect();
    let assignments = clustering(&seqs, config);
    data.iter_mut()
        .zip(assignments)
        .for_each(|(x, a)| x.cluster = a as u64);
    config.cluster_num
}

pub fn clustering_m(seqs: &[Vec<u8>], config: &ClusteringConfig) -> Vec<u8> {
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(config.seed);
    (0..4)
        .map(|_| clustering_from(seqs, config, &mut rng))
        .min_by_key(|x| x.0)
        .map(|x| x.1.iter().map(|&x| x as u8).collect::<Vec<_>>())
        .unwrap_or_else(|| vec![0; seqs.len()])
}

fn clustering_from<R: Rng>(
    seqs: &[Vec<u8>],
    config: &ClusteringConfig,
    rng: &mut R,
) -> (u32, Vec<usize>) {
    let mut assignments: Vec<_> = (0..seqs.len())
        .map(|_| rng.gen_range(0..config.cluster_num))
        .collect();
    let template =
        kiley::ternary_consensus(seqs, rng.gen(), config.repeat_number, config.band_size);
    let mut consensus = vec![template; config.cluster_num];
    consensus = polish_consensus(&assignments, seqs, &consensus, config);
    let mut dist: u32 = get_total_edit_distance(&assignments, seqs, &consensus);
    debug!("Dist:{}", dist);
    use rand::seq::SliceRandom;
    for _ in 0.. {
        for idx in rand::seq::index::sample(rng, seqs.len(), seqs.len()) {
            // Update assignments.
            let mut mins: Vec<_> = consensus
                .iter()
                .map(|c| edlib_sys::global_dist(c, &seqs[idx]))
                .enumerate()
                .collect();
            mins.shuffle(rng);
            let new_asn = mins.iter().min_by_key(|x| x.1).map(|x| x.0).unwrap();
            assignments[idx] = new_asn;
            consensus = polish_consensus(&assignments, seqs, &consensus, config);
        }
        let new_dist = get_total_edit_distance(&assignments, seqs, &consensus);
        debug!("Dist:{}", new_dist);
        if dist <= new_dist {
            break;
        }
        dist = new_dist;
    }
    (dist, assignments)
}

pub fn clustering(seqs: &[Vec<u8>], config: &ClusteringConfig) -> Vec<u8> {
    debug!("Start");
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(config.seed);
    let mut assignments: Vec<_> = (0..seqs.len())
        .map(|_| rng.gen_range(0..config.cluster_num))
        .collect();
    let template =
        kiley::ternary_consensus(seqs, rng.gen(), config.repeat_number, config.band_size);
    let mut consensus = vec![template; config.cluster_num];
    consensus = polish_consensus(&assignments, seqs, &consensus, config);
    let mut dist: u32 = get_total_edit_distance(&assignments, seqs, &consensus);
    debug!("Dist:{}", dist);
    use rand::seq::SliceRandom;
    for _ in 0.. {
        for idx in rand::seq::index::sample(&mut rng, seqs.len(), seqs.len()) {
            // Update assignments.
            let mut mins: Vec<_> = consensus
                .iter()
                .map(|c| edlib_sys::global_dist(c, &seqs[idx]))
                .enumerate()
                .collect();
            mins.shuffle(&mut rng);
            let new_asn = mins.iter().min_by_key(|x| x.1).map(|x| x.0).unwrap();
            assignments[idx] = new_asn;
            // Update parameters.
            consensus = polish_consensus(&assignments, seqs, &consensus, config);
        }
        let new_dist = get_total_edit_distance(&assignments, seqs, &consensus);
        debug!("Dist:{}", new_dist);
        if dist <= new_dist {
            break;
        }
        dist = new_dist;
    }

    assignments.iter().map(|&x| x as u8).collect()
}

fn polish_consensus(
    assignments: &[usize],
    seqs: &[Vec<u8>],
    consensus: &[Vec<u8>],
    config: &ClusteringConfig,
) -> Vec<Vec<u8>> {
    let mut slots = vec![vec![]; consensus.len()];
    for (&asn, seq) in assignments.iter().zip(seqs.iter()) {
        slots[asn].push(seq.as_slice());
    }
    // This unwarap might also dangerous.
    consensus
        .iter()
        .zip(slots)
        .map(|(c, xs)| kiley::polish_until_converge_banded(c, &xs, config.band_size))
        .collect()
}

fn get_total_edit_distance(assignments: &[usize], seqs: &[Vec<u8>], consensus: &[Vec<u8>]) -> u32 {
    assignments
        .iter()
        .zip(seqs.iter())
        .map(|(&asn, x)| edlib_sys::global_dist(&consensus[asn], x))
        .sum()
}

// impl ClusteringConfig {
//     pub fn default(chain_len: usize, cluster_num: usize) -> Self {
//         Self {
//             consensus_size: 10,
//             consensus_number: 10,
//             chain_len,
//             seed: 0,
//             cluster_num,
//             trial_number: 10,
//             repeat_number: 20,
//             match_score: 3,
//             mismat_score: -5,
//             gap_open: -3,
//             gap_extend: -1,
//             kmer_match: 3,
//             band_size: 10,
//         }
//     }
// }

// pub fn clustering_by_assemble(data: &mut Vec<ChunkedUnit>, config: &ClusteringConfig) -> usize {
//     let mut labels = vec![vec![]; data.len()];
//     let tn = config.trial_number as u64;
//     use rayon::prelude::*;
//     let assignments: Vec<_> = (0..tn)
//         .into_par_iter()
//         .map(|i| clustering_with(data, config, i + config.seed))
//         .collect();
//     // for assignment in (0..tn).map(|i| clustering_with(data, config, i + config.seed)) {
//     for assignment in assignments {
//         for (i, asn) in assignment.into_iter().enumerate() {
//             labels[i].push(asn);
//         }
//     }
//     for (i, d) in labels.iter().enumerate() {
//         trace!("{}\t{:?}", i, d);
//     }
//     assert!(labels.iter().all(|xs| xs.len() == config.trial_number));
//     let (asn, _) = ((config.cluster_num / 2).max(1)..=config.cluster_num * 2)
//         .map(|cl| {
//             let param = (cl as u64, config.trial_number, config.cluster_num);
//             em_clustering(&labels, cl, param)
//         })
//         .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
//         .unwrap();
//     let cluster_num = *asn.iter().max().unwrap();
//     for (a, x) in asn.into_iter().zip(data.iter_mut()) {
//         x.cluster = a;
//     }
//     cluster_num
// }

// fn clustering_with(data: &[ChunkedUnit], config: &ClusteringConfig, seed: u64) -> Vec<usize> {
//     let k = config.cluster_num;
//     let chain = config.chain_len;
//     let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
//     let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen::<usize>() % k).collect();
//     let templates: Vec<_> = {
//         let mut chunks = vec![vec![]; chain];
//         for read in data.iter() {
//             for c in read.chunks.iter() {
//                 chunks[c.pos].push(c.seq.as_slice());
//             }
//         }
//         // for (i, c) in chunks.iter().enumerate() {
//         //     debug!("{}\t{}", i, c.len());
//         // }
//         take_consensus(&chunks, &mut rng, config)
//     };
//     let mut aligner = Aligner::new(
//         config.gap_open,
//         config.gap_extend,
//         |x, y| {
//             if x == y {
//                 config.match_score
//             } else {
//                 config.mismat_score
//             }
//         },
//         config.kmer_match,
//         config.band_size,
//     );
//     let choises: Vec<_> = (0..k).collect();
//     let mut count = vec![0; k];
//     let mut score;
//     for &asn in assignments.iter() {
//         count[asn] += 1;
//     }
//     let mut changes = vec![];
//     let mut chunks = vec![vec![vec![]; chain]; k];
//     //for _i in 0..config.repeat_number {
//     loop {
//         for (i, t) in templates.iter().enumerate().filter(|x| !x.1.is_empty()) {
//             for slot in chunks.iter_mut().take(k) {
//                 slot[i].push(t.as_slice());
//             }
//         }
//         for (read, &asn) in data.iter().zip(assignments.iter()) {
//             for c in read.chunks.iter() {
//                 chunks[asn][c.pos].push(c.seq.as_slice());
//             }
//         }
//         let templates: Vec<Vec<_>> = chunks
//             .iter()
//             .map(|cs| take_consensus(cs, &mut rng, config))
//             .collect();
//         chunks
//             .iter_mut()
//             .for_each(|xs| xs.iter_mut().for_each(|x| x.clear()));
//         let (current_score, changed) = data
//             .iter()
//             .zip(assignments.iter_mut())
//             .map(|(read, asn)| {
//                 let scores: Vec<_> = templates
//                     .iter()
//                     .map(|t| {
//                         read.chunks
//                             .iter()
//                             .filter(|c| !t[c.pos].is_empty())
//                             .map(|c| aligner.global(&t[c.pos], &c.seq).score)
//                             .sum::<i32>()
//                     })
//                     .collect();
//                 let probs: Vec<_> = scores
//                     .iter()
//                     .map(|x| {
//                         scores
//                             .iter()
//                             .map(|y| ((y - x) as f64).exp())
//                             .sum::<f64>()
//                             .recip()
//                     })
//                     .collect();
//                 let new = *choises.choose_weighted(&mut rng, |&k| probs[k]).unwrap();
//                 let changed = (*asn != new) as usize;
//                 *asn = new;
//                 (*scores.iter().max().unwrap(), changed)
//             })
//             .fold((0, 0), |(x, y), (a, b)| (x + a, y + b));
//         changes.push(changed);
//         count = vec![0; k];
//         for &asn in assignments.iter() {
//             count[asn] += 1;
//         }
//         score = current_score;
//         if changed < data.len() / 10 {
//             break;
//         }
//     }
//     trace!("HISTORY\t{:?}", changes);
//     trace!("LK\t{}\t{:?}", score, count);
//     assignments
// }

// fn take_consensus<R: Rng>(chunks: &[Vec<&[u8]>], r: &mut R, c: &ClusteringConfig) -> Vec<Vec<u8>> {
//     chunks
//         .iter()
//         .map(|cs| match cs.is_empty() {
//             true => vec![],
//             false => consensus(cs, r, c),
//         })
//         .collect()
// }

// fn consensus<R: Rng>(cs: &[&[u8]], r: &mut R, c: &ClusteringConfig) -> Vec<u8> {
//     let mut cs = cs.to_vec();
//     let take_num = cs.len().min(c.consensus_size);
//     let cs: Vec<_> = (0..c.consensus_number)
//         .map(|_| {
//             cs.shuffle(r);
//             poa_hmm::POA::from_slice_default(&cs[..take_num]).consensus()
//         })
//         .collect();
//     poa_hmm::POA::from_vec_default(&cs).consensus()
// }
