use super::em_clustering;
use crate::local_clustering::eread::ChunkedUnit;
use bio::alignment::pairwise::banded::*;
use log::*;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
#[derive(Debug, Clone)]
pub struct ClusteringConfig {
    // How many draft-consensus we generate to get a consensus.
    pub consensus_size: usize,
    // How many reads we use to get a draft-consensus
    pub consensus_number: usize,
    // Chain length.
    pub chain_len: usize,
    // seed number
    pub seed: u64,
    // Number of cluster. We run assembly-based clustering with k.
    pub cluster_num: usize,
    // How many assembly-based clustering we exec.
    pub trial_number: usize,
    // How many iteration we exec in a assembly-based clustering.
    pub repeat_number: usize,
    // match score in alignment
    pub match_score: i32,
    pub mismat_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub kmer_match: usize,
    pub band_size: usize,
}

impl ClusteringConfig {
    pub fn default(chain_len: usize, cluster_num: usize) -> Self {
        Self {
            consensus_size: 10,
            consensus_number: 10,
            chain_len,
            seed: 0,
            cluster_num,
            trial_number: 10,
            repeat_number: 20,
            match_score: 3,
            mismat_score: -5,
            gap_open: -3,
            gap_extend: -1,
            kmer_match: 3,
            band_size: 10,
        }
    }
}

pub fn clustering_by_assemble(data: &mut Vec<ChunkedUnit>, config: &ClusteringConfig) -> usize {
    let mut labels = vec![vec![]; data.len()];
    let tn = config.trial_number as u64;
    use rayon::prelude::*;
    let assignments: Vec<_> = (0..tn)
        .into_par_iter()
        .map(|i| clustering_with(data, config, i + config.seed))
        .collect();
    // for assignment in (0..tn).map(|i| clustering_with(data, config, i + config.seed)) {
    for assignment in assignments {
        for (i, asn) in assignment.into_iter().enumerate() {
            labels[i].push(asn);
        }
    }
    for (i, d) in labels.iter().enumerate() {
        trace!("{}\t{:?}", i, d);
    }
    assert!(labels.iter().all(|xs| xs.len() == config.trial_number));
    let (asn, _) = ((config.cluster_num / 2).max(1)..=config.cluster_num * 2)
        .map(|cl| {
            let param = (cl as u64, config.trial_number, config.cluster_num);
            em_clustering(&labels, cl, param)
        })
        .max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
        .unwrap();
    let cluster_num = *asn.iter().max().unwrap();
    for (a, x) in asn.into_iter().zip(data.iter_mut()) {
        x.cluster = a;
    }
    cluster_num
}

fn clustering_with(data: &[ChunkedUnit], config: &ClusteringConfig, seed: u64) -> Vec<usize> {
    let k = config.cluster_num;
    let chain = config.chain_len;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen::<usize>() % k).collect();
    let templates: Vec<_> = {
        let mut chunks = vec![vec![]; chain];
        for read in data.iter() {
            for c in read.chunks.iter() {
                chunks[c.pos].push(c.seq.as_slice());
            }
        }
        // for (i, c) in chunks.iter().enumerate() {
        //     debug!("{}\t{}", i, c.len());
        // }
        take_consensus(&chunks, &mut rng, config)
    };
    let mut aligner = Aligner::new(
        config.gap_open,
        config.gap_extend,
        |x, y| {
            if x == y {
                config.match_score
            } else {
                config.mismat_score
            }
        },
        config.kmer_match,
        config.band_size,
    );
    let choises: Vec<_> = (0..k).collect();
    let mut count = vec![0; k];
    let mut score;
    for &asn in assignments.iter() {
        count[asn] += 1;
    }
    let mut changes = vec![];
    let mut chunks = vec![vec![vec![]; chain]; k];
    //for _i in 0..config.repeat_number {
    loop {
        for (i, t) in templates.iter().enumerate().filter(|x| !x.1.is_empty()) {
            for slot in chunks.iter_mut().take(k) {
                slot[i].push(t.as_slice());
            }
        }
        for (read, &asn) in data.iter().zip(assignments.iter()) {
            for c in read.chunks.iter() {
                chunks[asn][c.pos].push(c.seq.as_slice());
            }
        }
        let templates: Vec<Vec<_>> = chunks
            .iter()
            .map(|cs| take_consensus(cs, &mut rng, config))
            .collect();
        chunks
            .iter_mut()
            .for_each(|xs| xs.iter_mut().for_each(|x| x.clear()));
        let (current_score, changed) = data
            .iter()
            .zip(assignments.iter_mut())
            .map(|(read, asn)| {
                let scores: Vec<_> = templates
                    .iter()
                    .map(|t| {
                        read.chunks
                            .iter()
                            .filter(|c| !t[c.pos].is_empty())
                            .map(|c| aligner.global(&t[c.pos], &c.seq).score)
                            .sum::<i32>()
                    })
                    .collect();
                let probs: Vec<_> = scores
                    .iter()
                    .map(|x| {
                        scores
                            .iter()
                            .map(|y| ((y - x) as f64).exp())
                            .sum::<f64>()
                            .recip()
                    })
                    .collect();
                let new = *choises.choose_weighted(&mut rng, |&k| probs[k]).unwrap();
                let changed = (*asn != new) as usize;
                *asn = new;
                (*scores.iter().max().unwrap(), changed)
            })
            .fold((0, 0), |(x, y), (a, b)| (x + a, y + b));
        changes.push(changed);
        count = vec![0; k];
        for &asn in assignments.iter() {
            count[asn] += 1;
        }
        score = current_score;
        if changed < data.len() / 10 {
            break;
        }
    }
    trace!("HISTORY\t{:?}", changes);
    trace!("LK\t{}\t{:?}", score, count);
    assignments
}

fn take_consensus<R: Rng>(chunks: &[Vec<&[u8]>], r: &mut R, c: &ClusteringConfig) -> Vec<Vec<u8>> {
    chunks
        .iter()
        .map(|cs| match cs.is_empty() {
            true => vec![],
            false => consensus(cs, r, c),
        })
        .collect()
}

fn consensus<R: Rng>(cs: &[&[u8]], r: &mut R, c: &ClusteringConfig) -> Vec<u8> {
    let mut cs = cs.to_vec();
    let take_num = cs.len().min(c.consensus_size);
    let cs: Vec<_> = (0..c.consensus_number)
        .map(|_| {
            cs.shuffle(r);
            poa_hmm::POA::from_slice_default(&cs[..take_num]).consensus()
        })
        .collect();
    poa_hmm::POA::from_vec_default(&cs).consensus()
}
