use super::ChunkedUnit;
use super::ClusteringConfig;
use poa_hmm::*;
use rand::seq::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use rayon::prelude::*;
use std::collections::VecDeque;
const CAP: usize = 50;
const INIT: usize = 30;
#[derive(Debug)]
struct Model {
    cluster_num: usize,
    chain_len: usize,
    weights: Vec<usize>,
    // Cluster -> Position -> Sequence.
    sequences: Vec<Vec<VecDeque<Vec<u8>>>>,
}

impl Model {
    fn initial_model<R: rand::Rng, F: Fn(u8, u8) -> i32>(
        data: &[ChunkedUnit],
        rng: &mut R,
        chain_len: usize,
        c: &ClusteringConfig<F>,
    ) -> Self {
        let cluster_num = c.cluster_num;
        let weights = vec![CAP; cluster_num];
        let sequences: Vec<_> = (0..cluster_num)
            .map(|_| {
                let mut chains = vec![VecDeque::with_capacity(CAP + 1); chain_len];
                for read in (0..INIT).filter_map(|_| data.choose(rng)) {
                    for chunk in read.chunks.iter() {
                        chains[chunk.pos].push_back(chunk.seq.to_vec());
                    }
                }
                chains
            })
            .collect();
        Self {
            cluster_num,
            weights,
            sequences,
            chain_len,
        }
    }
    fn push_seq(&mut self, read: &ChunkedUnit) {
        self.weights[read.cluster] += 1;
        for chunk in read.chunks.iter() {
            self.sequences[read.cluster][chunk.pos].push_back(chunk.seq.to_vec());
            while self.sequences[read.cluster][chunk.pos].len() > CAP {
                self.sequences[read.cluster][chunk.pos].pop_front();
            }
        }
    }
    fn generate_poa<R: Rng, F: Fn(u8, u8) -> i32 + std::marker::Sync>(
        &self,
        rng: &mut R,
        c: &ClusteringConfig<F>,
    ) -> Vec<Vec<POA>> {
        let seeds: Vec<Vec<_>> = self
            .sequences
            .iter()
            .map(|sequences| (0..sequences.len()).map(|_| rng.gen::<u64>()).collect())
            .collect();
        self.sequences
            .par_iter()
            .zip(seeds.par_iter())
            .map(|(sequences, seeds)| {
                sequences
                    .par_iter()
                    .zip(seeds.par_iter())
                    .map(|(cs, s)| {
                        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(*s);
                        let mut cs: Vec<&[u8]> = cs.iter().map(|x| x.as_slice()).collect();
                        cs.shuffle(&mut rng);
                        // super::create_model::construct_poa(&None, &cs, c)
                        let &super::AlignmentParameters {
                            ins,
                            del,
                            ref score,
                        } = &c.alnparam;
                        let param = (ins, del, score);
                        POA::default().update(&cs, &vec![1.; cs.len()], param)
                    })
                    .collect()
            })
            .collect()
    }
}

fn update_assignment<R: rand::Rng, F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    read: &mut ChunkedUnit,
    models: &Model,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) {
    let mut weights = get_weights(read, &models, rng, c);
    let sum = weights.iter().sum::<f64>();
    weights.iter_mut().for_each(|x| *x /= sum);
    trace!("{:?}", weights);
    let (new_assignment, _max) = weights
        .iter()
        .enumerate()
        .max_by(|a, b| (a.1).partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or((0, &0.));
    read.cluster = new_assignment;
}

fn get_weights<R: rand::Rng, F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    read: &mut ChunkedUnit,
    models: &Model,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> Vec<f64> {
    let poas = models.generate_poa(rng, c);
    let likelihoods: Vec<(usize, Vec<_>)> = read
        .chunks
        .iter()
        .map(|chunk| {
            let lks: Vec<_> = poas
                .iter()
                .map(|ms| ms[chunk.pos].forward(&chunk.seq, &c.poa_config))
                .collect();
            (chunk.pos, lks)
        })
        .collect();
    (0..c.cluster_num)
        .map(|l| {
            (0..c.cluster_num)
                .map(|k| {
                    if k == l {
                        return 0.;
                    }
                    // let (i, j) = (l.max(k), l.min(k));
                    let prior = (models.weights[k] as f64 / models.weights[l] as f64).ln();
                    prior
                        + likelihoods
                            .iter()
                            .map(|&(_, ref lks)| (lks[k] - lks[l]))
                            .sum::<f64>()
                })
                .map(|lkdiff| lkdiff.exp())
                .sum::<f64>()
                .recip()
        })
        .collect()
}

pub fn clustering_by_kmeans<F: Fn(u8, u8) -> i32 + std::marker::Sync>(
    data: &mut Vec<ChunkedUnit>,
    chain_len: usize,
    c: &ClusteringConfig<F>,
    ref_unit: u64,
) -> f64 {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(ref_unit);
    let rng = &mut rng;
    data.iter_mut()
        .for_each(|x| x.cluster = rng.gen::<usize>() % c.cluster_num);
    let mut models = Model::initial_model(data, rng, chain_len, c);
    let mut count = 0;
    // let mut margin = 0.;
    //    while margin < 3. * data.len() as f64 || count == 0 {
    while count < 10 {
        for picked in rand::seq::index::sample(rng, data.len(), data.len()).into_iter() {
            trace!("Pick {}, {:?}", picked, models.weights);
            update_assignment(&mut data[picked], &models, rng, c);
            models.push_seq(&data[picked]);
        }
        let (_, _, _, margin) =
            super::variant_calling::get_variants(&data, chain_len, rng, c, c.variant_num);
        count += 1;
        trace!("Margin:{}", margin);
    }
    0.
    // margin
}
