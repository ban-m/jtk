use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
const SAMPLE_RATE: f64 = 0.02;
const INITIAL_BETA: f64 = 0.001;
const BETA_INCREASE: f64 = 1.02;
const MAX_BETA: f64 = 0.8;
const REPEAT_NUM: usize = 1;
const GIBBS_PRIOR: f64 = 0.02;
const STABLE_LIMIT: u32 = 8;
const VARIANT_FRACTION: f64 = 0.1;
#[derive(Debug, Clone)]
pub struct ClusteringConfig<F: Fn(u8, u8) -> i32> {
    pub threads: usize,
    pub cluster_num: usize,
    pub subchunk_length: usize,
    pub limit: u64,
    pub sample_rate: f64,
    pub alnparam: AlignmentParameters<F>,
    pub poa_config: poa_hmm::Config,
    pub seed: u64,
    pub id: u64,
    pub beta_increase: f64,
    pub initial_beta: f64,
    pub max_beta: f64,
    pub repeat_num: usize,
    pub gibbs_prior: f64,
    pub stable_limit: u32,
    pub variant_fraction: f64,
}

impl ClusteringConfig<fn(u8, u8) -> i32> {
    pub fn with_default(
        threads: usize,
        cluster_num: usize,
        subchunk_length: usize,
        limit: u64,
    ) -> Self {
        let seed = ((threads << cluster_num) * subchunk_length / 17) as u64 * limit;
        let id: u64 = thread_rng().gen::<u64>() % 100_000;
        Self {
            threads,
            cluster_num,
            subchunk_length,
            limit,
            alnparam: DEFAULT_ALN,
            poa_config: poa_hmm::DEFAULT_CONFIG,
            sample_rate: SAMPLE_RATE,
            seed,
            id,
            beta_increase: BETA_INCREASE,
            initial_beta: INITIAL_BETA,
            max_beta: MAX_BETA,
            repeat_num: REPEAT_NUM,
            gibbs_prior: GIBBS_PRIOR,
            stable_limit: STABLE_LIMIT,
            variant_fraction: VARIANT_FRACTION,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentParameters<F>
where
    F: Fn(u8, u8) -> i32,
{
    pub ins: i32,
    pub del: i32,
    pub score: F,
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        6
    } else {
        -12
    }
}

pub const DEFAULT_ALN: AlignmentParameters<fn(u8, u8) -> i32> = AlignmentParameters {
    ins: -10,
    del: -10,
    score,
};
