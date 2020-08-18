use haplotyper::clustering_by_kmeans;
use haplotyper::ClusteringConfig;
use poa_hmm::*;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.initial_beta = 0.0001;
    c.max_beta = 0.3;
    c.cluster_num = 2;
    c.beta_increase = 1.03;
    c.stable_limit = 6;
    c.repeat_num = 3;
    c.read_type = haplotyper::ReadType::CLR;
    c.poa_config = poa_hmm::DEFAULT_CONFIG;
    let profile = gen_sample::PROFILE.norm().mul(0.20);
    //let profile = gen_sample::CCS_PROFILE;
    let args: Vec<_> = std::env::args().collect();
    let (seed, test_num, clusters, errors, probs) = {
        let seed: usize = args[1].parse().unwrap();
        let test_num: usize = args[2].parse().unwrap();
        let clusters: usize = args[3].parse().unwrap();
        let errors: f64 = args[4].parse().unwrap();
        let mut probs: Vec<f64> = args[5..].iter().map(|e| e.parse().unwrap()).collect();
        let sum = probs.iter().sum::<f64>();
        probs.iter_mut().for_each(|e| *e /= sum);
        assert_eq!(clusters, probs.len());
        (seed, test_num, clusters, errors, probs)
    };
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let chain_len = 20;
    let mut templates = vec![];
    let len = 100;
    let p = &gen_sample::Profile {
        sub: errors / 6.,
        ins: errors / 6.,
        del: errors / 6.,
    };
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    for _ in 0..clusters {
        let seq: Vec<_> = template
            .iter()
            .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
            .collect();
        templates.push(seq);
    }
    use sandbox::generate_mul_data;
    let (dataset, answer) = generate_mul_data(&templates, test_num, &mut rng, &probs, &profile);
    let p = gen_sample::Profile {
        sub: 0.004,
        del: 0.004,
        ins: 0.004,
    };
    let mut dataset = (0..4 * dataset.len())
        .map(|i| {
            if i < dataset.len() {
                dataset[i].clone()
            } else {
                let mut d = dataset[i % dataset.len()].clone();
                d.chunks.iter_mut().for_each(|chunk| {
                    chunk.seq = gen_sample::introduce_randomness(&chunk.seq, &mut rng, &p);
                });
                d
            }
        })
        .collect();
    clustering_by_kmeans(&mut dataset, chain_len, &c, 0);
    use std::collections::HashMap;
    let mut result: HashMap<_, u32> = HashMap::new();
    for (data, answer) in dataset.iter().zip(answer) {
        *result.entry((data.cluster, answer)).or_default() += 1;
    }
    let mut result: Vec<_> = result.into_iter().collect();
    result.sort_by_key(|x| x.1);
    for ((pred, ans), count) in result {
        eprintln!("{}\t{}\t{}", pred, ans, count);
    }
    Ok(())
}
