use haplotyper::local_clustering::ClusteringConfig;
use kiley::gen_seq;
use log::*;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
const LK: f64 = 1.8;

fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.read_type = definitions::ReadType::CLR;
    c.limit = 60;
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
    let profile = gen_seq::PROFILE.norm().mul(errors);
    c.cluster_num = clusters;
    c.variant_num = 2;
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let len = 1000;
    let former_half = gen_seq::generate_seq(&mut rng, len);
    let latter_half = gen_seq::generate_seq(&mut rng, len);
    let rep: &[u8] = if rng.gen_bool(0.5) { b"AT" } else { b"A" };
    let mut templates = vec![];
    for i in 0..clusters {
        let rep_range = 10usize..20usize;
        let rep_num = rng.gen_range(rep_range.clone());
        debug!("{}\t{}", i, rep_num);
        let mut template = former_half.clone();
        for _ in 0..rep_num {
            template.extend_from_slice(rep);
        }
        template.extend_from_slice(&latter_half);
        templates.push(template);
    }
    use sandbox::generate_test_data;
    let (dataset, answer) = generate_test_data(&templates, test_num, &mut rng, &probs, &profile);
    let mut counts: HashMap<_, u32> = HashMap::new();
    for &ans in answer.iter() {
        *counts.entry(ans).or_default() += 1;
    }
    for k in 0..clusters {
        debug!("CL\t{}\t{}", k, counts[&k]);
    }
    {
        let start = std::time::Instant::now();
        let coverage = (dataset.len() / clusters) as f64;
        let clusters = clusters as u8;
        use haplotyper::local_clustering::kmeans::ClusteringConfig;
        let config = ClusteringConfig::new(100, clusters, coverage, LK, definitions::ReadType::CLR);
        use haplotyper::local_clustering::kmeans;
        let (preds, _, _, _) = kmeans::clustering(&dataset, &mut rng, &config).unwrap();
        let end = std::time::Instant::now();
        let score = haplotyper::local_clustering::rand_index(&preds, &answer);
        let time = (end - start).as_millis();
        let acc = preds
            .iter()
            .zip(answer.iter())
            .filter(|(x, y)| x == y)
            .count() as f64
            / preds.len() as f64;
        let acc = (1f64 - acc).max(acc);
        println!("RESULT\t{}\tNEW\t{}\t{}\t{}", seed, score, time, acc);
        log::debug!("\n{:?}\n{:?}", preds, answer);
    }
    Ok(())
}
