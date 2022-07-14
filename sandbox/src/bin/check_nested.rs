use haplotyper::local_clustering::ClusteringConfig;
use kiley::gen_seq;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.read_type = definitions::ReadType::CLR;
    c.limit = 60;
    let args: Vec<_> = std::env::args().collect();
    let seed: u64 = args[1].parse().unwrap();
    let test_num: usize = args[2].parse().unwrap();
    let profile = gen_seq::PROFILE.norm().mul(0.05);
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed);
    let len = 2000;
    let hap1: Vec<_> = gen_seq::generate_seq(&mut rng, len);
    let div_rate = gen_seq::Profile {
        sub: 0.05 / 3f64,
        ins: 0.05 / 3f64,
        del: 0.05 / 3f64,
    };
    let hap2 = gen_seq::introduce_randomness(&hap1, &mut rng, &div_rate);
    use sandbox::generate_test_data;
    let probs = vec![0.5, 0.5];
    let templates = vec![hap1.clone(), hap2];
    let (dataset, answer) = generate_test_data(&templates, test_num, &mut rng, &probs, &profile);
    let mut counts: HashMap<_, u32> = HashMap::new();
    for &ans in answer.iter() {
        *counts.entry(ans).or_default() += 1;
    }
    log::debug!("{}\t{}", seed, test_num);
    {
        let start = std::time::Instant::now();
        let coverage = (test_num / 4) as f64;
        let clusters = 4;
        use haplotyper::local_clustering::kmeans::ClusteringConfig;
        let config = ClusteringConfig::new(100, clusters, coverage, 1.8);
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
    for (i, (ans, read)) in answer.iter().zip(dataset.iter()).enumerate() {
        let dist = kiley::bialignment::edit_dist(&templates[0], read) as i32;
        let others: Vec<_> = templates
            .iter()
            .map(|x| kiley::bialignment::edit_dist(x, read) as i32 - dist)
            .collect();
        let others: Vec<_> = others.iter().map(|x| format!("{}", x)).collect();
        log::trace!("OPTIMAL\t{}\t{}\t{}", i, ans, others.join("\t"));
    }
    Ok(())
}
