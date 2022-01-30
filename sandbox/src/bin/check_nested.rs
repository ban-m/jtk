use haplotyper::local_clustering::ClusteringConfig;
use kiley::gen_seq;
use log::*;
// use poa_hmm::*;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.read_type = definitions::ReadType::CLR;
    c.limit = 60;
    let args: Vec<_> = std::env::args().collect();
    let (seed, test_num, clusters, errors) = {
        let seed: usize = args[1].parse().unwrap();
        let test_num: usize = args[2].parse().unwrap();
        let clusters: usize = args[3].parse().unwrap();
        let errors: f64 = args[4].parse().unwrap();
        (seed, test_num, clusters, errors)
    };
    let profile = gen_seq::PROFILE.norm().mul(errors);
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let len = 2000;
    let template: Vec<_> = gen_seq::generate_seq(&mut rng, len);
    let mut templates = vec![template.clone()];
    fn intro_var<R: Rng>(seq: &[u8], rng: &mut R, count: usize) -> Vec<u8> {
        let mut seq = seq.to_vec();
        for _ in 0..count {
            let var_pos = rng.gen_range(0..seq.len());
            let var_type = loop {
                let var_type = rng.gen::<usize>() % 9;
                if 4 <= var_type || b"ACGT"[var_type] != seq[var_pos] {
                    break var_type;
                }
            };
            match var_type {
                0..=3 => {
                    debug!("INTR\tSubs,{}", var_pos);
                    seq[var_pos] = b"ACGT"[var_type];
                }
                4..=7 => {
                    debug!("INTR\tIns,{}", var_pos);
                    seq.insert(var_pos, b"ACGT"[var_type - 4]);
                }
                _ => {
                    debug!("INTR\tDel,{}", var_pos);
                    seq.remove(var_pos);
                }
            }
        }
        seq
    }
    for i in (0..clusters).rev() {
        let var_num = 2usize.pow(i as u32);
        templates = templates
            .iter()
            .flat_map(|template| {
                let a1 = intro_var(template, &mut rng, var_num);
                // let a2 = intro_var(template, &mut rng, var_num);
                vec![template.clone(), a1]
            })
            .collect();
    }
    let clusters = templates.len();
    use sandbox::generate_test_data;
    let probs = vec![(clusters as f64).recip(); clusters];
    let (dataset, answer) = generate_test_data(&templates, test_num, &mut rng, &probs, &profile);
    let mut counts: HashMap<_, u32> = HashMap::new();
    for &ans in answer.iter() {
        *counts.entry(ans).or_default() += 1;
    }
    for k in 0..clusters {
        debug!("CL\t{}\t{}", k, counts[&(k as u8)]);
    }
    {
        let start = std::time::Instant::now();
        let coverage = (dataset.len() / clusters) as f64;
        let clusters = templates.len() as u8;
        use haplotyper::local_clustering::kmeans::ClusteringConfig;
        let config = ClusteringConfig::new(100, clusters, coverage, definitions::ReadType::CLR);
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
