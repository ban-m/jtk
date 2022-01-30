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
    // c.poa_config = poa_hmm::DEFAULT_CONFIG;
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
    let len = 2000;
    let template: Vec<_> = gen_seq::generate_seq(&mut rng, len);
    let mut templates = vec![template.clone()];
    for _ in 0..clusters - 1 {
        let mut seq = template.clone();
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
        };
        templates.push(seq);
    }
    use sandbox::generate_test_data;
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
        let clusters = clusters as u8;
        use haplotyper::local_clustering::kmeans::ClusteringConfig;
        let config = ClusteringConfig::new(100, clusters, coverage, definitions::ReadType::CLR);
        use haplotyper::local_clustering::kmeans;
        let (preds, gains, _, _) = kmeans::clustering(&dataset, &mut rng, &config).unwrap();
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
        for ((p, a), post) in preds.iter().zip(answer.iter()).zip(gains.iter()) {
            log::debug!("{}\t{}\t{}", p, a, vec2str(post));
        }
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

fn vec2str(xs: &[f64]) -> String {
    let xs: Vec<_> = xs.iter().map(|&x| format!("{:6.1}", x)).collect();
    xs.join(",")
}
