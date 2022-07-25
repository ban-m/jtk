use haplotyper::local_clustering::ClusteringConfig;
use kiley::gen_seq;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::collections::HashMap;
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
    assert_eq!(clusters, 2);
    let profile = gen_seq::PROFILE.norm().mul(errors);
    c.cluster_num = clusters;
    c.variant_num = 2;
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let len = 2000;
    let template: Vec<_> = gen_seq::generate_seq(&mut rng, len);
    let mut templates = vec![template.clone()];
    let (var_type, seq) = {
        let mut seq = template;
        let var_pos = rng.gen_range(0..seq.len());
        let var_type = loop {
            let var_type = rng.gen::<usize>() % 9;
            if 4 <= var_type || b"ACGT"[var_type] != seq[var_pos] {
                break var_type;
            }
        };
        log::debug!("POS\t{}\t{}", var_pos, var_type);
        match var_type {
            0..=3 => {
                seq[var_pos] = b"ACGT"[var_type];
                ("Sub", seq)
            }
            4..=7 => {
                seq.insert(var_pos, b"ACGT"[var_type - 4]);
                ("Ins", seq)
            }
            _ => {
                seq.remove(var_pos);
                ("Del", seq)
            }
        }
    };
    templates.push(seq);
    use sandbox::generate_test_data;
    let (dataset, answer) = generate_test_data(&templates, test_num, &mut rng, &probs, &profile);
    let mut counts: HashMap<_, u32> = HashMap::new();
    for &ans in answer.iter() {
        *counts.entry(ans).or_default() += 1;
    }
    {
        let start = std::time::Instant::now();
        use haplotyper::local_clustering::kmeans;
        let (preds, _, _, _) = kmeans::clustering(&dataset, &mut rng, clusters, 100).unwrap();
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
        println!(
            "RESULT\t{}\tNEW\t{}\t{}\t{}\t{}",
            seed, score, time, acc, var_type
        );
        log::debug!("\n{:?}\n{:?}", preds, answer);
    }
    Ok(())
}
