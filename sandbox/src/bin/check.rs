use haplotyper::ClusteringConfig;
use log::*;
use poa_hmm::*;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.poa_config = poa_hmm::DEFAULT_CONFIG;
    c.read_type = definitions::ReadType::CLR;
    c.limit = 60;
    let profile = gen_sample::PROFILE.norm().mul(0.15);
    let args: Vec<_> = std::env::args().collect();
    let (seed, test_num, clusters, _errors, probs) = {
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
    c.cluster_num = clusters;
    c.variant_num = 2;
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed as u64);
    let chain_len = 20;
    let len = 100;
    let template: Vec<_> = (0..chain_len)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect::<Vec<_>>();
    // let p = gen_sample::Profile {
    //     sub: errors / 3.,
    //     ins: errors / 3.,
    //     del: errors / 3.,
    // };
    let mut templates = vec![template.clone()];
    assert!(clusters > 1);
    for _ in 0..clusters - 1 {
        let var_pos = rng.gen_range(0..chain_len);
        let mut seq = template.clone();
        seq[var_pos] = match rng.gen::<u8>() % 3 {
            0 => {
                debug!("Ins,{}", var_pos);
                gen_sample::introduce_errors(&seq[var_pos], &mut rng, 0, 0, 1)
            }
            1 => {
                debug!("Del,{}", var_pos);
                gen_sample::introduce_errors(&seq[var_pos], &mut rng, 0, 1, 0)
            }
            2 => {
                debug!("Subs,{}", var_pos);
                gen_sample::introduce_errors(&seq[var_pos], &mut rng, 1, 0, 0)
            }
            _ => panic!(),
        };
        // let seq: Vec<_> = template
        //     .iter()
        //     .map(|e| gen_sample::introduce_randomness(e, &mut rng, &p))
        //     .collect();
        templates.push(seq);
    }
    use sandbox::generate_mul_data;
    let (mut dataset, answer) = generate_mul_data(&templates, test_num, &mut rng, &probs, &profile);
    dataset
        .iter_mut()
        .zip(answer.iter())
        .for_each(|(x, &ans)| x.cluster = ans as usize);
    // let unit = definitions::Unit {
    //     id: 0,
    //     seq: String::new(),
    //     cluster_num: 2,
    // };
    // let start = std::time::Instant::now();
    // haplotyper::clustering_by_kmeans(&mut dataset, chain_len, &c, &unit, 10);
    // let end = std::time::Instant::now();
    // let preds: Vec<_> = dataset.iter().map(|x| x.cluster as u8).collect();
    // let score = haplotyper::rand_index(&preds, &answer);
    // let time = (end - start).as_millis();
    // println!("RESULT\t{}\tOLD\t{}\t{}", seed, score, time);
    let reads: Vec<Vec<_>> = dataset
        .iter()
        .map(|x| x.chunks.iter().flat_map(|x| x.seq.clone()).collect())
        .collect();
    let start = std::time::Instant::now();
    let config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(50, 30, 3, 2, 20);
    // let preds =
    //     haplotyper::local_clustering::kmeans::clustering_rep(&reads, &mut rng, &config).unwrap();
    let preds = haplotyper::local_clustering::kmeans::clustering_varcall(&reads, &mut rng, &config)
        .unwrap();
    let end = std::time::Instant::now();
    let score = haplotyper::rand_index(&preds, &answer);
    let time = (end - start).as_millis();
    println!("RESULT\t{}\tNEW\t{}\t{}", seed, score, time);
    use std::collections::HashMap;
    let mut result: HashMap<_, u32> = HashMap::new();
    let templates: Vec<_> = templates
        .iter()
        .map(|xs| {
            xs.iter()
                .flat_map(std::convert::identity)
                .copied()
                .collect::<Vec<_>>()
        })
        .collect();
    let mut sum = 0;
    for (idx, ((pred, answer), read)) in preds.iter().zip(answer).zip(&reads).enumerate() {
        let ed0 = edlib_sys::global_dist(&templates[0], &read) as i32;
        let diff: Vec<_> = templates
            .iter()
            .skip(1)
            .map(|template| {
                let diff = edlib_sys::global_dist(read, template) as i32 - ed0;
                sum += diff;
                diff
            })
            .map(|e| format!("{}", e))
            .collect();
        let diff = diff.join("\t");
        eprintln!("D\t{}\t{}\t{}\t{}", idx, answer, pred, diff);
        *result.entry((pred, answer)).or_default() += 1;
    }
    eprintln!("{}", sum as f64 / reads.len() as f64);
    let tempdiff = kiley::gen_seq::introduce_errors(&templates[0], &mut rng, 1, 0, 0);
    let mut sum = 0;
    for (idx, read) in reads.iter().enumerate() {
        let ed0 = edlib_sys::global_dist(&templates[0], &read) as i32;
        let ed1 = edlib_sys::global_dist(&tempdiff, &read) as i32;
        sum += (ed1 - ed0);
        eprintln!("{}\t{}", idx, ed1 - ed0);
    }
    eprintln!("{}", sum as f64 / reads.len() as f64);
    let mut result: Vec<_> = result.into_iter().collect();
    result.sort_by_key(|x| x.1);
    for ((pred, ans), count) in result {
        eprintln!("{}\t{}\t{}", pred, ans, count);
    }
    Ok(())
}
