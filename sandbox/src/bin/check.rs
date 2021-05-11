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
        templates.push(seq);
    }
    use sandbox::generate_mul_data;
    let (mut dataset, answer) = generate_mul_data(&templates, test_num, &mut rng, &probs, &profile);
    // {
    // let unit = definitions::Unit {
    //     id: 0,
    //     seq: String::new(),
    //     cluster_num: clusters,
    // };
    //     let start = std::time::Instant::now();
    //     haplotyper::clustering_by_kmeans(&mut dataset, chain_len, &c, &unit, 10);
    //     let end = std::time::Instant::now();
    //     let preds: Vec<_> = dataset.iter().map(|x| x.cluster as u8).collect();
    //     let score = haplotyper::rand_index(&preds, &answer);
    //     let time = (end - start).as_millis();
    // let acc = preds
    //     .iter()
    //     .zip(answer.iter())
    //     .filter(|x, y| x == y)
    //     .count() as f64
    //     / preds.len() as f64;
    // let acc = (1f64 - acc).max(acc);
    //     println!("RESULT\t{}\tOLD\t{}\t{}\t{}", seed, score, time, acc);
    // }
    {
        let reads: Vec<Vec<_>> = dataset
            .iter()
            .map(|x| x.chunks.iter().flat_map(|x| x.seq.clone()).collect())
            .collect();
        let start = std::time::Instant::now();
        let clusters = clusters as u8;
        let config =
            haplotyper::local_clustering::kmeans::ClusteringConfig::new(100, 30, 3, clusters, 20);
        let (preds, _) =
            haplotyper::local_clustering::kmeans::clustering(&reads, &mut rng, &config).unwrap();
        let end = std::time::Instant::now();
        let score = haplotyper::rand_index(&preds, &answer);
        let time = (end - start).as_millis();
        let acc = preds
            .iter()
            .zip(answer.iter())
            .filter(|(x, y)| x == y)
            .count() as f64
            / preds.len() as f64;
        let acc = (1f64 - acc).max(acc);
        println!("RESULT\t{}\tNEW\t{}\t{}\t{}", seed, score, time, acc);
        // eprintln!("{:?}\n{:?}", preds, answer);
    }
    Ok(())
}
