use haplotyper::clustering_by_kmeans;
use haplotyper::ClusteringConfig;
use poa_hmm::*;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.cluster_num = 2;
    c.variant_num = 2;
    c.poa_config = poa_hmm::DEFAULT_CONFIG;
    c.read_type = haplotyper::ReadType::CLR;
    let profile = gen_sample::PROFILE.norm().mul(0.15);
    // c.read_type = haplotyper::ReadType::CCS;
    // let profile = gen_sample::CCS_PROFILE;
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
        use log::debug;
        use rand::Rng;
        let var_pos = rng.gen_range(0, chain_len);
        let mut seq = template.clone();
        seq[var_pos] = match rng.gen::<u8>() % 3 {
            0 => {
                debug!("Ins");
                gen_sample::introduce_errors(&seq[var_pos], &mut rng, 0, 0, 1)
            }
            1 => {
                debug!("Del");
                gen_sample::introduce_errors(&seq[var_pos], &mut rng, 0, 1, 0)
            }
            2 => {
                debug!("Subs");
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
    let unit = definitions::Unit {
        id: 0,
        seq: String::new(),
        cluster_num: 2,
    };
    clustering_by_kmeans(&mut dataset, chain_len, &c, &unit, 10);
    // let asn: Vec<_> = dataset.iter().map(|c| c.cluster).collect();
    // (m, asn)
    //     })
    //     .max_by(|x, y| x.0.partial_cmp(&y.0).unwrap())
    //     .unwrap();
    // dataset.iter_mut().zip(asn).for_each(|(x, y)| x.cluster = y);
    use std::collections::HashMap;
    let mut result: HashMap<_, u32> = HashMap::new();
    for (idx, (data, answer)) in dataset.iter().zip(answer).enumerate() {
        let ed0 = data
            .chunks
            .iter()
            .zip(templates[0].iter())
            .map(|(r, q)| bio_utils::alignments::edit_dist(&r.seq, q) as i32)
            .sum::<i32>();
        let diff: Vec<_> = templates
            .iter()
            .map(|template| {
                data.chunks
                    .iter()
                    .zip(template.iter())
                    .map(|(r, q)| bio_utils::alignments::edit_dist(&r.seq, q) as i32)
                    .sum::<i32>()
                    - ed0
            })
            .map(|e| format!("{}", e))
            .collect();
        let diff = diff.join("\t");
        eprintln!("D\t{}\t{}\t{}\t{}", idx, answer, data.cluster, diff);
        *result.entry((data.cluster, answer)).or_default() += 1;
    }
    let mut result: Vec<_> = result.into_iter().collect();
    result.sort_by_key(|x| x.1);
    for ((pred, ans), count) in result {
        eprintln!("{}\t{}\t{}", pred, ans, count);
    }
    Ok(())
}
