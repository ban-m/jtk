use haplotyper::ClusteringConfig;
use poa_hmm::*;
use rand_xoshiro::{Xoroshiro128PlusPlus, Xoshiro256StarStar};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let mut c = ClusteringConfig::default();
    c.cluster_num = 2;
    c.stable_limit = 6;
    c.variant_num = 2;
    c.read_type = haplotyper::ReadType::CLR;
    c.poa_config = poa_hmm::DEFAULT_CONFIG;
    let profile = gen_sample::PROFILE.norm().mul(0.20);
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
    //let mut templates = vec![];
    // let p = gen_sample::Profile {
    //     sub: errors / 6.,
    //     ins: errors / 6.,
    //     del: errors / 6.,
    // };
    let mut templates = vec![template.clone()];
    for _ in 0..clusters - 1 {
        use log::debug;
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
    let (mut dataset, _answer) =
        generate_mul_data(&templates, test_num, &mut rng, &probs, &profile);
    let ref_unit = 12232941;
    use rand::Rng;
    let mut rng: Xoshiro256StarStar = rand::SeedableRng::seed_from_u64(21 * ref_unit);
    dataset.iter_mut().for_each(|cs| {
        cs.cluster = rng.gen_range(0, c.cluster_num);
    });
    let dim = (2, chain_len);
    let (betas, _pos, _, _) =
        haplotyper::variant_calling::get_variants(&dataset, dim, &mut rng, &c, c.variant_num);
    for (i, bss) in betas.iter().enumerate() {
        for (j, bs) in bss.iter().enumerate() {
            eprintln!("{:?}", bs);
            let bs: Vec<_> = bs
                .iter()
                .enumerate()
                .filter(|&(_, &e)| e > 0.001)
                .map(|(i, e)| format!("{}:{:.3}", i, e))
                .collect();
            eprintln!("{}\t{}\t{}", i, j, bs.join(","));
        }
    }
    Ok(())
}
