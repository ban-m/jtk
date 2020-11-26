use log::*;
use poa_hmm::gen_sample;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
const CLUSTER_NUM: usize = 2;
fn main() {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("trace")).init();
    debug!("Start");
    let len = 100;
    let chain = 20;
    let seed = 923;
    let num = 20;
    let div = gen_sample::Profile {
        sub: 0.001 / 3.,
        ins: 0.001 / 3.,
        del: 0.001 / 3.,
    };
    let p = gen_sample::PROFILE;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let template: Vec<_> = (0..chain)
        .map(|_| gen_sample::generate_seq(&mut rng, len))
        .collect();
    let templates: Vec<Vec<_>> = (0..CLUSTER_NUM)
        .map(|_| {
            template
                .iter()
                .map(|t| gen_sample::introduce_randomness(t, &mut rng, &div))
                .collect()
        })
        .collect();
    for (i, t) in templates.iter().enumerate() {
        for (j, s) in templates.iter().enumerate().skip(i + 1) {
            let dist = t
                .iter()
                .zip(s)
                .map(|(a, b)| bio::alignment::distance::levenshtein(a, b))
                .sum::<u32>();
            eprintln!("DIST\t{}\t{}\t{}", i, j, dist);
        }
    }
    use sandbox::clustering;
    let reads = gen_reads(&templates, num, &mut rng, &p);
    let assignments = clustering(&reads, chain, CLUSTER_NUM, 21);
    println!("ID\tANSWER\tPRED");
    for (idx, asn) in assignments.iter().enumerate() {
        println!("{}\t{}\t{}", idx, idx / num, asn);
    }
}

fn gen_reads<R: rand::Rng>(
    templates: &[Vec<Vec<u8>>],
    num: usize,
    rng: &mut R,
    p: &gen_sample::Profile,
) -> Vec<(usize, Vec<Vec<u8>>)> {
    templates
        .iter()
        .enumerate()
        .flat_map(|(idx, t)| {
            (0..num)
                .map(|_| {
                    let seq: Vec<_> = t
                        .iter()
                        .map(|s| gen_sample::introduce_randomness(s, rng, p))
                        .collect();
                    (idx, seq)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}
