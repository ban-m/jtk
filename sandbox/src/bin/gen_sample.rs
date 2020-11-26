use poa_hmm::gen_sample;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
fn main() {
    let cluster_num = 2;
    let len = 100;
    let seed = 923;
    let num = 30;
    let p = gen_sample::PROFILE;
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let template = gen_sample::generate_seq(&mut rng, len);
    let templates: Vec<_> = (0..cluster_num)
        .map(|_| gen_sample::introduce_errors(&template, &mut rng, 1, 0, 0))
        .collect();
    let reads: Vec<_> = templates
        .iter()
        .enumerate()
        .flat_map(|(idx, t)| {
            (0..num)
                .map(|_| (idx, gen_sample::introduce_randomness(t, &mut rng, &p)))
                .collect::<Vec<_>>()
        })
        .collect();
    for (idx, read) in reads.iter().enumerate() {
        let line: String = read.1.iter().map(|&x| x as char).collect();
        println!(">{}\n{}", idx, line);
    }
}
