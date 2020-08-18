use poa_hmm::*;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let seed: u64 = args[1].parse().unwrap();
    let cov: usize = args[2].parse().unwrap();
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(seed);
    let template = gen_sample::generate_seq(&mut rng, 100);
    let p = gen_sample::PROFILE.norm().mul(0.2);
    let seqs: Vec<_> = (0..cov)
        .map(|_| gen_sample::introduce_randomness(&template, &mut rng, &p))
        .collect();
    let seqs: Vec<_> = seqs.iter().map(|x| x.as_slice()).collect();
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let node_num_thr = (max_len as f64 * 1.5).floor() as usize;
    let param = (-1, -1, &score);
    let model0 = POA::from_slice(&seqs, &vec![1.; cov], param);
    let model = seqs
        .into_iter()
        .fold(POA::default(), |x, y| {
            let res = if x.nodes().len() > node_num_thr {
                x.add(y, 1., param).remove_node(0.5)
            } else {
                x.add(y, 1., param)
            };
            res
        })
        .remove_node(0.7)
        .finalize();
    println!("DEFAULT:{}", model0);
    println!("TUNED:{}\n", model);
}

fn score(x: u8, y: u8) -> i32 {
    if x == y {
        1
    } else {
        -1
    }
}
