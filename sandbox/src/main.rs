use poa_hmm::*;
use rand::{seq::SliceRandom, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use poa_hmm::gen_sample::*;
fn main() {
    let bases = b"ACTG";
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1212132);
    let template: Vec<_> = (0..150)
        .filter_map(|_| bases.choose(&mut rng))
        .copied()
        .collect();
    let model1: Vec<Vec<_>> = (0..50)
        .map(|_| introduce_randomness(&template, &mut rng, &CCS_PROFILE))
        .collect();
    let model1: Vec<&[u8]> = model1.iter().map(|e| e.as_slice()).collect();
    let ps = (-1, -1, &|x, y| if x == y { 1 } else { -1 });
    let model = POA::generate_banded(&model1, ps, 10, 101010);
    println!("{}", model);
}
