use bio_utils::alignments::edit_dist;
use poa_hmm::*;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
const REFERENCE_LENGTH: usize = 1_000;
const QUERY_NUM: usize = 200;
fn main() {
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(239180);
    let reference = gen_sample::generate_seq(&mut rng, REFERENCE_LENGTH);
    let mut haplotypes = vec![reference.clone()];
    haplotypes.push(gen_sample::introduce_errors(&reference, &mut rng, 0, 0, 1));
    let profile = &gen_sample::PROFILE;
    println!("ID\tAnswer\tDist1\tDist2");
    for i in 0..QUERY_NUM {
        let answer = 2 * i / QUERY_NUM;
        let query = gen_sample::introduce_randomness(&haplotypes[answer], &mut rng, &profile);
        let dist1 = edit_dist(&haplotypes[0], &query);
        let dist2 = edit_dist(&haplotypes[1], &query);
        println!("{}\t{}\t{}\t{}", i, answer, dist1, dist2);
    }
}
