use definitions::*;
use kiley::gen_seq;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(std::fs::File::open(&args[1]).map(BufReader::new)?).unwrap();
    let hmm = haplotyper::model_tune::get_model(&ds).unwrap();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(394820394);
    let seq1 = gen_seq::generate_seq(&mut rng, 100);
    let seq2 = gen_seq::generate_seq(&mut rng, 100);
    let hap1 = vec![seq1.clone(), vec![b'A'; 4], seq2.clone()].concat();
    let hap2 = vec![seq1.clone(), vec![b'A'; 6], seq2.clone()].concat();
    let prof = gen_seq::PROFILE;
    for i in 0..100 {
        let read = match i % 2 == 0 {
            true => gen_seq::introduce_randomness(&hap1, &mut rng, &prof),
            false => gen_seq::introduce_randomness(&hap2, &mut rng, &prof),
        };
        let lk1 = hmm.likelihood(&hap1, &read, 20);
        let lk2 = hmm.likelihood(&hap2, &read, 20);
        println!("{i}\t{}\t{}", i % 2, lk1 - lk2);
    }
    Ok(())
}
