use poa_hmm::*;
use rand::seq::SliceRandom;
use rand_xoshiro::Xoroshiro128PlusPlus;
use rayon::prelude::*;
fn main() {
    for j in 90..95 {
        use rand::SeedableRng;
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(2 * j as u64);
        let template_a = gen_sample::generate_seq(&mut rng, 100);
        //let template_b = gen_sample::introduce_errors(&template_a, &mut rng, 0, 0, 1);
        let template_b = gen_sample::introduce_errors(&template_a, &mut rng, 0, 0, 1);
        eprintln!("{}:{}", j, String::from_utf8_lossy(&template_a));
        eprintln!("{}:{}", j, String::from_utf8_lossy(&template_b));
        let profile = gen_sample::PROFILE.norm().mul(0.2);
        let num = 40;
        let data_a: Vec<_> = (0..num)
            .map(|_| gen_sample::introduce_randomness(&template_a, &mut rng, &profile))
            .collect();
        let data_b: Vec<_> = (0..num)
            .map(|_| gen_sample::introduce_randomness(&template_b, &mut rng, &profile))
            .collect();
        for i in 0..10 {
            let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(i as u64 * 2);
            let query = gen_sample::introduce_randomness(&template_a, &mut rng, &profile);
            let rep_num = 1000;
            use bio_utils::alignments::edit_dist;
            let ed1 = edit_dist(&query, &template_b) as i32 - edit_dist(&query, &template_a) as i32;
            let ed2 = edit_dist(&query, &template_a) as i32 - edit_dist(&query, &template_b) as i32;
            (0..rep_num).into_par_iter().for_each(|t| {
                let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(8 * t as u64);
                let data_a: Vec<_> = data_a.clone();
                let data_b: Vec<_> = data_b.clone();
                let m1 = model(data_a, &mut rng);
                let m2 = model(data_b, &mut rng);
                let lk1 = m1.forward(&query, &DEFAULT_CONFIG);
                let lk2 = m2.forward(&query, &DEFAULT_CONFIG);
                println!("{}\t{}\t{}\t{}\t{}", j, i, lk1, true, ed1);
                println!("{}\t{}\t{}\t{}\t{}", j, i, lk2, false, ed2);
            });
        }
    }
}
fn model<R: rand::Rng>(mut xs: Vec<Vec<u8>>, rng: &mut R) -> POA {
    // let mut xs: Vec<_> = (0..2 * xs.len())
    //     .map(|i| xs[i % xs.len()].as_slice())
    //     .collect();
    // xs.shuffle(rng);
    //POA::from_slice_default(&xs)
    let xs: Vec<_> = (0..10)
        .map(|_| {
            xs.shuffle(rng);
            POA::from_vec_default(&xs).consensus()
        })
        .collect();
    POA::from_vec_default(&xs)
}

// fn logsumexp(xs: &[f64]) -> f64 {
//     if xs.is_empty() {
//         return 0.;
//     }
//     let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
//     let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
//     assert!(sum >= 0., "{:?}->{}", xs, sum);
//     max + sum
// }
