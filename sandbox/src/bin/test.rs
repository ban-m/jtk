use poa_hmm::*;
use rand_xoshiro::Xoroshiro128PlusPlus;
use rayon::prelude::*;
fn main() {
    for j in 90..100 {
        use rand::SeedableRng;
        //let j = 0;
        // let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(12320);
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(2 * j as u64);
        let template_a = gen_sample::generate_seq(&mut rng, 100);
        let template_b = gen_sample::introduce_errors(&template_a, &mut rng, 1, 0, 0);
        eprintln!("{}", String::from_utf8_lossy(&template_a));
        eprintln!("{}", String::from_utf8_lossy(&template_b));
        let profile = gen_sample::PROFILE.norm().mul(0.2);
        let num = 20;
        let data_a: Vec<_> = (0..num)
            .map(|_| gen_sample::introduce_randomness(&template_a, &mut rng, &profile))
            .collect();
        let data_b: Vec<_> = (0..num)
            .map(|_| gen_sample::introduce_randomness(&template_b, &mut rng, &profile))
            .collect();
        (0..1000).into_par_iter().for_each(|i| {
            let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(i as u64 * 2);
            let template = if i < 500 { &template_a } else { &template_b };
            let query = gen_sample::introduce_randomness(template, &mut rng, &profile);
            let mut lks = vec![vec![]; 2];
            let rep_num = 10;
            for _ in 0..rep_num {
                use rand::seq::SliceRandom;
                let mut data_a: Vec<_> = data_a.clone();
                let mut data_b = data_b.clone();
                data_a.shuffle(&mut rng);
                data_b.shuffle(&mut rng);
                let m1 = POA::from_vec_default(&data_a);
                let m2 = POA::from_vec_default(&data_b);
                let lk1 = m1.forward(&query, &DEFAULT_CONFIG);
                let lk2 = m2.forward(&query, &DEFAULT_CONFIG);
                //println!("{}\t{}\t{}\t{}", i < 500, w1, w2, j);
                // assert!((1. - w1 - w2).abs() < 0.001);
                lks[0].push(lk1);
                lks[1].push(lk2);
            }
            let lks: Vec<_> = lks
                .iter()
                .map(|xs| logsumexp(xs) - (rep_num as f64).ln())
                .collect();
            // lks[0] /= rep_num as f64;
            // lks[1] /= rep_num as f64;
            let w1 = ((lks[1] - lks[0]).exp() + 1.).recip();
            let w2 = ((lks[0] - lks[1]).exp() + 1.).recip();
            // let ed = bio_utils::alignments::edit_dist(&query, &template_b) as i32
            //     - bio_utils::alignments::edit_dist(&query, &template_a) as i32;
            //println!("{}\t{}\t{}\t{}", i < 500, lks[0], lks[1], j);
            println!("{}\t{}\t{}\t{}", i < 500, w1, w2, j);
        });
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}
