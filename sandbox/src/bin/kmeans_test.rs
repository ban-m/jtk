use kiley::alignment::bialignment::edit_dist;
use log::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let len = 2000;
    let num_hap = 2;
    let coverage = 20;
    let config = haplotyper::local_clustering::kmeans::ClusteringConfig::new(30, 20, 10, 2, 30);
    for seed in 0..10u64 {
        let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(seed);
        let template = kiley::gen_seq::generate_seq(&mut rng, len);
        let alt_hap = match rng.gen::<u8>() % 3 {
            0 => kiley::gen_seq::introduce_errors(&template, &mut rng, 0, 1, 0),
            1 => kiley::gen_seq::introduce_errors(&template, &mut rng, 1, 0, 0),
            _ => kiley::gen_seq::introduce_errors(&template, &mut rng, 0, 0, 1),
        };
        let haps = vec![template.clone(), alt_hap];
        for (i, h) in haps.iter().enumerate() {
            for (j, l) in haps.iter().take(i).enumerate() {
                let d = kiley::alignment::bialignment::edit_dist(h, l);
                println!("{}\t{}\t{}", i, j, d);
            }
        }
        let p = kiley::gen_seq::PROFILE;
        let (reads, answers): (Vec<_>, Vec<_>) = (0..coverage * num_hap)
            .map(|x| {
                let answer = (x / coverage) as u8;
                let read =
                    kiley::gen_seq::introduce_randomness(&haps[answer as usize], &mut rng, &p);
                (read, answer)
            })
            .unzip();
        let assignments =
            haplotyper::local_clustering::kmeans::clustering(&reads, &mut rng, &config);
        let score = haplotyper::rand_index(&answers, &assignments);
        for ((r, a), p) in reads.iter().zip(answers.iter()).zip(assignments.iter()) {
            let true_dist = edit_dist(r, &haps[0]) as i32 - edit_dist(r, &haps[1]) as i32;
            println!("{}\t{}\t{}", a, true_dist, p);
        }
        println!("{}\n", score);
    }
    Ok(())
}
