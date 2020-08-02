use poa_hmm::POA;
use poa_hmm::*;
use rand::seq::*;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
fn super_long_consensus_rand() {
    use rayon::prelude::*;
    let bases = b"ACTG";
    let coverage = 100usize;
    let start = 20usize;
    let len = 2000;
    let result: Vec<_> = (start..coverage)
        .into_par_iter()
        .map(|cov| {
            let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(1_234_567);
            let template1: Vec<_> = (0..len)
                .filter_map(|_| bases.choose(&mut rng))
                .copied()
                .collect();
            let seqs: Vec<_> = (0..cov)
                .map(|_| {
                    gen_sample::introduce_randomness(&template1, &mut rng, &gen_sample::PROFILE)
                })
                .collect();
            let normal_cons = POA::from_vec_default(&seqs).consensus();
            let dist = edit_dist(&normal_cons, &template1);
            let consensus = consensus(&seqs);
            let dist2 = edit_dist(&consensus, &template1);
            eprintln!("{}\t{}", dist, dist2);
            (dist, dist2)
        })
        .collect();
    let normal = result.iter().map(|x| x.0).sum::<u32>();
    let multi = result.iter().map(|x| x.1).sum::<u32>();
    eprintln!("{}\t{}", normal, multi);
}

fn edit_dist(x1: &[u8], x2: &[u8]) -> u32 {
    let mut dp = vec![vec![0; x2.len() + 1]; x1.len() + 1];
    for (i, row) in dp.iter_mut().enumerate() {
        row[0] = i as u32;
    }
    for j in 0..=x2.len() {
        dp[0][j] = j as u32;
    }
    for (i, x1_b) in x1.iter().enumerate() {
        for (j, x2_b) in x2.iter().enumerate() {
            let m = if x1_b == x2_b { 0 } else { 1 };
            dp[i + 1][j + 1] = (dp[i][j + 1] + 1).min(dp[i + 1][j] + 1).min(dp[i][j] + m);
        }
    }
    dp[x1.len()][x2.len()]
}

fn consensus(seqs: &[Vec<u8>]) -> Vec<u8> {
    if seqs.len() <= 10 {
        POA::from_vec_default(&seqs).consensus()
    } else {
        let subchunks = 3 * seqs.len() / 10;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3343);
        let subseq: Vec<_> = (0..subchunks)
            .map(|_| {
                let subchunk: Vec<_> = seqs
                    .choose_multiple(&mut rng, 10)
                    .map(|e| e.as_slice())
                    .collect();
                POA::from_slice_default(&subchunk).consensus()
            })
            .collect();
        POA::from_vec_default(&subseq).consensus()
    }
}

fn main() {
    super_long_consensus_rand();
}
