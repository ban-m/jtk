const BASES: &[char] = &['A', 'C', 'G', 'T'];
use rand::seq::SliceRandom;
fn main() {
    use std::time::Instant;
    let mut rng = rand::thread_rng();
    let units: Vec<Vec<_>> = (0..10)
        .map(|_| {
            (0..2000)
                .filter_map(|_| BASES.choose(&mut rng))
                .copied()
                .collect()
        })
        .collect();
    let mut score = 0;
    let start = Instant::now();
    let k: u32 = 7;
    let kmer_vec: Vec<_> = units
        .iter()
        .map(|unit| {
            let mut fingerprint = vec![false; 4usize.pow(k)];
            for kmer in unit.windows(k as usize) {
                fingerprint[to_index(kmer)] = true;
            }
            fingerprint
        })
        .collect();
    println!("Finger:{:?}", Instant::now() - start);
    let mut success = 0;
    for (idx, unit1) in units.iter().enumerate() {
        for (jdx, unit2) in units.iter().enumerate() {
            let share_k_mer = kmer_vec[idx]
                .iter()
                .zip(kmer_vec[jdx].iter())
                .filter(|(&x, &y)| x && y)
                .count();
            if share_k_mer < 500 {
                success += 1;
                continue;
            }
            let aln = alignment(unit1, unit2);
            if aln > 1000 {
                score += aln;
            }
        }
    }
    let end = Instant::now();
    println!("{}", success);
    println!("{:?}", end - start);
    println!("{}", score);
}

fn to_index(kmer: &[char]) -> usize {
    kmer.iter().fold(0, |x, y| match y {
        'A' => (x << 2),
        'C' => (x << 2) + 1,
        'G' => (x << 2) + 2,
        'T' => (x << 2) + 3,
        _ => panic!(),
    })
}

fn alignment(seq1: &[char], seq2: &[char]) -> i32 {
    let mut dp = vec![vec![0; seq2.len() + 1]; seq1.len() + 1];
    for i in 1..seq1.len() + 1 {
        for j in 1..seq2.len() + 1 {
            let m = if seq1[i - 1] == seq2[j - 1] { 1 } else { -1 };
            dp[i][j] = (dp[i - 1][j] - 1)
                .max(dp[i][j - 1] - 1)
                .max(dp[i - 1][j - 1] + m);
        }
    }
    dp[seq1.len()][seq2.len()].max(0)
}
