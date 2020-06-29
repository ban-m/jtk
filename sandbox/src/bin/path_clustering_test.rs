use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

#[derive(Clone, Copy, Debug)]
struct TestConfig {
    cl: usize,
    num: usize,
    fail: f64,
    skip: f64,
    max_len: usize,
    min_len: usize,
    unit_len: usize,
}
use std::collections::{HashMap, HashSet};
fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> (Vec<Vec<usize>>, Vec<usize>) {
    let TestConfig {
        cl,
        num,
        fail,
        skip,
        max_len,
        min_len,
        unit_len,
    } = conf;
    let mut answer = vec![];
    let mut reads = vec![];
    for i in 0..num {
        let cluster = (i % cl) as u64;
        let len = r.gen::<usize>() % (max_len - min_len) + min_len;
        let start = r.gen::<usize>() % (unit_len - len);
        let units: Vec<_> = if r.gen_bool(0.5) {
            (start..=start + len).collect()
        } else {
            let start = start + len;
            (start - len..=start).rev().collect()
        };
        let mut read = vec![];
        for unit in units {
            if r.gen_bool(skip) {
                continue;
            } else if r.gen_bool(fail) {
                let cluster = r.gen::<u64>() % (cl + 1) as u64;
                read.push((unit as u64, cluster));
            } else {
                read.push((unit as u64, cluster));
            }
        }
        answer.push(cluster as usize);
        reads.push(read);
    }
    let units: HashSet<_> = reads.iter().flat_map(|r| r.iter().copied()).collect();
    let units: HashMap<(u64, u64), usize> =
        units.into_iter().enumerate().map(|(x, y)| (y, x)).collect();
    let reads: Vec<Vec<_>> = reads
        .into_iter()
        .map(|read| read.into_iter().map(|x| units[&x]).collect())
        .collect();
    (reads, answer)
}
fn main() {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
    let conf = TestConfig {
        cl: 2,
        num: 200,
        fail: 0.0,
        skip: 0.0,
        max_len: 20,
        min_len: 10,
        unit_len: 50,
    };
    let (reads, answer) = gen_dataset(&mut rng, conf);
    let init: Vec<_> = answer
        .iter()
        .map(|&ans| {
            if ans == 0 {
                ans
            } else if rng.gen_bool(0.5) {
                1
            } else {
                2
            }
        })
        .collect();
    let preds = haplotyper::path_clustering(&reads, &init);
    assert_eq!(preds.len(), answer.len());
    let correct = answer
        .iter()
        .zip(preds.iter())
        .filter(|&(ans, pred)| ans == pred)
        .count();
    eprintln!("{}/{}", correct, reads.len());
    for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
        eprintln!("{}\t{}\t{}", ans, assign, read.len());
    }
    let correct = correct.max(reads.len() - correct);
    assert!(correct > reads.len() * 8 / 10);
}
