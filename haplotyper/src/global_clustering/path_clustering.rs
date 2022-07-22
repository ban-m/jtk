#![allow(dead_code)]
use rand::Rng;
pub struct Graph {
    // TODO: Maybe it should be the sum of the read length?
    haploid_reads: Vec<u32>,
    // Ploidy. For diploid, 2.
    ploidy: usize,
    // Choises.
    haplotypes: Vec<usize>,
    cluster_num: Vec<usize>,
    //Prior.
    mock_count: u32,
    // Haplotype -> Position -> Cluster
    counts: Vec<Vec<Vec<u32>>>,
    // Total counts: Haplotype -> Position
    total_counts: Vec<Vec<u32>>,
    reads: Vec<Vec<(usize, usize)>>,
    assignments: Vec<usize>,
}

pub trait ToPath {
    // Into (Position, Cluster)
    fn to_path(&self) -> Vec<(usize, usize)>;
}

impl ToPath for definitions::EncodedRead {
    fn to_path(&self) -> Vec<(usize, usize)> {
        self.nodes
            .iter()
            .map(|n| (n.unit as usize, n.cluster as usize))
            .collect()
    }
}

impl Graph {
    pub fn new<R: Rng, T: ToPath>(
        reads: &[T],
        // If some position(unit) is missing, it should be zero.
        cluster_num: &[usize],
        ploidy: usize,
        mock_count: u32,
        r: &mut R,
    ) -> Self {
        let reads: Vec<Vec<_>> = reads
            .iter()
            .map(ToPath::to_path)
            .map(|r| r.into_iter().map(|(u, c)| (u as usize, c)).collect())
            .collect();
        let assignments: Vec<usize> = reads.iter().map(|_| r.gen::<usize>() % ploidy).collect();
        let mut counts: Vec<Vec<Vec<u32>>> = (0..ploidy)
            .map(|_| cluster_num.iter().map(|&c| vec![mock_count; c]).collect())
            .collect();
        let mut haploid_reads = vec![mock_count; ploidy];
        let mut total_counts: Vec<Vec<u32>> = (0..ploidy)
            .map(|_| cluster_num.iter().map(|&c| mock_count * c as u32).collect())
            .collect();
        for (read, &a) in reads.iter().zip(assignments.iter()) {
            haploid_reads[a] += 1;
            for &(unit, cluster) in read.iter() {
                counts[a][unit][cluster] += 1;
                total_counts[a][unit] += 1;
            }
        }
        let haplotypes: Vec<_> = (0..ploidy).collect();
        Self {
            haplotypes,
            haploid_reads,
            ploidy,
            cluster_num: cluster_num.to_vec(),
            mock_count,
            counts,
            total_counts,
            reads,
            assignments,
        }
    }
    pub fn sanity_check(&self) {
        for (total, counts) in self.total_counts.iter().zip(self.counts.iter()) {
            for (&t, c) in total.iter().zip(counts) {
                assert_eq!(t, c.iter().sum::<u32>());
            }
        }
    }
    pub fn gibbs_sampling<R: Rng>(&mut self, rng: &mut R, iter_num: usize) -> (f64, &[usize]) {
        let mut idx = 0;
        for _ in 0..iter_num {
            self.sanity_check();
            self.update_ith(idx, rng);
            idx += 1;
            idx %= self.reads.len();
        }
        (self.lk(), self.assignments())
    }
    fn update_ith<R: Rng>(&mut self, idx: usize, rng: &mut R) {
        // Remove the i-th read.
        let ans = self.assignments[idx];
        self.haploid_reads[ans] -= 1;
        for &(unit, cluster) in self.reads[idx].iter() {
            self.counts[ans][unit][cluster] -= 1;
            self.total_counts[ans][unit] -= 1;
        }
        // Compute the i-th read's likelihood.
        let lks: Vec<_> = self.lk_of_the_ith(idx);
        // Select cluster.
        let lk = logsumexp(&lks);
        let weights: Vec<_> = lks.iter().map(|x| (x - lk).exp()).collect();
        assert!((1. - weights.iter().sum::<f64>()).abs() < 0.001);
        use rand::seq::SliceRandom;
        let k = *self
            .haplotypes
            .choose_weighted(rng, |&k| weights[k])
            .unwrap();
        // Add the i-th read.
        for &(unit, cluster) in self.reads[idx].iter() {
            self.counts[k][unit][cluster] += 1;
            self.total_counts[k][unit] += 1;
        }
        self.haploid_reads[k] += 1;
        self.assignments[idx] = k;
    }
    fn lk_of_the_ith(&self, i: usize) -> Vec<f64> {
        let total_hap = self.haploid_reads.iter().sum::<u32>();
        let read = &self.reads[i];
        self.haploid_reads
            .iter()
            .zip(self.counts.iter())
            .zip(self.total_counts.iter())
            .map(|((&f, counts), tot)| {
                let hap_frac = (f as f64 / total_hap as f64).ln();
                let read_lk = read
                    .iter()
                    .map(|&(u, c)| (counts[u][c] as f64 / tot[u] as f64).ln())
                    .sum::<f64>();
                hap_frac + read_lk
            })
            .collect()
    }
    pub fn lk(&self) -> f64 {
        (0..self.reads.len())
            .map(|i| crate::misc::logsumexp(&self.lk_of_the_ith(i)))
            .sum::<f64>()
    }
    pub fn assignments(&self) -> &[usize] {
        &self.assignments
    }
}

impl std::fmt::Display for Graph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let fractions: Vec<_> = self
            .haploid_reads
            .iter()
            .map(|&c| c as f64 / self.reads.len() as f64)
            .enumerate()
            .map(|(i, x)| format!("{}:{:.3}", i, x))
            .collect();
        writeln!(f, "{}", fractions.join(","))?;
        for (k, cs) in self.counts.iter().enumerate() {
            writeln!(f, "{}-th haplotype", k)?;
            for (p, c) in cs.iter().enumerate() {
                writeln!(f, "{}\t{:?}", p, c)?;
            }
        }
        write!(f, "MockCount:{}, LK:{}", self.mock_count, self.lk())
    }
}

// struct GNode {
//     edges: Vec<GEdge>,
//     total: f64,
// }

// impl GNode {
//     fn new() -> Self {
//         Self {
//             edges: vec![],
//             total: 0.,
//         }
//     }
//     fn register(&mut self, to: usize) {
//         self.total += 1.;
//         if let Some(edge) = self.edges.iter_mut().find(|e| e.to == to) {
//             edge.weight += 1.;
//         } else {
//             self.edges.push(GEdge { to, weight: 1. });
//         }
//     }
// }

// struct GEdge {
//     to: usize,
//     weight: f64,
// }

// pub fn path_clustering(reads: &[Vec<usize>], init: &[usize]) -> Vec<usize> {
//     let cluster_num = *init.iter().max().unwrap() + 1;
//     let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(100);
//     let mut cluster = init.to_vec();
//     assert!(!reads.is_empty());
//     debug!("Start path clustering");
//     let mut count = 0;
//     //let stable_thr = (reads.len() / 100).max(4) as u32;
//     let stable_thr = 2;
//     let len = reads.len();
//     while count < 3 {
//         let change_num = rand::seq::index::sample(&mut rng, len, len)
//             .into_iter()
//             .map(|idx| {
//                 let models: Vec<_> = (0..cluster_num)
//                     .map(|cl| Graph::new(&reads, &cluster, cl, &mut rng, idx))
//                     .collect();
//                 let scores: Vec<_> = models.iter().map(|m| m.score(&reads[idx])).collect();
//                 let (argmax, _) = scores
//                     .iter()
//                     .enumerate()
//                     .max_by(|x, y| match (x.1).partial_cmp(&y.1) {
//                         Some(res) => res,
//                         None => panic!("{}\t{:?}", line!(), scores),
//                     })
//                     .unwrap();
//                 let changed = 1 - (cluster[idx] == argmax) as u32;
//                 cluster[idx] = argmax;
//                 changed
//             })
//             .sum::<u32>();
//         count += (change_num < stable_thr) as u32;
//         count *= (change_num < stable_thr) as u32;
//         eprintln!("ChangeNum:{},{}", change_num, stable_thr);
//     }
//     cluster
// }
// #[cfg(test)]
// mod tests {
//     #[derive(Clone, Copy, Debug)]
//     struct TestConfig {
//         cl: usize,
//         num: usize,
//         fail: f64,
//         skip: f64,
//         max_len: usize,
//         min_len: usize,
//         unit_len: usize,
//     }
//     use super::*;
//     use std::collections::{HashMap, HashSet};
//     fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> (Vec<Vec<usize>>, Vec<usize>) {
//         let TestConfig {
//             cl,
//             num,
//             fail,
//             skip,
//             max_len,
//             min_len,
//             unit_len,
//         } = conf;
//         let mut answer = vec![];
//         let mut reads = vec![];
//         for i in 0..num {
//             let cluster = (i % cl) as u64;
//             let len = r.gen::<usize>() % (max_len - min_len) + min_len;
//             let start = r.gen::<usize>() % (unit_len - len);
//             let units: Vec<_> = if r.gen_bool(0.5) {
//                 (start..=start + len).collect()
//             } else {
//                 let start = start + len;
//                 (start - len..=start).rev().collect()
//             };
//             let mut read = vec![];
//             for unit in units {
//                 if r.gen_bool(skip) {
//                     continue;
//                 } else if r.gen_bool(fail) {
//                     let cluster = r.gen::<u64>() % (cl + 1) as u64;
//                     read.push((unit as u64, cluster));
//                 } else {
//                     read.push((unit as u64, cluster));
//                 }
//             }
//             answer.push(cluster as usize);
//             reads.push(read);
//         }
//         let units: HashSet<_> = reads.iter().flat_map(|r| r.iter().copied()).collect();
//         let units: HashMap<(u64, u64), usize> =
//             units.into_iter().enumerate().map(|(x, y)| (y, x)).collect();
//         let reads: Vec<Vec<_>> = reads
//             .into_iter()
//             .map(|read| read.into_iter().map(|x| units[&x]).collect())
//             .collect();
//         (reads, answer)
//     }
//     fn path_clustering_test() {
//         let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
//         let conf = TestConfig {
//             cl: 2,
//             num: 200,
//             fail: 0.0,
//             skip: 0.0,
//             max_len: 20,
//             min_len: 10,
//             unit_len: 50,
//         };
//         let (reads, answer) = gen_dataset(&mut rng, conf);
//         let init: Vec<_> = answer
//             .iter()
//             .map(|&ans| {
//                 if ans == 0 {
//                     ans
//                 } else if rng.gen_bool(0.5) {
//                     1
//                 } else {
//                     2
//                 }
//             })
//             .collect();
//         let preds = path_clustering(&reads, &init);
//         assert_eq!(preds.len(), answer.len());
//         let correct = answer
//             .iter()
//             .zip(preds.iter())
//             .filter(|&(ans, pred)| ans == pred)
//             .count();
//         eprintln!("{}/{}", correct, reads.len());
//         for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
//             eprintln!("{}\t{}\t{}", ans, assign, read.len());
//         }
//         let correct = correct.max(reads.len() - correct);
//         assert!(correct > reads.len() * 8 / 10);
//     }
// }
