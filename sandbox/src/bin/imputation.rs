#![allow(dead_code)]
use definitions::*;
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
struct DataConfig {
    trim: bool,
    cluster_error: f64,
    unit_number: usize,
    min_cluster_num: usize,
    max_cluster_num: usize,
    read_number: usize,
}

impl DataConfig {
    fn create<R: Rng>(&self, rng: &mut R) -> (Vec<Unit>, Vec<EncodedRead>, Vec<usize>, usize) {
        let units: Vec<_> = rand::seq::index::sample(rng, self.unit_number * 2, self.unit_number)
            .into_iter()
            .map(|_| {
                let id: u64 = rng.gen_range(0, self.unit_number * 2) as u64;
                let seq = String::new();
                let cluster_num = rng.gen_range(self.min_cluster_num, self.max_cluster_num + 1);
                Unit {
                    id,
                    seq,
                    cluster_num,
                }
            })
            .collect();
        let (reads, answer): (Vec<_>, Vec<_>) = (0..self.read_number as u64)
            .map(|id| self.gen_reads(id, rng, &units))
            .unzip();
        let focal_unit = units.len() / 2;
        let num_cluster = units[focal_unit].cluster_num;
        (units, reads, answer, num_cluster)
    }
    fn gen_reads<R: Rng>(&self, id: u64, rng: &mut R, units: &[Unit]) -> (EncodedRead, usize) {
        let focal_unit = units.len() / 2;
        let num_cluster = units[focal_unit].cluster_num;
        let answer = rng.gen_range(0, num_cluster);
        let (start, end) = if self.trim {
            let start = rng.gen_range(0, focal_unit - 1);
            let end = rng.gen_range(focal_unit + 1, units.len() + 1);
            (start, end)
        } else {
            (0, units.len())
        };
        let mut nodes: Vec<_> = units[start..end]
            .iter()
            .map(|ref_unit| {
                let cluster = if rng.gen_bool(self.cluster_error) {
                    rng.gen_range(0, ref_unit.cluster_num) as u64
                } else {
                    (answer % ref_unit.cluster_num) as u64
                };
                Node {
                    position_from_start: 0,
                    unit: ref_unit.id,
                    cluster,
                    seq: String::new(),
                    is_forward: true,
                    cigar: vec![],
                }
            })
            .collect();
        if rng.gen_bool(0.5) {
            nodes.reverse();
        }
        let read = EncodedRead {
            original_length: 0,
            leading_gap: 0,
            trailing_gap: 0,
            id,
            edges: vec![],
            nodes,
        };
        (read, answer)
    }
}

fn main() {
    for seed in 0..100 {
        let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
        let config = DataConfig {
            trim: true,
            cluster_error: 0.1,
            unit_number: 20,
            min_cluster_num: 2,
            max_cluster_num: 2,
            read_number: 20,
        };
        let (units, reads, answer, cluster_num) = config.create(&mut rng);
        let reads: Vec<_> = reads.iter().collect();
        use haplotyper::em_correction::clustering;
        let focal = units[units.len() / 2].id;
        let config = haplotyper::em_correction::Config::new(10, 2190, cluster_num, focal, 5);
        let preds = clustering(&units, &reads, &config);
        use std::collections::HashMap;
        let mut count: HashMap<_, u32> = HashMap::new();
        for ((_read, ans), pred) in reads.iter().zip(answer).zip(preds) {
            *count.entry((ans, pred)).or_default() += 1;
            // let line: Vec<_> = read
            //     .nodes
            //     .iter()
            //     .map(|n| format!("{:2}-{:2}", n.unit, n.cluster))
            //     .collect();
            // println!("READ\t{}\t{}\t{}\t{}", read.id, ans, pred, line.join(" "));
        }
        print!("{}\t", seed);
        for pred in 0..cluster_num {
            for ans in 0..cluster_num {
                print!("{}\t", count.get(&(ans, pred)).unwrap_or(&0));
            }
        }
        println!();
    }
}
