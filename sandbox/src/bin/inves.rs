#![allow(unused_imports)]
use definitions::*;
use haplotyper::assemble::ditch_graph::Focus;
use haplotyper::assemble::*;
use haplotyper::determine_units::DetermineUnit;
use nalgebra::RealField;
use rand::SeedableRng;
use rand_distr::CauchyError;
use rand_xoshiro::Xoshiro256StarStar;
use sandbox::IS_MOCK;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    let units: Vec<u64> = args[2..].iter().map(|x| x.parse().unwrap()).collect();
    let num_cluster: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    println!("{units:?}");
    for (i, &unit1) in units.iter().enumerate() {
        for &unit2 in units.iter().skip(i + 1) {
            let mut occs = vec![vec![0; num_cluster[&unit2]]; num_cluster[&unit1]];
            for read in ds.encoded_reads.iter() {
                for node1 in read.nodes.iter().filter(|n| n.unit == unit1) {
                    for node2 in read.nodes.iter().filter(|n| n.unit == unit2) {
                        occs[node1.cluster as usize][node2.cluster as usize] += 1;
                    }
                }
            }
            if occs.iter().flatten().sum::<u32>() == 0 {
                continue;
            }
            println!("{unit1}\t{unit2}");
            for row in occs.iter() {
                let row: Vec<_> = row.iter().map(|x| format!("{x}")).collect();
                println!("\t{}", row.join("\t"));
            }
        }
    }
    Ok(())
}

// pub fn estimate_minimum_gain(
//     hmm: &kiley::hmm::guided::PairHiddenMarkovModel,
//     error_rate: ErrorRate,
// ) {
//     const SEED: u64 = 23908;
//     const SAMPLE_NUM: usize = 1000;
//     const SEQ_NUM: usize = 500;
//     const LEN: usize = 100;
//     const BAND: usize = 25;
//     const FRAC: f64 = 0.25;
//     const PICK: usize = (SAMPLE_NUM as f64 * FRAC) as usize;
//     let prof = kiley::gen_seq::Profile {
//         sub: error_rate.mismatch,
//         del: error_rate.del,
//         ins: error_rate.ins,
//     };
//     let mut medians: Vec<_> = (0..SAMPLE_NUM)
//         .map(|seed| {
//             let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(SEED + seed as u64);
//             let hap1 = kiley::gen_seq::generate_seq(&mut rng, LEN);
//             let hap2 = match seed % 2 == 0 {
//                 true => kiley::gen_seq::introduce_errors(&hap1, &mut rng, 0, 0, 1),
//                 false => kiley::gen_seq::introduce_errors(&hap1, &mut rng, 0, 1, 0),
//             };
//             let mut lks: Vec<_> = (0..SEQ_NUM)
//                 .map(|_| {
//                     let read = kiley::gen_seq::introduce_randomness(&hap1, &mut rng, &prof);
//                     let lk_base = hmm.likelihood(&hap1, &read, BAND);
//                     let lk_diff = hmm.likelihood(&hap2, &read, BAND);
//                     lk_base - lk_diff
//                 })
//                 .collect();
//             for (i, x) in lks.iter().enumerate() {
//                 println!("{seed}\t{i}\t{x}");
//             }
//             *lks.select_nth_unstable_by(SEQ_NUM / 2, |x, y| x.partial_cmp(y).unwrap())
//                 .1
//         })
//         .collect();
//     medians.sort_by(|x, y| x.partial_cmp(y).unwrap());
//     eprintln!("{}", medians[PICK]);
// }
