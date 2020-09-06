use definitions::*;
use log::*;
use rand::{Rng, SeedableRng};
use rand_distr::ChiSquared;
use rand_distr::Distribution;
use rand_xoshiro::Xoshiro256PlusPlus;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    debug!("Started");
    let dataset: DataSet = serde_json::de::from_reader(ds).unwrap();
    let to_be_used = select_units(&dataset);
    let mut to_be_used: Vec<_> = to_be_used.into_iter().collect();
    to_be_used.sort_by_key(|x| x.0);
    for (unit_id, is_ok) in to_be_used {
        println!("{}\t{}", unit_id, is_ok);
    }
    Ok(())
}

use std::collections::HashMap;
fn select_units(dataset: &DataSet) -> HashMap<u64, bool> {
    dataset
        .selected_chunks
        .iter()
        .map(|unit| (unit.id, is_ok_unit(unit, dataset)))
        .collect()
}

fn is_ok_unit(unit: &Unit, dataset: &DataSet) -> bool {
    debug!("UNIT\t{}", unit.id);
    let mut connections: HashMap<_, Vec<_>> = HashMap::new();
    for read in dataset.encoded_reads.iter() {
        let read = read.nodes.iter().enumerate();
        for (idx, node) in read.clone().filter(|(_, n)| n.unit == unit.id) {
            let cluster = node.cluster;
            for (_, n) in read.clone().filter(|&(i, _)| i != idx) {
                connections
                    .entry(n.unit)
                    .or_default()
                    .push((cluster, n.cluster));
            }
        }
    }
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(43289);
    if connections.is_empty() {
        true
    } else {
        connections.iter().any(|(unit, pairs)| {
            let s: u64 = rng.gen();
            debug!("PAIR\t{}", unit);
            is_significant(&pairs, s, 2, 2)
        })
    }
}

fn is_significant(clusters: &[(u64, u64)], seed: u64, k1: u64, k2: u64) -> bool {
    let mut dump: HashMap<_, u32> = HashMap::new();
    for &cl in clusters.iter() {
        *dump.entry(cl).or_default() += 1;
    }
    debug!("DUMP\t{:?}", dump);
    let degree_of_freedom = (k1 - 1) * (k2 - 1);
    assert!(degree_of_freedom > 0);
    let k1_fraction: Vec<_> = (0..k1)
        .map(|c1| clusters.iter().filter(|&&(x, _)| x == c1).count())
        .collect();
    let k2_fraction: Vec<_> = (0..k2)
        .map(|c2| clusters.iter().filter(|&&(_, y)| y == c2).count())
        .collect();
    if k1_fraction.iter().any(|&x| x == 0) || k2_fraction.iter().any(|&x| x == 0) {
        return false;
    }
    let len = clusters.len() as f64;
    let chi_squared = (0..k1 as usize)
        .map(|c1| {
            (0..k2 as usize)
                .map(|c2| {
                    let c_12 = clusters
                        .iter()
                        .filter(|&&(x, y)| x == c1 as u64 && y == c2 as u64)
                        .count();
                    let expected = k1_fraction[c1] * k2_fraction[c2];
                    let expected = expected as f64 / len;
                    (c_12 as f64 - expected).powi(2i32) / expected
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    let sample_num = 10_000;
    let distr = ChiSquared::new(degree_of_freedom as f64).unwrap();
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let hit = distr
        .sample_iter(&mut rng)
        .take(sample_num)
        .filter(|&x| x > chi_squared)
        .count();
    debug!("CHI\t{}\t{:.2}", hit, chi_squared);
    hit < sample_num * 1 / 100
}
