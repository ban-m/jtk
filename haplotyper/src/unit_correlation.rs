use definitions::*;
use rand::SeedableRng;
use rand_distr::ChiSquared;
use rand_distr::Distribution;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;
use std::collections::HashMap;
pub fn calc_p_values(ds: &DataSet, thr: u32) -> HashMap<u64, f64> {
    let mut pair_units: HashMap<_, HashMap<_, u32>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for (i, n) in read.nodes.iter().enumerate() {
            for m in read.nodes.iter().skip(i + 1) {
                let (n, m) = (n.unit.min(m.unit), n.unit.max(m.unit));
                *(pair_units.entry(n).or_default().entry(m).or_default()) += 1;
            }
        }
    }
    let pvalues: Vec<_> = pair_units
        .into_par_iter()
        .map(|(id, connections)| {
            let mut pvalues = vec![];
            for (&jd, _) in connections.iter().filter(|&(_, &count)| thr < count) {
                let mut pairs = vec![];
                for read in ds.encoded_reads.iter() {
                    let read = read.nodes.iter().enumerate();
                    for (idx, node) in read.clone().filter(|(_, n)| n.unit == id) {
                        let cluster = node.cluster;
                        for (_, n) in read.clone().filter(|&(i, n)| i != idx && n.unit == jd) {
                            pairs.push((cluster, n.cluster));
                        }
                    }
                }
                let p_value = calc_p_value(&pairs, id * jd);
                pvalues.push((id, jd, p_value));
            }
            pvalues
        })
        .collect();
    let mut min_p_value: HashMap<u64, f64> = HashMap::new();
    for pvalue in pvalues {
        for (id, jd, pvalue) in pvalue {
            if let Some(p) = min_p_value.get_mut(&id) {
                *p = p.min(pvalue);
            } else {
                min_p_value.insert(id, pvalue);
            }
            if let Some(p) = min_p_value.get_mut(&jd) {
                *p = p.min(pvalue);
            } else {
                min_p_value.insert(jd, pvalue);
            }
        }
    }
    min_p_value
}

// Return true if the unit is informative, or containing variants and
// propery clustered, assumed.
pub fn select_informative_units(dataset: &DataSet, p_value: f64) -> HashMap<u64, bool> {
    dataset
        .selected_chunks
        .iter()
        .map(|unit| (unit.id, is_ok_unit(unit, dataset, p_value)))
        .collect()
}

fn is_ok_unit(unit: &definitions::Unit, dataset: &DataSet, p_value: f64) -> bool {
    trace!("UNIT\t{}", unit.id);
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
    if connections.is_empty() {
        // It is OK to have isolated node.
        true
    } else if connections.iter().all(|(_, pairs)| pairs.len() < 5) {
        true
    } else {
        connections
            .iter()
            .filter(|(_, pairs)| 5 <= pairs.len())
            .any(|(&p, pairs)| {
                let pair_unit = dataset.selected_chunks.iter().find(|u| u.id == p).unwrap();
                let seed = unit.id * pair_unit.id;
                calc_p_value(&pairs, seed) < p_value
            })
    }
}

pub fn calc_p_value(clusters: &[(u64, u64)], seed: u64) -> f64 {
    let mut cluster1: HashMap<_, u32> = HashMap::new();
    let mut cluster2: HashMap<_, u32> = HashMap::new();
    let mut dump: HashMap<_, u32> = HashMap::new();
    for &cl in clusters.iter() {
        *dump.entry(cl).or_default() += 1;
        *cluster1.entry(cl.0).or_default() += 1;
        *cluster2.entry(cl.1).or_default() += 1;
    }
    let (k1, k2) = (cluster1.len(), cluster2.len());
    trace!("DUMP\t{:?}", dump);
    let degree_of_freedom = (k1 - 1) * (k2 - 1);
    if degree_of_freedom == 0 {
        return 1f64;
    }
    assert!(degree_of_freedom > 0);
    let len = clusters.len() as f64;
    let chi_squared = cluster1
        .iter()
        .map(|(&key1, c1)| {
            cluster2
                .iter()
                .map(|(&key2, c2)| {
                    let c_12 = *dump.get(&(key1, key2)).unwrap_or(&0);
                    let expected = (c1 * c2) as f64 / len;
                    (c_12 as f64 - expected).powi(2i32) / expected
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    // let chi_squared = (0..k1 as usize)
    //     .map(|c1| {
    //         (0..k2 as usize)
    //             .map(|c2| {
    //                 let c_12 = clusters
    //                     .iter()
    //                     .filter(|&&(x, y)| x == c1 as u64 && y == c2 as u64)
    //                     .count();
    //                 let expected = k1_fraction[c1] * k2_fraction[c2];
    //                 let expected = expected as f64 / len;
    //                 (c_12 as f64 - expected).powi(2i32) / expected
    //             })
    //             .sum::<f64>()
    //     })
    //     .sum::<f64>();
    let sample_num = 100_000;
    let distr = ChiSquared::new(degree_of_freedom as f64).unwrap();
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let hit = distr
        .sample_iter(&mut rng)
        .take(sample_num)
        .filter(|&x| x > chi_squared)
        .count();
    trace!("CHI\t{}\t{:.2}", hit, chi_squared);
    hit as f64 / sample_num as f64
}
