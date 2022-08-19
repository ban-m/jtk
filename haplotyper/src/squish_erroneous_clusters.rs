use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct SquishConfig {
    frac: f64,
    count_thr: usize,
}

impl SquishConfig {
    pub fn new(frac: f64, count_thr: usize) -> Self {
        Self { frac, count_thr }
    }
}

impl std::default::Default for SquishConfig {
    fn default() -> Self {
        Self {
            frac: 0.01,
            count_thr: 10,
        }
    }
}

pub trait SquishErroneousClusters {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig);
}

impl SquishErroneousClusters for DataSet {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig) {
        let unit_class = classify_units(self, config);
        self.selected_chunks
            .iter_mut()
            .filter(|n| matches!(unit_class.get(&n.id), Some(RelClass::Suspicious)))
            .for_each(|n| n.cluster_num = 1);
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .filter(|n| matches!(unit_class.get(&n.unit), Some(RelClass::Suspicious)))
            .for_each(|n| {
                n.cluster = 0;
                n.posterior = vec![0f64];
            });
    }
}

#[derive(Debug, Clone, Copy)]
pub enum RelClass {
    Stiff,
    Isolated,
    Suspicious,
}

impl std::fmt::Display for RelClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RelClass::Stiff => write!(f, "Stiff"),
            RelClass::Isolated => write!(f, "Isolated"),
            RelClass::Suspicious => write!(f, "Suspicious"),
        }
    }
}

fn classify_units(ds: &DataSet, config: &SquishConfig) -> HashMap<u64, RelClass> {
    let mut unit_pairs: HashMap<_, usize> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for (i, n1) in read.nodes.iter().enumerate() {
            for n2 in read.nodes.iter().skip(i + 1) {
                let key = (n1.unit.min(n2.unit), n1.unit.max(n2.unit));
                *unit_pairs.entry(key).or_default() += 1;
            }
        }
    }
    let units: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|n| (n.id, n.cluster_num))
        .collect();
    unit_pairs.retain(|_, val| config.count_thr < *val);
    unit_pairs.retain(|(u1, u2), _| 1 < units[u1] && 1 < units[u2]);
    trace!("REL\tUnit1\tUnit2\tRel\tCount\tTotal");
    let adj_rand_indices: Vec<_> = unit_pairs
        .par_iter()
        .map(|(&(u1, u2), _)| {
            let (cl1, cl2) = (units[&u1], units[&u2]);
            let (rel, count, _total) = check_correl(ds, (u1, cl1), (u2, cl2));
            // trace!("REL\t{u1}\t{u2}\t{rel:.3}\t{count}\t{total}");
            // trace!("REL\t{u2}\t{u1}\t{rel:.3}\t{count}\t{total}");
            (u1, u2, (rel, count))
        })
        .collect();
    let mut unit_to_relvector: HashMap<_, Vec<_>> = HashMap::new();
    for &(u1, u2, (rel, _)) in adj_rand_indices.iter() {
        unit_to_relvector.entry(u1).or_default().push((u2, rel));
        unit_to_relvector.entry(u2).or_default().push((u1, rel));
    }
    let thr = get_thr(&unit_to_relvector, config);
    debug!("SUPRESS\tTHR\t{thr:.3}");
    let stiff_units: HashSet<_> = unit_to_relvector
        .iter()
        .filter(|(_, rels)| rels.iter().any(|&(_, rel)| thr < rel))
        .map(|x| x.0)
        .collect();
    ds.selected_chunks
        .iter()
        .map(|c| {
            let rels = unit_to_relvector.get(&c.id);
            // let max = rels
            //     .and_then(|xs| {
            //         xs.iter()
            //             .map(|x| x.1)
            //             .max_by(|x, y| x.partial_cmp(y).unwrap())
            //     })
            //     .unwrap_or(0f64);
            // let len = rels.map(|x| x.len()).unwrap_or(0);
            // debug!("SQI\t{}\t{}\t{max}\t{len}", c.id, c.cluster_num);
            let touch_stiff = rels.map(|rels| rels.iter().any(|(to, _)| stiff_units.contains(to)));
            if stiff_units.contains(&c.id) || 2 < c.copy_num {
                (c.id, RelClass::Stiff)
            } else if touch_stiff == Some(true) {
                let max = rels
                    .unwrap()
                    .iter()
                    .map(|x| x.1)
                    .max_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap();
                let len = rels.unwrap().len();
                trace!("SUSPICOUS\t{}\t{}\t{:.3}\t{}", c.id, c.copy_num, max, len);
                (c.id, RelClass::Suspicious)
            } else {
                (c.id, RelClass::Isolated)
            }
        })
        .collect()
}

fn get_thr(unit_to_relvector: &HashMap<u64, Vec<(u64, f64)>>, config: &SquishConfig) -> f64 {
    let mut max_indices: Vec<f64> = unit_to_relvector
        .values()
        .filter_map(|values| {
            values
                .iter()
                .map(|x| x.1)
                .max_by(|x, y| x.partial_cmp(y).unwrap())
        })
        .collect();
    max_indices.sort_by(|x, y| x.partial_cmp(y).unwrap());
    debug!("{:?},{}", &max_indices[..10], max_indices.len());
    max_indices[(config.frac * max_indices.len() as f64).ceil() as usize]
}

fn check_correl(
    ds: &DataSet,
    (unit1, cl1): (u64, usize),
    (unit2, cl2): (u64, usize),
) -> (f64, usize, usize) {
    let (mut c1, mut c2) = (vec![], vec![]);
    for read in ds.encoded_reads.iter() {
        let node1 = read
            .nodes
            .iter()
            .filter(|n| n.unit == unit1)
            .map(|n| n.cluster as usize)
            .min();
        let node2 = read
            .nodes
            .iter()
            .filter(|n| n.unit == unit2)
            .map(|n| n.cluster as usize)
            .min();
        if let (Some(n1), Some(n2)) = (node1, node2) {
            c1.push(n1);
            c2.push(n2);
        }
    }

    if c1.is_empty() {
        return (0f64, c1.len(), c2.len());
    }
    let c1_is_same = c1.iter().all(|&x| x == c1[0]);
    let c2_is_same = c2.iter().all(|&x| x == c2[0]);
    let rel_value = match (c1_is_same && c2_is_same, cl1 == 1 && cl2 == 1) {
        (true, true) => 0f64,
        (true, false) => 1f64,
        (false, _) => crate::misc::adjusted_rand_index(&c1, &c2),
    };
    if rel_value.is_nan() {
        panic!("\n{:?}\n{:?}", c1, c2);
    }
    (rel_value, c1.len(), c1.len())
    // let mut occs = vec![];
    // let mut count = 0;
    // for read in ds.encoded_reads.iter() {
    //     let node1 = read
    //         .nodes
    //         .iter()
    //         .filter(|n| n.unit == unit1)
    //         .map(|n| n.cluster as u32)
    //         .min();
    //     let node2 = read
    //         .nodes
    //         .iter()
    //         .filter(|n| n.unit == unit2)
    //         .map(|n| n.cluster as u32)
    //         .min();
    //     if let (Some(n1), Some(n2)) = (node1, node2) {
    //         occs.push((n1, n2));
    //         count += 1;
    //     }
    // }
    // let rel_value = crate::misc::cramers_v(&occs, (cl1, cl2));
    // (rel_value, count, occs.len())
}
