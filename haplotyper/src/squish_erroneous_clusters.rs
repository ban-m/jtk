use definitions::*;
use log::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Copy)]
pub struct SquishConfig {
    thr: f64,
    count_thr: u32,
}

impl SquishConfig {
    pub fn new(thr: f64, count_thr: u32) -> Self {
        Self { thr, count_thr }
    }
}

impl std::default::Default for SquishConfig {
    fn default() -> Self {
        Self {
            thr: 0.8,
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

pub fn classify_units(ds: &DataSet, config: &SquishConfig) -> HashMap<u64, RelClass> {
    let mut unit_pairs: HashMap<_, u32> = HashMap::new();
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
    let cramers_vs: Vec<_> = unit_pairs
        .par_iter()
        .map(|(&(u1, u2), _)| {
            let (cl1, cl2) = (units[&u1], units[&u2]);
            let (rel, count, total) = check_correl(ds, (u1, cl1), (u2, cl2));
            trace!("REL\t{u1}\t{u2}\t{rel:.3}\t{count}\t{total}");
            trace!("REL\t{u2}\t{u1}\t{rel:.3}\t{count}\t{total}");
            (u1, u2, (rel, count))
        })
        .collect();
    let mut unit_to_relvector: HashMap<_, Vec<_>> = HashMap::new();
    for &(u1, u2, (rel, _)) in cramers_vs.iter() {
        unit_to_relvector.entry(u1).or_default().push((u2, rel));
        unit_to_relvector.entry(u2).or_default().push((u1, rel));
    }
    let stiff_units: HashSet<_> = unit_to_relvector
        .iter()
        .filter(|(_, rels)| rels.iter().any(|&(_, rel)| config.thr < rel))
        .map(|x| x.0)
        .collect();
    ds.selected_chunks
        .iter()
        .map(|c| {
            let rels = unit_to_relvector.get(&c.id);
            let touch_stiff = rels.map(|rels| rels.iter().any(|(to, _)| stiff_units.contains(to)));
            if stiff_units.contains(&c.id) {
                (c.id, RelClass::Stiff)
            } else if touch_stiff == Some(true) {
                (c.id, RelClass::Suspicious)
            } else {
                (c.id, RelClass::Isolated)
            }
        })
        .collect()
}

fn check_correl(
    ds: &DataSet,
    (unit1, cl1): (u64, usize),
    (unit2, cl2): (u64, usize),
) -> (f64, usize, usize) {
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
    let mut occs = vec![];
    let mut count = 0;
    for read in ds.encoded_reads.iter() {
        let node1 = read
            .nodes
            .iter()
            .filter(|n| n.unit == unit1)
            .map(|n| n.cluster as u32 + 1)
            .min()
            .unwrap_or(0);
        let node2 = read
            .nodes
            .iter()
            .filter(|n| n.unit == unit2)
            .map(|n| n.cluster as u32 + 1)
            .min()
            .unwrap_or(0);
        if node1 != 0 || node2 != 0 {
            occs.push((node1, node2));
        }
        count += (node1 != 0 && node2 != 0) as usize;
    }
    let rel_value = crate::misc::cramers_v(&occs, (cl1 + 1, cl2 + 1));
    (rel_value, count, occs.len())
}
