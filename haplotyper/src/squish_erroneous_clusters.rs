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
            thr: 0.5,
            count_thr: 5,
        }
    }
}

pub trait SquishErroneousClusters {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig);
}

impl SquishErroneousClusters for DataSet {
    fn squish_erroneous_clusters(&mut self, config: &SquishConfig) {
        let mut unit_pairs: HashMap<_, u32> = HashMap::new();
        for read in self.encoded_reads.iter() {
            for (i, n1) in read.nodes.iter().enumerate() {
                for n2 in read.nodes.iter().skip(i + 1) {
                    let key = (n1.unit.min(n2.unit), n1.unit.max(n2.unit));
                    *unit_pairs.entry(key).or_default() += 1;
                }
            }
        }
        let units: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|n| (n.id, n.cluster_num))
            .collect();
        unit_pairs.retain(|_, val| config.count_thr < *val);
        unit_pairs.retain(|(u1, u2), _| 1 < units[u1] && 1 < units[u2]);
        let cramers_vs: Vec<_> = unit_pairs
            .par_iter()
            .filter_map(|(&(u1, u2), _)| check_correl(self, u1, u2).map(|x| (u1, u2, x)))
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
        let suspic_units: HashSet<_> = unit_to_relvector
            .iter()
            .filter(|(u, _)| !stiff_units.contains(u))
            .filter(|(_, rels)| {
                rels.iter()
                    .any(|(to, rel)| stiff_units.contains(to) && *rel < config.thr)
            })
            .map(|x| x.0)
            .collect();
        debug!("Squish\t{}", suspic_units.len());
        self.selected_chunks
            .iter_mut()
            .filter(|n| suspic_units.contains(&n.id))
            .for_each(|n| n.cluster_num = 1);
        self.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .filter(|n| suspic_units.contains(&n.unit))
            .for_each(|n| {
                n.cluster = 0;
                n.posterior = vec![0f64];
            });
    }
}

fn check_correl(ds: &DataSet, unit1: u64, unit2: u64) -> Option<(f64, usize)> {
    let mut occs = vec![];
    for read in ds.encoded_reads.iter() {
        for node1 in read.nodes.iter().filter(|n| n.unit == unit1) {
            for node2 in read.nodes.iter().filter(|n| n.unit == unit2) {
                occs.push((node1, node2));
            }
        }
    }
    corel(&occs)
}

fn corel(pairs: &[(&Node, &Node)]) -> Option<(f64, usize)> {
    use crate::misc::cramers_v;
    let occs: Vec<_> = pairs
        .iter()
        .map(|(n1, n2)| (n1.cluster as u32, n2.cluster as u32))
        .collect();
    let (first_slot_len, second_slot_len) = occs
        .iter()
        .fold((0, 0), |(fst, snd), &(x, y)| (fst.max(x), snd.max(y)));
    if first_slot_len == 0 || second_slot_len == 0 {
        None
    } else {
        Some((cramers_v(&occs), pairs.len()))
    }
}
