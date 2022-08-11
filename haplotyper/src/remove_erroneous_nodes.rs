pub const MEDIAN_FRAC: usize = 4;
pub const IMPROVE_THR: f64 = 3f64;
pub trait RemoveErroneousNodes {
    fn remove_erroneous_nodes(&mut self);
}

use crate::assemble::copy_number::CoverageCalibrator;
use definitions::*;
use std::collections::{HashMap, HashSet};

fn normalize_node(from: &Node, to: &Node) -> (u64, u64) {
    match from.unit <= to.unit {
        true => (from.unit, to.unit),
        false => (to.unit, from.unit),
    }
}

fn edge_coverage_lengths(ds: &DataSet) -> HashMap<(u64, u64), (usize, i64)> {
    let mut counts: HashMap<_, (usize, i64)> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for edge in read.edges.iter() {
            let key = (edge.from.min(edge.to), edge.from.max(edge.to));
            let entry = counts.entry(key).or_default();
            entry.0 += 1;
            entry.1 += edge.offset;
        }
    }
    counts
}

fn edge_calibed_coverage(ds: &DataSet) -> HashMap<(u64, u64), f64> {
    let edge_coverage_lengths = edge_coverage_lengths(ds);
    let lens: Vec<_> = ds.encoded_reads.iter().map(|r| r.original_length).collect();
    let calib = CoverageCalibrator::new(&lens);
    edge_coverage_lengths
        .into_iter()
        .map(|(key, (obs, totlen))| {
            let len = (totlen / obs as i64).max(0) as usize;
            let calibed = calib.calib(obs, len);
            (key, calibed)
        })
        .collect()
}

fn median_of_coverage<T>(coverages: &HashMap<T, f64>) -> f64 {
    let mut covs: Vec<_> = coverages.values().copied().collect();
    let len = covs.len() / MEDIAN_FRAC;
    let (_, median, _) = covs.select_nth_unstable_by(len, |x, y| x.partial_cmp(y).unwrap());
    *median
}

fn enumerate_edges_to_remove(ds: &DataSet) -> HashMap<(u64, u64), u64> {
    let edge_coverages = edge_calibed_coverage(ds);
    let median = median_of_coverage(&edge_coverages);
    debug!("MEDIAN\t{}", median);
    edge_coverages
        .iter()
        .filter(|(_, &cov)| cov < median / 4f64)
        .filter_map(|(&key, &cov)| {
            // Search reads with (from,fd,to,td) occurence.
            let mut former_neighbor = HashSet::new();
            let mut later_neighbor = HashSet::new();
            for read in ds.encoded_reads.iter() {
                for (i, from) in read.nodes.iter().enumerate().take(read.nodes.len() - 1) {
                    let to = &read.nodes[i + 1];
                    if normalize_node(from, to) == key {
                        if let Some(next) = read.nodes.get(i + 2) {
                            match from.unit <= to.unit {
                                true => former_neighbor.insert(next.unit),
                                false => later_neighbor.insert(next.unit),
                            };
                        }
                        if let Some(prev) = read.nodes.get(i - 1) {
                            match from.unit <= to.unit {
                                true => later_neighbor.insert(prev.unit),
                                false => former_neighbor.insert(prev.unit),
                            };
                        }
                    }
                }
            }
            for next_unit in former_neighbor {
                let probe = if key.0 <= next_unit {
                    (key.0, next_unit)
                } else {
                    (next_unit, key.0)
                };
                let mod_cov = *edge_coverages.get(&probe).unwrap_or(&0f64);
                if IMPROVE_THR * cov < mod_cov {
                    let (from, to) = key;
                    debug!(
                        "REMOVING\t{}\t{}\t{}\t{}\t{}\t{}",
                        from, to, to, next_unit, cov, mod_cov
                    );
                    return Some((key, key.1));
                }
            }
            for prev_unit in later_neighbor {
                let probe = if key.1 <= prev_unit {
                    (key.1, prev_unit)
                } else {
                    (prev_unit, key.1)
                };
                let mod_cov = *edge_coverages.get(&probe).unwrap_or(&0f64);
                if IMPROVE_THR * cov < mod_cov {
                    let (from, to) = key;
                    debug!(
                        "REMOVING\t{}\t{}\t{}\t{}\t{}\t{}",
                        from, to, from, prev_unit, cov, mod_cov
                    );
                    return Some((key, key.0));
                }
            }
            None
        })
        .collect()
}

fn remove_node_idx(read: &EncodedRead, to_remove: &HashMap<(u64, u64), u64>) -> Vec<usize> {
    read.nodes
        .iter()
        .enumerate()
        .take(read.nodes.len() - 1)
        .filter_map(|(i, from)| {
            let to = &read.nodes[i + 1];
            to_remove
                .get(&normalize_node(from, to))
                .map(|&removing| match removing == from.unit {
                    true => i,
                    false => i + 1,
                })
        })
        .collect()
}

impl RemoveErroneousNodes for DataSet {
    fn remove_erroneous_nodes(&mut self) {
        let to_remove = enumerate_edges_to_remove(self);
        for read in self.encoded_reads.iter_mut() {
            let remove_idx = remove_node_idx(read, &to_remove);
            for (offset, i) in remove_idx.into_iter().enumerate() {
                read.remove(i - offset);
            }
        }
    }
}
