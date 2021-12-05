pub const MEDIAN_FRAC: usize = 4;
pub const IMPROVE_THR: f64 = 3f64;
pub trait RemoveErroneousNodes {
    fn remove_erroneous_nodes(&mut self);
}

use crate::assemble::copy_number::CoverageCalibrator;
use definitions::*;
use std::collections::{HashMap, HashSet};
impl RemoveErroneousNodes for DataSet {
    #[allow(clippy::needless_collect)]
    fn remove_erroneous_nodes(&mut self) {
        let mut counts: HashMap<_, (usize, i64)> = HashMap::new();
        let normalize = |from: &Node, to: &Node| -> (u64, u64) {
            match from.unit <= to.unit {
                true => (from.unit, to.unit),
                false => (to.unit, from.unit),
            }
        };
        for read in self.encoded_reads.iter() {
            for (i, from) in read.nodes.iter().enumerate().take(read.nodes.len() - 1) {
                let to = &read.nodes[i + 1];
                let entry = counts.entry(normalize(from, to)).or_default();
                entry.0 += 1;
                entry.1 += to.position_from_start as i64
                    - (from.position_from_start - from.query_length()) as i64;
            }
        }
        // Convert to the "calibrated occurance".
        let lens: Vec<_> = self
            .encoded_reads
            .iter()
            .map(|r| r.original_length)
            .collect();
        let calib = CoverageCalibrator::new(&lens);
        let calib_coverage: HashMap<_, f64> = counts
            .iter()
            .map(|(&(f, t), &(obs, totlen))| {
                let len = (totlen / obs as i64).max(0) as usize;
                let calibed = calib.calib(obs, len);
                ((f, t), calibed)
            })
            .collect();
        let median = {
            let mut covs: Vec<_> = calib_coverage.values().collect();
            let len = covs.len() / MEDIAN_FRAC;
            let (_, median, _) = covs.select_nth_unstable_by(len, |x, y| x.partial_cmp(y).unwrap());
            *median
        };
        debug!("MEDIAN\t{}", median);
        // (edge)->Option<The element to be removed>.
        let to_remove: HashMap<_, u64> = calib_coverage
            .iter()
            .filter(|(_, &cov)| cov < median / 4f64)
            .filter_map(|(&key, &cov)| {
                // Search reads with (from,fd,to,td) occurence.
                let mut former_neighbor = HashSet::new();
                let mut later_neighbor = HashSet::new();
                for read in self.encoded_reads.iter() {
                    for (i, from) in read.nodes.iter().enumerate().take(read.nodes.len() - 1) {
                        let to = &read.nodes[i + 1];
                        if normalize(from, to) == key {
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
                    let mod_cov = *calib_coverage.get(&probe).unwrap_or(&0f64);
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
                    let mod_cov = *calib_coverage.get(&probe).unwrap_or(&0f64);
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
            .collect();
        for read in self.encoded_reads.iter_mut() {
            let remove_idx: Vec<_> = read
                .nodes
                .iter()
                .enumerate()
                .take(read.nodes.len() - 1)
                .filter_map(|(i, from)| {
                    let to = &read.nodes[i + 1];
                    to_remove.get(&normalize(from, to)).map(|&removing| {
                        match removing == from.unit {
                            true => i,
                            false => i + 1,
                        }
                    })
                })
                .collect();
            for (offset, i) in remove_idx.into_iter().enumerate() {
                read.remove(i - offset);
            }
        }
    }
}
