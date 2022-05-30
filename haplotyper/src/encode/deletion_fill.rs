//! Filling deletion.
// TODO:Curretly the alignment can not handle to re-encode missing tips. It reduce some of the reads which should be encoded otherwise...
use crate::ALN_PARAMETER;
use definitions::*;
use rayon::prelude::*;
// use log::*;
use std::collections::{HashMap, HashSet};
// identity would be increased by this value when evaluating the edges.
const EDGE_BOUND: f64 = 0.5;
// Evaluating length of each side.
const EDGE_LEN: usize = 100;
const INS_THR: usize = 2;
// Tuple of unit and cluster.
#[derive(Debug, Clone)]
pub struct CorrectDeletionConfig {
    re_clustering: bool,
    sim_thr: f64,
}

impl CorrectDeletionConfig {
    pub fn new(re_clustering: bool, sim_thr: f64) -> Self {
        Self {
            re_clustering,
            sim_thr,
        }
    }
}

pub trait CorrectDeletion {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig);
}

impl CorrectDeletion for DataSet {
    fn correct_deletion(&mut self, config: &CorrectDeletionConfig) {
        let find_new_units = correct_unit_deletion(self, config.sim_thr);
        if config.re_clustering {
            // Log original assignments.
            let original_assignments = log_original_assignments(self);
            // Log previous copy number
            let prev_copy_numbers: HashMap<u64, _> = self
                .selected_chunks
                .iter()
                .map(|c| (c.id, c.copy_num))
                .collect();
            // Re-estimate copy number
            // Erase cluster
            self.encoded_reads
                .iter_mut()
                .flat_map(|r| r.nodes.iter_mut())
                .for_each(|n| n.cluster = 0);
            use crate::multiplicity_estimation::*;
            let seed = (231043290490.0 * config.sim_thr).round() as u64;
            let config = MultiplicityEstimationConfig::new(1, seed, self.coverage, None);
            self.estimate_multiplicity(&config);
            // Retain all the units changed their copy numbers.
            let chainged_units: HashSet<_> = self
                .selected_chunks
                .iter()
                .filter_map(|c| match prev_copy_numbers.get(&c.id) {
                    None => None,
                    Some(&prev) if c.copy_num == prev => None,
                    Some(_) => Some(c.id),
                })
                .collect();
            // Merge these two.
            let selection: HashSet<_> = find_new_units.union(&chainged_units).copied().collect();
            // Recover the original assignments on the retained units.
            self.encoded_reads
                .iter_mut()
                .zip(original_assignments)
                .for_each(|(read, (id, log))| {
                    assert_eq!(read.id, id);
                    recover_original_assignments(read, &log, &selection);
                });
            // Reclustering.
            crate::local_clustering::local_clustering_selected(self, &selection);
        }
    }
}
// Logging the original assignment into a vector.
fn log_original_assignments(ds: &DataSet) -> Vec<(u64, Vec<(u64, u64)>)> {
    ds.encoded_reads
        .iter()
        .map(|r| {
            let xs: Vec<_> = r.nodes.iter().map(|u| (u.unit, u.cluster)).collect();
            (r.id, xs)
        })
        .collect()
}

// Recover the previous clustering. Note that sometimes the node is added
// so that the length of the read is different from the logged one.
// But it does not removed!
fn recover_original_assignments(read: &mut EncodedRead, log: &[(u64, u64)], except: &HashSet<u64>) {
    let mut read = read.nodes.iter_mut();
    for &(unit, cluster) in log {
        for node in &mut read {
            if node.unit == unit {
                if !except.contains(&node.unit) {
                    node.cluster = cluster;
                }
                break;
            }
        }
    }
}

/**
The second argument is the vector of (index,unit_id) of the previous failed trials.
for example, if failed_trials[i][0] = (j,id), then, we've already tried to encode the id-th unit after the j-th
position of the i-th read, and failed it.
If we can encode some position in the i-th read, the failed trials would be erased, as it change the
condition of the read, making it possible to encode an unit previously failed to encode.
sim_thr is the similarity threshold.
This function corrects "unit-deletions". To do that,
it first align other reads onto a read to be corrected in unit resolution, detecting putative insertions.
Then, it tries to encode these putative insertions in base-pair resolution.
Note that, in the first - unit resolution - alignment, there's no distinction between clusters.
However, in the second alignment, it tries to encode the putative region by each cluster's representative.
Of course, if there's only one cluster for a unit, then, it just tries to encode by that unit.
Auto-tune the similarity threshold.
 */
pub fn correct_unit_deletion(ds: &mut DataSet, fallback: f64) -> HashSet<u64> {
    const OUTER_LOOP: usize = 3;
    let mut find_new_node = HashSet::new();
    for t in 0..OUTER_LOOP {
        remove_weak_edges(ds);
        let (read_error_rate, unit_error_rate, standard_dev) =
            estimate_error_rate_dev(ds, fallback);
        debug!("ErrorRateSTDDev\t{}\t{}", t, standard_dev);
        let (new_nodes, is_updated) =
            filling_until_stable(ds, (&read_error_rate, &unit_error_rate, standard_dev));
        find_new_node.extend(new_nodes);
        if !is_updated {
            break;
        }
    }
    find_new_node
}

pub fn remove_weak_edges(ds: &mut DataSet) {
    let edge_counts = {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for read in ds.encoded_reads.iter() {
            for (e, w) in read.edges.iter().zip(read.nodes.windows(2)) {
                assert_eq!(e.from, w[0].unit);
                assert_eq!(e.to, w[1].unit);
                let entry1 = (w[0].unit, w[0].is_forward, w[1].unit, w[1].is_forward);
                let entry2 = (w[1].unit, !w[1].is_forward, w[0].unit, !w[0].is_forward);
                *counts.entry(entry1).or_default() += 1;
                *counts.entry(entry2).or_default() += 1;
            }
        }
        counts
    };
    let (thr, penalty) = {
        let mut counts: Vec<_> = edge_counts.values().copied().collect();
        let index = counts.len() / 2;
        let median = *counts.select_nth_unstable(index).1;
        debug!("RMEDGE\tMedian\t{median}");
        ((median / 8).max(2), (median / 3).max(1))
    };
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let deleted_nodes: usize = ds
        .encoded_reads
        .iter_mut()
        .map(|read| {
            let len = read.nodes.len();
            match len {
                0..=1 => 0,
                2 => {
                    let (f, t) = (&read.nodes[0], &read.nodes[1]);
                    let (fmat, flen) = f.aln_info(&units[&f.unit]);
                    let (tmat, tlen) = t.aln_info(&units[&t.unit]);
                    let query = (f.unit, f.is_forward, t.unit, t.is_forward);
                    if edge_counts[&query] < thr {
                        let erroneous = (tmat * flen < fmat * tlen) as usize;
                        let unit = read.nodes[erroneous].unit;
                        let count = edge_counts[&query];
                        trace!("RMEDGE\tRemoved\t{unit}\t{len}\t{count}");
                        read.remove(erroneous);
                        1
                    } else {
                        0
                    }
                }
                _ => {
                    let weak_nodes = check_weak_edges(read, penalty, &edge_counts);
                    for &idx in weak_nodes.iter() {
                        let (unit, forward) = (read.nodes[idx].unit, read.nodes[idx].is_forward);
                        let count = {
                            let back = &read.nodes[idx - 1];
                            let back = edge_counts[&(unit, forward, back.unit, back.is_forward)];
                            let front = &read.nodes[idx + 1];
                            let front = edge_counts[&(unit, forward, front.unit, front.is_forward)];
                            back + front
                        };
                        trace!("RMEDGE\tRemoved\t{unit}\t{len}\t{count}",);
                        read.remove(idx);
                    }
                    let mut removed = weak_nodes.len();
                    // Remove the edge nodes.
                    removed += try_remove_head(read, thr, &edge_counts);
                    removed += try_remove_tail(read, thr, &edge_counts);
                    removed
                }
            }
        })
        .sum();
    debug!("RMEDGE\tTotal\t{deleted_nodes}");
}

fn check_weak_edges(
    read: &EncodedRead,
    penalty: u32,
    edge_counts: &HashMap<(u64, bool, u64, bool), u32>,
) -> Vec<usize> {
    assert!(2 < read.nodes.len(), "{}", read.nodes.len());
    let nodes = read.nodes.iter().enumerate();
    let edges: Vec<Vec<_>> = nodes
        .map(|(i, to)| {
            let nodes = read.nodes.iter().enumerate();
            nodes
                .take(i)
                .filter_map(|(j, from)| {
                    let query = (from.unit, from.is_forward, to.unit, to.is_forward);
                    edge_counts.get(&query).map(|&c| (j, c))
                })
                .collect()
        })
        .collect();
    for eds in edges.iter().skip(1) {
        assert!(!eds.is_empty());
    }
    max_chain(read.nodes.len(), &edges, penalty)
}

// Input:Number of nodes, edges between nodes, deletion penalty.
// Output: The indices of the nodes dropped from the input, in reverse order.
// The output would be maximize the penalty of
// retained.windows(2).map(weight).sum() + removed.len() * penalty.
fn max_chain(len: usize, edges: &[Vec<(usize, u32)>], penalty: u32) -> Vec<usize> {
    let penalty = -(penalty as i32);
    assert!(2 < len);
    let mut max_when_in: Vec<(i32, usize)> = Vec::with_capacity(len);
    max_when_in.push((0, 0));
    for (i, edges) in edges.iter().enumerate().skip(1) {
        let next = edges
            .iter()
            .map(|&(from, count)| {
                let total = (i - from) as i32 * penalty + max_when_in[from].0 + count as i32;
                (total, from)
            })
            .max_by_key(|x| x.0)
            .unwrap();
        max_when_in.push(next);
    }
    let mut is_taken = vec![false; len];
    let mut pos = len - 1;
    while 0 < pos {
        is_taken[pos] = true;
        pos = max_when_in[pos].1;
    }
    is_taken[pos] = true;
    is_taken
        .iter()
        .enumerate()
        .filter_map(|(i, is_taken)| (!is_taken).then(|| i))
        .rev()
        .collect()
}

fn try_remove_head(
    read: &mut EncodedRead,
    thr: u32,
    edge_counts: &HashMap<(u64, bool, u64, bool), u32>,
) -> usize {
    let f = &read.nodes[0];
    let t = &read.nodes[1];
    let query = (f.unit, f.is_forward, t.unit, t.is_forward);
    let count = edge_counts[&query];
    if count < thr {
        trace!("RMEDGE\tRemoved\t{}\tH\t{count}", read.nodes[0].unit);
        read.remove(0);
        1
    } else {
        0
    }
}

fn try_remove_tail(
    read: &mut EncodedRead,
    thr: u32,
    edge_counts: &HashMap<(u64, bool, u64, bool), u32>,
) -> usize {
    let len = read.nodes.len();
    let f = &read.nodes[len - 2];
    let t = &read.nodes[len - 1];
    let query = (f.unit, f.is_forward, t.unit, t.is_forward);
    let count = edge_counts[&query];
    if count < thr {
        trace!("RMEDGE\tRemoved\t{}\tT\t{count}", read.nodes[len - 1].unit);
        read.remove(len - 1);
        1
    } else {
        0
    }
}

fn filling_until_stable(
    ds: &mut DataSet,
    (read_error_rate, unit_error_rate, stddev): (&[f64], &[Vec<f64>], f64),
) -> (HashSet<u64>, bool) {
    let raw_seq: HashMap<_, _> = ds
        .raw_reads
        .iter()
        .map(|read| (read.id, read.seq()))
        .collect();
    const INNER_LOOP: usize = 15;
    let mut find_new_node = HashSet::new();
    let mut current: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
    let representative = take_consensus_sequence(ds);
    // i->vector of failed index and units.
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|x| (x.id, x)).collect();
    let mut failed_trials = vec![vec![]; ds.encoded_reads.len()];
    let mut is_updated = vec![true; ds.encoded_reads.len()];
    for i in 0..INNER_LOOP {
        let read_skeltons: Vec<_> = ds.encoded_reads.iter().map(ReadSkelton::new).collect();
        assert_eq!(ds.encoded_reads.len(), failed_trials.len());
        assert_eq!(ds.encoded_reads.len(), is_updated.len());
        let newly_encoded_units: Vec<_> = if log_enabled!(log::Level::Trace) {
            trace!("Tracing. Single thread mode.");
            let reads = ds
                .encoded_reads
                .iter_mut()
                .zip(failed_trials.iter_mut())
                .zip(is_updated.iter_mut())
                .filter(|((r, _), is_updated)| 0 < r.nodes.len() && **is_updated);
            reads
                .flat_map(|((read, fails), is_updated)| {
                    let error_rate = read_error_rate[read.id as usize];
                    let seq = raw_seq[&read.id];
                    let read = (read, seq, error_rate, is_updated);
                    let units = (&units, unit_error_rate, &representative);
                    correct_deletion_error(read, fails, units, stddev, &read_skeltons)
                })
                .collect()
        } else {
            let reads = ds
                .encoded_reads
                .par_iter_mut()
                .zip(failed_trials.par_iter_mut())
                .zip(is_updated.par_iter_mut())
                .filter(|((r, _), is_updated)| 0 < r.nodes.len() && **is_updated);
            reads
                .flat_map(|((read, fails), is_updated)| {
                    let error_rate = read_error_rate[read.id as usize];
                    let seq = raw_seq[&read.id];
                    let read = (read, seq, error_rate, is_updated);
                    let units = (&units, unit_error_rate, &representative);
                    correct_deletion_error(read, fails, units, stddev, &read_skeltons)
                })
                .collect()
        };
        find_new_node.extend(newly_encoded_units);
        let after: usize = ds.encoded_reads.iter().map(|x| x.nodes.len()).sum();
        debug!("Filled:{}\t{}", current, after);
        if after == current && i == 0 {
            debug!("Filled\tBREAK\tOuter");
            return (find_new_node, false);
        } else if after == current {
            debug!("Filled\tBREAK\tInner");
            break;
        }
        current = after;
    }
    ds.encoded_reads.retain(|r| !r.nodes.is_empty());
    (find_new_node, true)
}

// Error Rate of the reads, error rate of the units, and the sqrt of the median of the squared error.
pub fn estimate_error_rate_dev(ds: &DataSet, fallback: f64) -> (Vec<f64>, Vec<Vec<f64>>, f64) {
    type Read = (usize, Vec<(usize, usize, f64)>);
    fn residual(errors: &[Read], reads: &[f64], units: &[Vec<f64>]) -> f64 {
        let residual: f64 = errors
            .iter()
            .map(|&(readid, ref errors)| -> f64 {
                let read = reads[readid];
                let data_error: f64 = errors
                    .iter()
                    .map(|&(unit, cluster, error)| (error - read - units[unit][cluster]).powi(2))
                    .sum();
                data_error
            })
            .sum();
        let reg_term: f64 = units.iter().flatten().map(|x| x * x).sum();
        residual + reg_term
    }
    let max_read_id = ds.raw_reads.iter().map(|r| r.id).max().unwrap() as usize;
    let max_unit_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap() as usize;
    let (mut unit_error_rate, unit_counts) = {
        let mut cluster_num = vec![0; max_unit_id + 1];
        for u in ds.selected_chunks.iter() {
            cluster_num[u.id as usize] = u.cluster_num;
        }
        let unit_errors: Vec<_> = cluster_num.iter().map(|&x| vec![0f64; x]).collect();
        let mut counts: Vec<_> = cluster_num.iter().map(|&x| vec![0; x]).collect();
        for read in ds.encoded_reads.iter() {
            for node in read.nodes.iter() {
                let (unit, cluster) = (node.unit as usize, node.cluster as usize);
                counts[unit][cluster] += 1;
            }
        }
        (unit_errors, counts)
    };
    let mut read_error_rate = vec![0f64; max_read_id + 1];
    for read in ds.encoded_reads.iter() {
        read_error_rate[read.id as usize] = fallback;
    }
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let errors: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|read| {
            let errors: Vec<_> = read
                .nodes
                .iter()
                .map(|node| {
                    let aln = node.recover(units[&node.unit]).1;
                    let errors = aln.iter().filter(|&&x| x != b'|').count();
                    let error = errors as f64 / aln.len() as f64;
                    (node.unit as usize, node.cluster as usize, error)
                })
                .collect::<Vec<_>>();
            (read.id as usize, errors)
        })
        .collect();
    let mut current_resid = residual(&errors, &read_error_rate, &unit_error_rate);
    loop {
        // Re-estimation of unit error rate
        unit_error_rate.iter_mut().flatten().for_each(|x| *x = 0f64);
        for &(readid, ref errors) in errors.iter() {
            let read_error = read_error_rate[readid];
            for &(unit, cluster, error) in errors.iter() {
                unit_error_rate[unit][cluster] += error - read_error;
            }
        }
        for (resid, counts) in unit_error_rate.iter_mut().zip(unit_counts.iter()) {
            for (err, count) in resid.iter_mut().zip(counts) {
                *err = err.max(0f64) / (*count as f64 + 0.1f64);
            }
        }
        // Re-estimate read error rate
        for &(readid, ref errors) in errors.iter() {
            let residual: f64 = errors
                .iter()
                .map(|&(unit, cluster, error)| error - unit_error_rate[unit][cluster])
                .sum();
            read_error_rate[readid] = residual / errors.len() as f64;
            // let len = errors.len() as f64 + 1f64;
            // read_error_rate[readid] = (residual + fallback) / len;
        }
        let resid = residual(&errors, &read_error_rate, &unit_error_rate);
        if (current_resid - resid).abs() < 0.001 {
            break;
        }
        current_resid = resid;
    }
    let mut residuals: Vec<f64> = errors
        .iter()
        .flat_map(|(readid, errors)| {
            let readerror = read_error_rate[*readid as usize];
            errors
                .iter()
                .map(|&(unit, cluster, error)| {
                    let expect = unit_error_rate[unit][cluster] + readerror;
                    error - expect
                })
                .map(|residual| residual.powi(2))
                .collect::<Vec<f64>>()
        })
        .collect();
    let idx = residuals.len() / 2;
    let median = residuals
        .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap())
        .1
        .sqrt();
    (read_error_rate, unit_error_rate, median)
}

pub fn estimate_error_rate(ds: &DataSet, fallback: f64) -> (Vec<f64>, Vec<f64>, f64) {
    fn residual(errors: &[(u64, Vec<(u64, f64)>)], reads: &[f64], units: &[f64]) -> f64 {
        let residual: f64 = errors
            .iter()
            .map(|(readid, errors)| -> f64 {
                let read = reads[*readid as usize];
                let data_error: f64 = errors
                    .iter()
                    .map(|(unitid, error)| (error - read - units[*unitid as usize]).powi(2))
                    .sum();
                data_error
            })
            .sum();
        let reg_term: f64 = units.iter().map(|x| x * x).sum();
        residual + reg_term
    }
    let max_read_id = ds.raw_reads.iter().map(|r| r.id).max().unwrap() as usize;
    let max_unit_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap() as usize;
    let mut read_error_rate = vec![0f64; max_read_id + 1];
    let mut unit_error_rate = vec![0f64; max_unit_id + 1];
    for read in ds.encoded_reads.iter() {
        read_error_rate[read.id as usize] = fallback;
    }
    let unit_counts: Vec<_> = {
        let mut counts = vec![0; max_unit_id + 1];
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            counts[node.unit as usize] += 1;
        }
        counts
    };
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let errors: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|read| {
            let errors: Vec<_> = read
                .nodes
                .iter()
                .map(|node| {
                    let aln = node.recover(units[&node.unit]).1;
                    let errors = aln.iter().filter(|&&x| x != b'|').count();
                    let error = errors as f64 / aln.len() as f64;
                    (node.unit, error)
                })
                .collect::<Vec<_>>();
            (read.id, errors)
        })
        .collect();
    let mut current_resid = residual(&errors, &read_error_rate, &unit_error_rate);
    loop {
        // Re-estimate read error rate
        for (readid, errors) in errors.iter() {
            let residual: f64 = errors
                .iter()
                .map(|&(unitid, error)| error - unit_error_rate[unitid as usize])
                .sum();
            read_error_rate[*readid as usize] = residual / errors.len() as f64;
        }
        // Re-estimation of unit error rate
        unit_error_rate.iter_mut().for_each(|x| *x = 0f64);
        for (readid, errors) in errors.iter() {
            let read_error = read_error_rate[*readid as usize];
            for &(unitid, error) in errors.iter() {
                unit_error_rate[unitid as usize] += error - read_error;
            }
        }
        unit_error_rate
            .iter_mut()
            .zip(unit_counts.iter())
            .filter(|&(_, &count)| 0 < count)
            .for_each(|(x, c)| {
                *x = *x / (*c as f64 + 1f64);
            });
        let resid = residual(&errors, &read_error_rate, &unit_error_rate);
        if (current_resid - resid).abs() < 0.001 {
            break;
        }
        current_resid = resid;
    }
    let mut residuals: Vec<f64> = errors
        .iter()
        .flat_map(|(readid, errors)| {
            let readerror = read_error_rate[*readid as usize];
            errors
                .iter()
                .map(|(unitid, error)| {
                    let expect = unit_error_rate[*unitid as usize] + readerror;
                    error - expect
                })
                .map(|residual| residual.powi(2))
                .collect::<Vec<f64>>()
        })
        .collect();
    let idx = residuals.len() / 2;
    let median = residuals
        .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap())
        .1
        .sqrt();
    // Fix read error rate
    for (&readid, len) in errors.iter().map(|(id, es)| (id, es.len() as f64)) {
        let error = read_error_rate[readid as usize];
        read_error_rate[readid as usize] = (error * len + fallback) / (len + 1f64);
    }
    (read_error_rate, unit_error_rate, median)
}

// fn estimate_upper_error_rate(ds: &DataSet, fallback: f64) -> (Vec<f64>, Vec<f64>, f64) {
//     let ref_chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
//     let error_rates: Vec<Vec<_>> = ds
//         .encoded_reads
//         .par_iter()
//         .map(|r| {
//             r.nodes
//                 .iter()
//                 .map(|n| {
//                     let (_, aln, _) = n.recover(ref_chunks[&n.unit]);
//                     let error = aln.iter().filter(|&&b| b != b'|').count();
//                     error as f64 / aln.len() as f64
//                 })
//                 .collect()
//         })
//         .collect();
//     let (sd_sum, sd_num) = error_rates
//         .iter()
//         .filter(|errs| errs.len() > 2)
//         .map(|errs| {
//             let (sum, sumsq) = errs
//                 .iter()
//                 .fold((0f64, 0f64), |(sum, sumsq), x| (sum + x, sumsq + x * x));
//             let mean = sum / errs.len() as f64;
//             let var = sumsq / errs.len() as f64 - mean * mean;
//             assert!(0f64 <= var);
//             var.sqrt()
//         })
//         .fold((0f64, 0), |(sdsum, num), sd| (sdsum + sd, num + 1));
//     let sd_mean = sd_sum / sd_num as f64;
//     debug!("MEAN of SD\t{:.3}", sd_mean);
//     let error_rate: Vec<_> = error_rates
//         .into_par_iter()
//         .map(|mut errors| match errors.len() {
//             0..=2 => fallback,
//             x => {
//                 let median = errors
//                     .select_nth_unstable_by(x / 2, |x, y| x.partial_cmp(y).unwrap())
//                     .1;
//                 *median + 4f64 * sd_mean
//             }
//         })
//         .collect();
//     (error_rate, vec![], 0f64)
// }

// Take consensus of each cluster of each unit, return the consensus seuqneces.
// UnitID->(clsuterID, its consensus).
fn take_consensus_sequence(ds: &DataSet) -> HashMap<u64, Vec<(u64, Vec<u8>)>> {
    fn polish(xs: &[&[u8]], unit: &Unit, band: usize) -> Vec<u8> {
        kiley::bialignment::guided::polish_until_converge(unit.seq(), xs, band)
    }
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|u| (u.id, u)).collect();
    let mut bucket: HashMap<u64, Vec<_>> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        bucket
            .entry(node.unit)
            .or_default()
            .push((node.cluster, node.seq()));
    }
    bucket
        .par_iter()
        .filter(|&(_, xs)| !xs.is_empty())
        .map(|(&unit, bucket)| {
            let ref_unit = &ref_units[&unit];
            let mut clusters: HashMap<_, Vec<&[u8]>> = HashMap::new();
            for (cl, seq) in bucket {
                clusters.entry(*cl).or_default().push(seq);
            }
            let band = ds.read_type.band_width(ref_unit.seq().len());
            assert!(!clusters.is_empty());
            let representative: Vec<_> = if clusters.len() == 1 {
                clusters
                    .iter()
                    .map(|(&cl, _)| (cl, ref_unit.seq().to_vec()))
                    .collect()
            } else {
                clusters
                    .iter()
                    .filter(|(_, xs)| !xs.is_empty())
                    .map(|(&cl, xs)| (cl, polish(xs, ref_unit, band)))
                    .collect()
            };
            (unit, representative)
        })
        .collect()
}

#[inline]
fn abs(x: usize, y: usize) -> usize {
    x.max(y) - x.min(y)
}

// Aligment offset. We align [s-offset..e+offset] region to the unit.
// const OFFSET: usize = 150;
const OFFSET_FACTOR: f64 = 0.1;
// returns the ids of the units newly encoded.
// Maybe each (unit,cluster) should corresponds to a key...?
type UnitInfo<'a> = (
    &'a HashMap<u64, &'a Unit>,
    &'a [Vec<f64>],
    &'a HashMap<u64, Vec<(u64, Vec<u8>)>>,
);
pub fn correct_deletion_error(
    (read, seq, read_error, is_changed): (&mut EncodedRead, &[u8], f64, &mut bool),
    failed_trials: &mut Vec<(usize, LightNode)>,
    unitinfo: UnitInfo,
    stddev: f64,
    reads: &[ReadSkelton],
) -> Vec<u64> {
    let pileups = get_pileup(read, reads);
    let nodes = &read.nodes;
    let mut inserts = vec![];
    let ins_thr = INS_THR.min(nodes.len());
    for (idx, pileup) in pileups.iter().enumerate() {
        let mut head_cand = pileup.check_insertion_head(nodes, ins_thr, idx);
        head_cand.retain(|node, _| !failed_trials.contains(&(idx, *node)));
        let head_best =
            try_encoding_head(nodes, &head_cand, idx, unitinfo, seq, read_error, stddev);
        match head_best {
            Some((head_node, _)) => inserts.push((idx, head_node)),
            None => failed_trials.extend(head_cand.into_iter().map(|(n, _)| (idx, n))),
        }
        let mut tail_cand = pileup.check_insertion_tail(nodes, ins_thr, idx);
        tail_cand.retain(|node, _| !failed_trials.contains(&(idx, *node)));
        let tail_best =
            try_encoding_tail(nodes, &tail_cand, idx, unitinfo, seq, read_error, stddev);
        match tail_best {
            Some((tail_node, _)) => inserts.push((idx, tail_node)),
            None => failed_trials.extend(tail_cand.into_iter().map(|x| (idx, x.0))),
        }
    }
    *is_changed = !inserts.is_empty();
    let new_inserts: Vec<_> = inserts.iter().map(|(_, n)| n.unit).collect();
    if !inserts.is_empty() {
        failed_trials.clear();
        for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
            read.nodes.insert(idx + accum_inserts, node);
        }
    }
    *is_changed |= remove_highly_erroneous(read, read_error, unitinfo, stddev);
    if *is_changed && !read.nodes.is_empty() {
        let mut nodes = Vec::with_capacity(read.nodes.len());
        nodes.append(&mut read.nodes);
        use super::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| (n.unit, n.position_from_start));
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
    new_inserts
}

fn remove_highly_erroneous(
    read: &mut EncodedRead,
    read_error: f64,
    (units, unit_error_rate, _): UnitInfo,
    stddev: f64,
) -> bool {
    let orig_len = read.nodes.len();
    read.nodes.retain(|node| {
        let (_, aln, _) = node.recover(units[&node.unit]);
        let diff = aln.iter().filter(|&&x| x != b'|').count();
        let error_rate = diff as f64 / aln.len() as f64;
        let expected = read_error + unit_error_rate[node.unit as usize][node.cluster as usize];
        let threshold = (expected + THR * stddev).max(0f64);
        // if threshold < error_rate {
        //     let unit_error = unit_error_rate[node.unit as usize][node.cluster as usize];
        //     let (unit,cluster) = (node.unit,node.cluster);
        //     let (xr,ar,yr) = node.recover(&units[&node.unit]);
        //     for ((xr,ar),yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)){
        //         eprintln!("{}",std::str::from_utf8(xr).unwrap());
        //         eprintln!("{}",std::str::from_utf8(ar).unwrap());
        //         eprintln!("{}\n",std::str::from_utf8(yr).unwrap());
        //     }
        // }
        error_rate < threshold
    });
    read.nodes.len() != orig_len
}

const THR: f64 = 9f64;
fn try_encoding_head(
    nodes: &[Node],
    head_cand: &HashMap<LightNode, usize>,
    idx: usize,
    (units, unit_error_rate, consensi): UnitInfo,
    seq: &[u8],
    read_error: f64,
    stddev: f64,
) -> Option<(Node, i32)> {
    head_cand
        .iter()
        .filter(|&(_, &start_position)| start_position < seq.len())
        .filter_map(|(node, &start_position)| {
            let (uid, cluster) = (node.unit, node.cluster);
            let unit = *units.get(&uid)?;
            let (_, cons) = consensi.get(&uid)?.iter().find(|&&(cl, _)| cl == cluster)?;
            let offset = (OFFSET_FACTOR * cons.len() as f64).ceil() as usize;
            let end_position = (start_position + cons.len() + offset).min(seq.len());
            let start_position = start_position.saturating_sub(offset);
            // let end_position = (start_position + cons.len() + 2 * OFFSET).min(seq.len());
            let is_the_same_encode = match nodes.get(idx) {
                Some(node) => {
                    node.unit == uid && abs(node.position_from_start, start_position) < cons.len()
                }
                None => false,
            };
            assert!(start_position < end_position);
            if is_the_same_encode {
                return None;
            }
            let expected = read_error + unit_error_rate[uid as usize][cluster as usize];
            let error_rate_bound = expected + THR * stddev;
            let position = (start_position, end_position, node.is_forward);
            let unit_info = (unit, cluster, cons.as_slice());
            encode_node(seq, position, unit_info, error_rate_bound)
        })
        .max_by_key(|x| x.1)
}

fn try_encoding_tail(
    nodes: &[Node],
    tail_cand: &HashMap<LightNode, usize>,
    idx: usize,
    (units, unit_error_rate, consensi): UnitInfo,
    seq: &[u8],
    read_error: f64,
    stddev: f64,
) -> Option<(Node, i32)> {
    tail_cand
        .iter()
        .filter_map(|(node, &end_position)| {
            let (uid, cluster) = (node.unit, node.cluster);
            let unit = *units.get(&uid)?;
            let (_, cons) = consensi.get(&uid)?.iter().find(|&&(cl, _)| cl == cluster)?;
            let offset = (OFFSET_FACTOR * cons.len() as f64).ceil() as usize;
            let start_position = end_position.min(seq.len()).saturating_sub(offset);
            let end_position = (end_position + offset).min(seq.len());
            // let end_position = end_position.min(seq.len());
            // let start_position = end_position.saturating_sub(cons.len() + 2 * OFFSET);
            assert!(start_position < end_position);
            let is_the_same_encode = match nodes.get(idx) {
                Some(node) => {
                    node.unit == uid && abs(node.position_from_start, start_position) < cons.len()
                }
                None => false,
            };
            if is_the_same_encode {
                return None;
            }
            assert!(start_position < end_position);
            let positions = (start_position, end_position, node.is_forward);
            let unit_info = (unit, cluster, cons.as_slice());
            let expected = read_error + unit_error_rate[uid as usize][cluster as usize];
            let error_rate_bound = expected + THR * stddev;
            encode_node(seq, positions, unit_info, error_rate_bound)
        })
        .max_by_key(|x| x.1)
}

// Try to Encode Node. Return Some(node) if the alignment is good.
// Return also the alignment score of the encoding.
// The match score is 2, mism is -6, gap open is -5, and gap ext is -1.
fn encode_node(
    query: &[u8],
    (start, end, is_forward): (usize, usize, bool),
    (unit, cluster, unitseq): (&Unit, u64, &[u8]),
    sim_thr: f64,
) -> Option<(Node, i32)> {
    // Initial filter.
    // If the query is shorter than the unitseq,
    // at least we need the edit operations to fill the gaps.
    // This is lower bound of the sequence identity.
    let edit_dist_lower_bound = unitseq.len().saturating_sub(end - start);
    let diff_lower_bound = edit_dist_lower_bound as f64 / unitseq.len() as f64;
    if sim_thr < diff_lower_bound {
        return None;
    }
    // Tune the query...
    let mut query = if is_forward {
        query[start..end].to_vec()
    } else {
        bio_utils::revcmp(&query[start..end])
    };
    query.iter_mut().for_each(u8::make_ascii_uppercase);
    let (seq, trim_head, trim_tail, kops, score) =
        fine_mapping(&query, (unit, cluster, unitseq), sim_thr)?;
    let ops = super::compress_kiley_ops(&kops);
    let cl = unit.cluster_num;
    let position_from_start = match is_forward {
        true => start + trim_head,
        false => start + trim_tail,
    };
    // I think we should NOT make likelihood gain to some biased value,
    // as 1. if the alignment gives the certaintly, then we can impute the clustering by the alignment,
    // 2. if `cluster` assignment is just by chance,
    // then we just should not introduce any bias into the likelihood gain.
    let mut node = Node::new(unit.id, is_forward, seq, ops, position_from_start, cl);
    node.cluster = cluster;
    Some((node, score))
}

// Sequence, trimed base from the head, trimed base from the tail, ops, score.
type FineMapping<'a> = (&'a [u8], usize, usize, Vec<kiley::Op>, i32);
fn fine_mapping<'a>(
    orig_query: &'a [u8],
    (unit, cluster, unitseq): (&Unit, u64, &[u8]),
    sim_thr: f64,
) -> Option<FineMapping<'a>> {
    use kiley::bialignment::guided::infix_guided;
    fn edlib_op_to_kiley_op(ops: &[u8]) -> Vec<kiley::Op> {
        use kiley::Op::*;
        ops.iter()
            .map(|&op| [Match, Ins, Del, Mismatch][op as usize])
            .collect()
    }
    let (query, trim_head, trim_tail, ops, band) = {
        // TODO:Is this correct?
        let mode = edlib_sys::AlignMode::Global;
        let task = edlib_sys::AlignTask::Alignment;
        let alignment = edlib_sys::edlib_align(unitseq, orig_query, mode, task);
        let band = ((orig_query.len() as f64 * sim_thr * 0.3).ceil() as usize).max(10);
        let ops = edlib_op_to_kiley_op(&alignment.operations.unwrap());
        // Align twice, to get an accurate alignment.
        let (_, ops) = infix_guided(orig_query, unitseq, &ops, band, ALN_PARAMETER);
        let (_, mut ops) = infix_guided(orig_query, unitseq, &ops, band, ALN_PARAMETER);
        // Reverse ops
        for op in ops.iter_mut() {
            *op = match *op {
                kiley::Op::Ins => kiley::Op::Del,
                kiley::Op::Del => kiley::Op::Ins,
                x => x,
            }
        }
        let (trim_head, trim_tail) = trim_head_tail_insertion(&mut ops);
        let query = &orig_query[trim_head..orig_query.len() - trim_tail];
        (query, trim_head, trim_tail, ops, band)
    };
    let (below_dissim, info) = {
        let mat_num = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
        let identity = mat_num as f64 / ops.len() as f64;
        let (head_identity, tail_identity) = edge_identity(unitseq, query, &ops, EDGE_LEN);
        let (rlen, qlen) = (unitseq.len(), query.len());
        let id = unit.id;
        let orig_len = orig_query.len();
        let info = format!("{id}\t{cluster}\t{identity:.2}\t{rlen}\t{qlen}\t{orig_len}");
        let iden_bound = 1f64 - sim_thr;
        let is_ok = iden_bound < identity && EDGE_BOUND < head_identity.min(tail_identity);
        (is_ok, info)
    };
    if log_enabled!(log::Level::Trace) {
        if below_dissim {
            trace!("FILLDEL\t{}\tOK", info);
        } else {
            trace!("FILLDEL\t{}\tNG", info);
            if log_enabled!(log::Level::Trace) {
                let (xr, ar, yr) = kiley::recover(unitseq, query, &ops);
                for ((xr, ar), yr) in xr.chunks(200).zip(ar.chunks(200)).zip(yr.chunks(200)) {
                    eprintln!("ALN\t{}", String::from_utf8_lossy(xr));
                    eprintln!("ALN\t{}", String::from_utf8_lossy(ar));
                    eprintln!("ALN\t{}", String::from_utf8_lossy(yr));
                }
            }
        }
    }
    (below_dissim).then(|| {
        let (score, new_ops) = infix_guided(unit.seq(), query, &ops, band, ALN_PARAMETER);
        (query, trim_head, trim_tail, new_ops, score)
    })
}

#[allow(dead_code)]
fn edge_identity(unit: &[u8], _: &[u8], ops: &[kiley::Op], len: usize) -> (f64, f64) {
    let (mut head_aln_len, mut head_match) = (0, 0);
    let (mut tail_aln_len, mut tail_match) = (0, 0);
    let head_eval_end = len.min(unit.len());
    let tail_eval_start = unit.len().saturating_sub(len);
    let mut rpos = 0;
    for &op in ops {
        match op {
            kiley::Op::Mismatch | kiley::Op::Match => rpos += 1,
            kiley::Op::Ins => {}
            kiley::Op::Del => rpos += 1,
        }
        if rpos < head_eval_end {
            head_aln_len += 1;
            head_match += (kiley::Op::Match == op) as usize;
        }
        if tail_eval_start <= rpos {
            tail_aln_len += 1;
            tail_match += (kiley::Op::Match == op) as usize;
        }
    }
    let head_identity = head_match as f64 / head_aln_len as f64;
    let tail_identity = tail_match as f64 / tail_aln_len as f64;
    (head_identity, tail_identity)
}

// Triming the head/tail insertion, re-calculate the start and end position.
fn trim_head_tail_insertion(ops: &mut Vec<kiley::Op>) -> (usize, usize) {
    let mut tail_ins = 0;
    while ops.last() == Some(&kiley::Op::Ins) {
        ops.pop();
        tail_ins += 1;
    }
    ops.reverse();
    let mut head_ins = 0;
    while ops.last() == Some(&kiley::Op::Ins) {
        ops.pop();
        head_ins += 1;
    }
    ops.reverse();
    (head_ins, tail_ins)
}

fn check_alignment_by_unitmatch(units: &[(u64, u64, bool)], query: &ReadSkelton) -> Option<bool> {
    fn count_match(units: &[(u64, u64, bool)], query: &[(u64, u64, bool)]) -> usize {
        let mut r_ptr = units.iter().peekable();
        let mut q_ptr = query.iter().peekable();
        let mut match_num = 0;
        while r_ptr.peek().is_some() && q_ptr.peek().is_some() {
            match r_ptr.peek().unwrap().cmp(q_ptr.peek().unwrap()) {
                std::cmp::Ordering::Less => r_ptr.next(),
                std::cmp::Ordering::Equal => {
                    match_num += 1;
                    r_ptr.next();
                    q_ptr.next()
                }
                std::cmp::Ordering::Greater => q_ptr.next(),
            };
        }
        match_num
    }
    let mut keys: Vec<_> = query.nodes.iter().map(|n| n.key()).collect();
    keys.sort_unstable();
    let forward_match = count_match(units, &keys);
    keys.iter_mut().for_each(|x| x.2 = !x.2);
    keys.sort_unstable();
    let reverse_match = count_match(units, &keys);
    let min_match = MIN_MATCH.min(units.len());
    (min_match <= forward_match.max(reverse_match)).then(|| reverse_match <= forward_match)
}

// Align read skeltons to read, return the pileup sumamries.
// i-> insertions before the i-th nodes.
// The coverage of the last slot is always zero.
pub fn get_pileup(read: &EncodedRead, reads: &[ReadSkelton]) -> Vec<Pileup> {
    assert!(!read.nodes.is_empty());
    let mut pileups = vec![Pileup::new(); read.nodes.len() + 1];
    let skelton = ReadSkelton::new(read);
    let mut units_in_read: Vec<_> = skelton.nodes.iter().map(|n| n.key()).collect();
    units_in_read.sort_unstable();
    for query in reads.iter() {
        let is_forward = match check_alignment_by_unitmatch(&units_in_read, query) {
            Some(is_forward) => is_forward,
            None => continue,
        };
        let id = read.id;
        let aln = match alignment(id, &skelton, query, is_forward) {
            Some(res) => res,
            None => continue,
        };
        let mut q_ptr = SkeltonIter::new(query, is_forward);
        let mut pileups = pileups.iter_mut();
        let mut current_pu = pileups.next().unwrap();
        let mut position = 0;
        for op in aln {
            // These unwraps are safe.
            match op {
                Op::Ins(l) if position == 0 => {
                    // Retain only the last insertion...
                    current_pu.add_tail(q_ptr.nth(l - 1).unwrap());
                }
                Op::Ins(l) if position == pileups.len() - 1 => {
                    // Retain only the first insertion...
                    current_pu.add_head(q_ptr.next().unwrap());
                    for _ in 0..l - 1 {
                        q_ptr.next().unwrap();
                    }
                }
                Op::Ins(l) => {
                    current_pu.add_head(q_ptr.next().unwrap());
                    for _ in 0..l - 1 {
                        current_pu.add_tail(q_ptr.next().unwrap());
                    }
                }
                Op::Del(l) => {
                    current_pu = pileups.nth(l - 1).unwrap();
                    position += l;
                }
                Op::Match(l) => {
                    q_ptr.nth(l - 1);
                    position += l;
                    for _ in 0..l {
                        current_pu.coverage += 1;
                        current_pu = pileups.next().unwrap();
                    }
                }
            }
        }
    }
    pileups
}

// Maybe we should tune this.
// For example, is it ok to use these parameters to treat long repeats?
// Maybe OK, as we confirm these candidate by alignment.
// Minimum required units to be matched.
const MIN_MATCH: usize = 2;
// Minimum required alignment score.
const SCORE_THR: i32 = 1;
fn alignment(_: u64, read: &ReadSkelton, query: &ReadSkelton, dir: bool) -> Option<Vec<Op>> {
    // let mut query = query.clone();
    let (score, ops) = match dir {
        true => pairwise_alignment_gotoh(read, query),
        false => {
            let query = query.rev();
            pairwise_alignment_gotoh(read, &query)
        }
    };
    let match_num = get_match_units(&ops);
    let min_match = MIN_MATCH.min(read.nodes.len()).min(query.nodes.len());
    // let que: Vec<_> = query
    //     .nodes
    //     .iter()
    //     .map(|x| format!("{}-{}", x.unit, x.cluster))
    //     .collect();
    // println!(
    //     "QUE\t{match_num}\t{score}\t{dir}\t{}\t{}",
    //     query.id,
    //     que.join("\t")
    // );
    (min_match <= match_num && SCORE_THR <= score && is_proper(&ops)).then(|| ops)
}

// Return true if the alignment is proper dovetail.
fn is_proper(ops: &[Op]) -> bool {
    ops.windows(2)
        .all(|xs| !matches!(xs, &[Op::Ins(_), Op::Del(_)] | &[Op::Del(_), Op::Ins(_)]))
}

const MIN_ALN: i32 = -10000000;
fn score(x: &LightNode, y: &LightNode) -> i32 {
    if x.unit != y.unit || x.is_forward != y.is_forward {
        MIN_ALN
    } else if x.cluster == y.cluster {
        1
    } else {
        -1
    }
}

fn pairwise_alignment_gotoh(read: &ReadSkelton, query: &ReadSkelton) -> (i32, Vec<Op>) {
    let (read, query) = (&read.nodes, &query.nodes);
    let (row_num, col_num) = (read.len() + 1, query.len() + 1);
    let mut dp = vec![0; row_num * col_num * 3];
    let read_row = col_num * 3;
    // Initialize.
    for i in 0..read.len() + 1 {
        dp[read_row * i] = MIN_ALN;
        dp[read_row * i + 1] = MIN_ALN;
    }
    for j in 0..query.len() + 1 {
        dp[3 * j] = MIN_ALN;
        dp[3 * j + 2] = MIN_ALN;
    }
    dp[0] = 0;
    // Filling DP Table.
    for (i, x) in read.iter().enumerate() {
        for (j, y) in query.iter().enumerate() {
            let (i, j) = (i + 1, j + 1);
            let fill_pos = read_row * i + 3 * j;
            let prev_match = read_row * (i - 1) + 3 * (j - 1);
            dp[fill_pos] = dp[prev_match..prev_match + 3].iter().max().unwrap() + score(x, y);
            let prev_ins = read_row * i + 3 * (j - 1);
            dp[fill_pos + 1] = (dp[prev_ins] - 1).max(dp[prev_ins + 1]);
            let prev_del = read_row * (i - 1) + 3 * j;
            dp[fill_pos + 2] = (dp[prev_del] - 1).max(dp[prev_del + 2]);
        }
    }
    let (mut r_pos, mut q_pos, mut state, dist) = (0..read.len() + 1)
        .map(|i| (i, query.len()))
        .chain((0..query.len() + 1).map(|j| (read.len(), j)))
        .filter_map(|(i, j)| {
            let position = read_row * i + 3 * j;
            dp[position..position + 3]
                .iter()
                .enumerate()
                .max_by_key(|x| x.1)
                .map(|(state, &score)| (i, j, state, score))
        })
        .max_by_key(|x| x.3)
        .unwrap();
    let mut ops = Vec::with_capacity(row_num.max(col_num) + 2);
    if read.len() != r_pos {
        ops.push(Op::Del(read.len() - r_pos));
    }
    if query.len() != q_pos {
        ops.push(Op::Ins(query.len() - q_pos));
    }
    while 0 < r_pos && 0 < q_pos {
        let current_pos = read_row * r_pos + 3 * q_pos + state;
        let current_dist = dp[current_pos];
        if state == 0 {
            let dist = current_dist - score(&read[r_pos - 1], &query[q_pos - 1]);
            let prev_pos = read_row * (r_pos - 1) + 3 * (q_pos - 1);
            let (new_state, _) = dp[prev_pos..prev_pos + 3]
                .iter()
                .enumerate()
                .find(|&(_, &score)| score == dist)
                .unwrap();
            state = new_state;
            ops.push(Op::Match(1));
            r_pos -= 1;
            q_pos -= 1;
        } else if state == 1 {
            let prev_pos = read_row * r_pos + 3 * (q_pos - 1);
            state = if current_dist == dp[prev_pos] - 1 {
                0
            } else {
                //assert_eq!(current_dist, dp[prev_pos + 1]);
                1
            };
            ops.push(Op::Ins(1));
            q_pos -= 1;
        } else {
            let prev_pos = read_row * (r_pos - 1) + 3 * q_pos;
            state = if current_dist == dp[prev_pos] - 1 {
                0
            } else {
                // assert_eq!(current_dist, dp[prev_pos + 2]);
                2
            };
            ops.push(Op::Del(1));
            r_pos -= 1;
        }
    }
    assert!(r_pos == 0 || q_pos == 0);
    if r_pos != 0 {
        ops.push(Op::Del(r_pos));
    }
    if q_pos != 0 {
        ops.push(Op::Ins(q_pos));
    }
    ops.reverse();
    let ops = compress_operations(ops);
    (dist, ops)
}

fn compress_operations(ops: Vec<Op>) -> Vec<Op> {
    assert!(!ops.is_empty());
    let mut current_op = ops[0];
    let mut compressed = Vec::with_capacity(ops.len());
    for &op in ops.iter().skip(1) {
        match (op, current_op) {
            (Op::Match(l), Op::Match(m)) => {
                current_op = Op::Match(l + m);
            }
            (Op::Ins(l), Op::Ins(m)) => {
                current_op = Op::Ins(l + m);
            }
            (Op::Del(l), Op::Del(m)) => {
                current_op = Op::Del(l + m);
            }
            (x, _) => {
                compressed.push(current_op);
                current_op = x;
            }
        }
    }
    compressed.push(current_op);
    compressed
}

fn get_match_units(ops: &[Op]) -> usize {
    ops.iter()
        .map(|op| match op {
            Op::Match(l) => *l,
            _ => 0,
        })
        .sum::<usize>()
}

#[derive(Debug, Clone)]
pub struct Pileup {
    // insertion at the beggining of this node
    head_inserted: Vec<LightNode>,
    // insertion at the last of this node
    tail_inserted: Vec<LightNode>,
    coverage: usize,
}

// TODO: Maybe we should care abount cluster...?
impl Pileup {
    // Return the maximum insertion from the same unit, the same direction.
    fn insertion_head(&self) -> HashMap<LightNode, usize> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.head_inserted.iter() {
            *count.entry(node.clone()).or_default() += 1;
        }
        count
    }
    fn insertion_tail(&self) -> HashMap<LightNode, usize> {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in self.tail_inserted.iter() {
            *count.entry(node.clone()).or_default() += 1;
        }
        count
    }
    fn new() -> Self {
        Self {
            head_inserted: Vec::with_capacity(5),
            tail_inserted: Vec::with_capacity(5),
            coverage: 0,
        }
    }
    fn information_head(&self, node: &LightNode) -> (Option<isize>, Option<isize>) {
        Self::summarize(&self.head_inserted, node)
    }
    fn information_tail(&self, node: &LightNode) -> (Option<isize>, Option<isize>) {
        Self::summarize(&self.tail_inserted, node)
    }
    fn summarize<'a>(inserts: &[LightNode], target: &LightNode) -> (Option<isize>, Option<isize>) {
        let inserts = inserts.iter().filter(|&node| node == target);
        let (mut prev_count, mut prev_total) = (0, 0);
        let (mut after_count, mut after_total) = (0, 0);
        for node in inserts {
            if let Some(len) = node.prev_offset {
                prev_count += 1;
                prev_total += len;
            }
            if let Some(len) = node.after_offset {
                after_count += 1;
                after_total += len;
            }
        }
        let prev_offset = (prev_count != 0).then(|| prev_total / prev_count);
        let after_offset = (after_count != 0).then(|| after_total / after_count);
        (prev_offset, after_offset)
    }
    fn add_head(&mut self, node: LightNode) {
        self.head_inserted.push(node);
    }
    fn add_tail(&mut self, node: LightNode) {
        self.tail_inserted.push(node);
    }
    pub fn check_insertion_head(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> HashMap<LightNode, usize> {
        let mut inserts = self.insertion_head();
        inserts.retain(|_, num| threshold <= *num);
        inserts.retain(|node, num| {
            let (prev_offset, _) = self.information_head(&node);
            let start_position = nodes[idx - 1].position_from_start + nodes[idx - 1].query_length();
            match prev_offset {
                Some(x) => {
                    *num = (start_position as isize + x) as usize;
                    true
                }
                None => false,
            }
        });
        inserts
    }
    pub fn check_insertion_tail(
        &self,
        nodes: &[Node],
        threshold: usize,
        idx: usize,
    ) -> HashMap<LightNode, usize> {
        let end_position = match nodes.get(idx) {
            Some(res) => res.position_from_start as isize,
            None => return HashMap::new(),
        };
        let mut inserts = self.insertion_tail();
        inserts.retain(|_, num| threshold <= *num);
        inserts.retain(|node, num| {
            let (_, after_offset) = self.information_tail(node);
            match after_offset {
                Some(x) => {
                    *num = (end_position - x).max(0) as usize;
                    true
                }
                None => false,
            }
        });
        inserts
    }
}

#[derive(Clone)]
pub struct ReadSkelton {
    id: u64,
    nodes: Vec<LightNode>,
}

impl std::fmt::Debug for ReadSkelton {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for node in self.nodes.iter() {
            write!(f, "{:?}:", node)?;
        }
        Ok(())
    }
}

impl ReadSkelton {
    pub fn new(read: &EncodedRead) -> Self {
        Self::from_rich_nodes(read.id, &read.nodes)
    }
    pub fn from_rich_nodes(id: u64, nodes: &[Node]) -> Self {
        // Convert the nodes into (start_position, end_position)s
        let summaries: Vec<_> = nodes
            .iter()
            .map(|node| {
                let start = node.position_from_start;
                let end = start + node.query_length();
                (start as isize, end as isize)
            })
            .collect();
        let nodes: Vec<_> = nodes
            .iter()
            .enumerate()
            .map(|(i, n)| {
                let prev_end = if i == 0 { None } else { summaries.get(i - 1) };
                let prev_offset = prev_end.map(|x| summaries[i].0 as isize - x.1 as isize);
                let after_offset = summaries.get(i + 1).map(|x| x.0 - summaries[i].1);
                LightNode {
                    prev_offset,
                    unit: n.unit,
                    cluster: n.cluster,
                    is_forward: n.is_forward,
                    after_offset,
                }
            })
            .collect();
        // let sets: HashSet<_> = nodes.iter().map(|n| n.unit).collect();
        ReadSkelton { id, nodes }
    }
    fn rev(&self) -> Self {
        let nodes: Vec<_> = self.nodes.iter().rev().map(LightNode::rev).collect();
        Self { id: self.id, nodes }
    }
}

#[derive(Clone, Copy)]
pub struct LightNode {
    // How long should be add to the last position of the previous node to
    // get the start position of this node.
    // None if this is the first node.
    prev_offset: Option<isize>,
    unit: u64,
    cluster: u64,
    is_forward: bool,
    // Almost the same as prev_offset. The distance between the last postion of this node to
    // the start position of the next node.
    // None if this is the last node.
    after_offset: Option<isize>,
}

impl std::cmp::PartialEq for LightNode {
    fn eq(&self, other: &Self) -> bool {
        self.unit == other.unit
            && self.cluster == other.cluster
            && self.is_forward == other.is_forward
    }
}

impl std::cmp::Eq for LightNode {}

impl std::hash::Hash for LightNode {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.unit.hash(state);
        self.cluster.hash(state);
        self.is_forward.hash(state);
    }
}

impl std::fmt::Debug for LightNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (u, c) = (self.unit, self.cluster);
        let dir = if self.is_forward { '+' } else { '-' };
        write!(f, "{u}-{c}({dir})")
    }
}

impl LightNode {
    fn key(&self) -> (u64, u64, bool) {
        (self.unit, self.cluster, self.is_forward)
    }
    fn rev(
        &Self {
            prev_offset,
            unit,
            cluster,
            is_forward,
            after_offset,
        }: &Self,
    ) -> Self {
        Self {
            prev_offset: after_offset,
            unit,
            cluster,
            is_forward: !is_forward,
            after_offset: prev_offset,
        }
    }
}

struct SkeltonIter<'a> {
    inner: &'a ReadSkelton,
    index: usize,
    is_forward: bool,
}

impl<'a> SkeltonIter<'a> {
    fn new(read: &'a ReadSkelton, is_forward: bool) -> Self {
        let mut it = Self {
            inner: read,
            is_forward,
            index: 0,
        };
        if !is_forward {
            it.index = read.nodes.len();
        }
        it
    }
}

impl<'a> std::iter::Iterator for SkeltonIter<'a> {
    type Item = LightNode;
    fn next(&mut self) -> Option<Self::Item> {
        if self.is_forward {
            self.index += 1;
            self.inner.nodes.get(self.index - 1).cloned()
        } else if 0 < self.index {
            self.index -= 1;
            self.inner.nodes.get(self.index).map(LightNode::rev)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod deletion_fill {
    use super::*;
    #[test]
    fn aln_test_gotoh() {
        let into_reads = |nodes: Vec<u64>| {
            let nodes: Vec<_> = nodes
                .into_iter()
                .map(|unit| LightNode {
                    prev_offset: None,
                    unit,
                    cluster: 0,
                    is_forward: true,
                    after_offset: None,
                })
                .collect();
            ReadSkelton { id: 0, nodes }
        };
        let read = into_reads(vec![69, 148, 318, 0]);
        let query = into_reads(vec![69, 221, 286, 148, 318]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 2, "{:?}", ops);
        let read = into_reads(vec![0]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 1, "{:?}", ops);
        let read = into_reads(vec![0, 4]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 1, "{:?}", ops);
        let (score, ops) = pairwise_alignment_gotoh(&query, &read);
        assert_eq!(score, 1, "{:?}", ops);
        let read = into_reads(vec![0, 1]);
        let query = into_reads(vec![0, 1, 2, 3, 4]);
        let (score, ops) = pairwise_alignment_gotoh(&read, &query);
        assert_eq!(score, 2, "{:?}", ops);
        let (score, ops) = pairwise_alignment_gotoh(&query, &read);
        assert_eq!(score, 2, "{:?}", ops);
    }
    // #[test]
    // fn aln_test_gotoh_2() {
    //     let into_read = |nodes: Vec<(u64, u64, bool)>| {
    //         let nodes: Vec<_> = nodes
    //             .into_iter()
    //             .map(|(unit, cluster, is_forward)| LightNode {
    //                 prev_offset: None,
    //                 unit,
    //                 cluster,
    //                 is_forward,
    //                 after_offset: None,
    //             })
    //             .collect();
    //         ReadSkelton { nodes }
    //     };
    //     let read = vec![(0, 0, true), (2, 0, true), (3, 0, true), (4, 0, true)];
    //     let read = into_read(read);
    //     let query = vec![(1, 0, true), (2, 0, true), (3, 0, true)];
    //     let query = into_read(query);
    //     let (score, ops) = pairwise_alignment_gotoh(&read, &query);
    //     assert_eq!(score, 2);
    //     use definitions::Op;
    //     let answer = vec![Op::Ins(1), Op::Del(1), Op::Match(2), Op::Del(1)];
    //     assert_eq!(ops, answer);
    // }
    #[test]
    fn max_chain_test() {
        let penalty = 3;
        let input = vec![vec![], vec![(0, 1)], vec![(0, 5), (1, 1)]];
        let removed = max_chain(3, &input, penalty);
        assert_eq!(removed, vec![1]);
        let input = vec![vec![], vec![(0, 1)], vec![(0, 2), (1, 1)]];
        let removed = max_chain(3, &input, penalty);
        assert_eq!(removed, vec![]);
        let input = vec![vec![], vec![(0, 1)], vec![(1, 1)]];
        let removed = max_chain(3, &input, penalty);
        assert_eq!(removed, vec![]);
        let input = vec![
            vec![],                          //0
            vec![(0, 1)],                    // 1
            vec![(0, 20), (1, 1)],           // 2
            vec![(2, 0)],                    // 3
            vec![(0, 100), (2, 90), (3, 0)], // 4
        ];
        let removed = max_chain(5, &input, penalty);
        assert_eq!(removed, vec![3, 1]);
    }
}
