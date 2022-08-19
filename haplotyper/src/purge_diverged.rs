//! This modeule includes methods to newly define units
//! From the result of local clustering.
use definitions::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Default)]
pub struct PurgeDivConfig {}

impl PurgeDivConfig {
    pub fn new() -> Self {
        Self {}
    }
}

#[derive(Debug, Clone, Default)]
pub struct PurgeLargeDelConfig {
    indel_size: usize,
    occupy_fraction: f64,
}

impl PurgeLargeDelConfig {
    pub fn new(indel_size: usize, occupy_fraction: f64) -> Self {
        Self {
            indel_size,
            occupy_fraction,
        }
    }
}

pub trait PurgeDivergent {
    fn purge(&mut self, config: &PurgeDivConfig);
    // Return the ids of the units.
    fn purge_largeindel(&mut self, config: &PurgeLargeDelConfig) -> HashSet<u64>;
}

use rayon::prelude::*;
const THR: f64 = 0.1;
// const SD_SIGMA: f64 = 8f64;
impl PurgeDivergent for DataSet {
    fn purge(&mut self, config: &PurgeDivConfig) {
        let prev = self.encoded_reads.len();
        let mut purged_cluster = HashSet::new();
        purged_cluster.extend(purge_diverged_nodes(self, THR, config));
        debug!("PD\tEncodedRead\t{}\t{}\t0", prev, self.encoded_reads.len());
        re_cluster(self, &purged_cluster);
    }
    fn purge_largeindel(&mut self, config: &PurgeLargeDelConfig) -> HashSet<u64> {
        let prev = self.encoded_reads.len();
        let purged_cluster = purge_large_deletion_nodes(self, config);
        self.encoded_reads.retain(|r| !r.nodes.is_empty());
        debug!(
            "PLI\tEncodedRead\t{}\t{}\t1",
            prev,
            self.encoded_reads.len()
        );
        purged_cluster
    }
}

const ACCEPT_RATE: f64 = 0.5;
const DEL_WEIGHT: i64 = 2;
const INS_WEIGHT: i64 = 0;
const MATCH_WEIGHT: i64 = 1;
fn purge_large_deletion_nodes(ds: &mut DataSet, config: &PurgeLargeDelConfig) -> HashSet<u64> {
    let mut indel_size_distr: HashMap<(u64, u64), Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for (idx, node) in read.nodes.iter().enumerate() {
            let xs = node.cigar.iter().map(|&op| match op {
                definitions::Op::Match(l) => -(l as i64) * MATCH_WEIGHT,
                definitions::Op::Del(l) => l as i64 * DEL_WEIGHT,
                definitions::Op::Ins(l) => l as i64 * INS_WEIGHT,
            });
            let del_size = crate::misc::max_region(xs) / DEL_WEIGHT;
            indel_size_distr
                .entry((node.unit, node.cluster))
                .or_default()
                .push((read.id, idx, del_size));
        }
    }
    let indel_size = config.indel_size;
    crate::misc::update_coverage(ds);
    let accept_size = (indel_size as f64 * ACCEPT_RATE).ceil() as usize;
    indel_size_distr.retain(|(unit, cluster), distr| {
        let dip_coverage = 2f64 * ds.coverage.unwrap();
        let lower_thr = (dip_coverage / 2f64 * config.occupy_fraction).floor() as usize;
        let upper_thr = (distr.len() as f64 - lower_thr as f64).max(0f64).ceil() as usize;
        let num = distr
            .iter()
            .map(|x| x.2)
            .filter(|&x| indel_size < x as usize)
            .count();
        if lower_thr < num {
            let sum: i64 = distr
                .iter()
                .map(|x| x.2)
                .filter(|&x| (accept_size as i64) < x)
                .sum();
            let mean = sum / num as i64;
            let cov = distr.len();
            let retain = (lower_thr..upper_thr).contains(&num);
            debug!("PLI\t{unit}\t{cluster}\t{mean}\t{num}\t{cov}\t{retain}");
        }
        if (lower_thr..upper_thr).contains(&num) {
            distr.retain(|&(_, _, size)| accept_size < size as usize);
            true
        } else {
            false
        }
    });
    let mut read_bucket: HashMap<_, Vec<_>> = HashMap::new();
    for (&(unit, cluster), distr) in indel_size_distr.iter() {
        for &(readid, idx, _) in distr.iter() {
            read_bucket
                .entry(readid)
                .or_default()
                .push((unit, cluster, idx));
        }
    }
    read_bucket
        .values_mut()
        .for_each(|bucket| bucket.sort_by_key(|x| std::cmp::Reverse(x.2)));
    for read in ds.encoded_reads.iter_mut() {
        let bucket = match read_bucket.get(&read.id) {
            Some(bucket) => bucket,
            None => continue,
        };
        for &(unit, cluster, idx) in bucket {
            let node = &read.nodes[idx];
            assert_eq!((unit, cluster), (node.unit, node.cluster));
            read.remove(idx);
        }
    }
    indel_size_distr.keys().map(|x| x.0).collect()
}

// Purge node with very high error rate.
pub fn purge_erroneous_nodes(ds: &mut DataSet, _config: &PurgeDivConfig) -> HashSet<u64> {
    const SAFE_MARGIN: f64 = 5f64;
    let fallback = crate::determine_units::calc_sim_thr(ds, 0.5);
    use crate::estimate_error_rate::estimate_error_rate;
    let errors = estimate_error_rate(ds, fallback);
    let sigma_of_er = errors.median_of_sqrt_err;
    debug!("PD\tPurgeErroneousNode\t{SAFE_MARGIN}\t{sigma_of_er}",);
    let mut units_of_removed_nodes = HashSet::new();
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for read in ds.encoded_reads.iter_mut() {
        let read_error = errors.read(read.id);
        read.nodes.retain(|node| {
            let unit_error = errors.unit((node.unit, node.cluster));
            let aln = node.recover(units[&node.unit]).1;
            let errors = aln.iter().filter(|&&x| x != b'|').count();
            let error = errors as f64 / aln.len() as f64;
            let to_retain = error < read_error + unit_error + SAFE_MARGIN * sigma_of_er;
            if !to_retain {
                units_of_removed_nodes.insert(node.unit);
            }
            to_retain
        });
    }
    ds.encoded_reads.retain(|r| !r.nodes.is_empty());
    units_of_removed_nodes
}

fn re_cluster(ds: &mut DataSet, selection: &HashSet<u64>) {
    let copy_number: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.copy_num))
        .collect();

    {
        // Preserve clustering.
        let preserve: Vec<_> = ds
            .encoded_reads
            .iter()
            .map(|r| {
                let cls: Vec<_> = r.nodes.iter().map(|n| n.cluster).collect();
                (r.id, cls, r.nodes.len())
            })
            .collect();
        ds.encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
            .for_each(|n| n.cluster = 0);
        use crate::multiplicity_estimation::*;
        let multip_config = MultiplicityEstimationConfig::new(230493, None);
        ds.estimate_multiplicity(&multip_config);
        for (r, (id, cls, len)) in ds.encoded_reads.iter_mut().zip(preserve) {
            assert_eq!(r.id, id);
            assert_eq!(r.nodes.len(), len);
            r.nodes
                .iter_mut()
                .zip(cls)
                .for_each(|(n, cl)| n.cluster = cl);
        }
    }
    let changed_units: HashSet<_> = ds
        .selected_chunks
        .iter()
        .filter_map(|c| {
            copy_number
                .get(&c.id)
                .filter(|&&copy_num| copy_num != c.copy_num)
                .map(|_| c.id)
        })
        .chain(selection.iter().copied())
        .collect();
    debug!("PD\tReClustering\t{}", changed_units.len());
    crate::local_clustering::local_clustering_selected(ds, &changed_units);
}

fn purge_diverged_nodes(ds: &mut DataSet, thr: f64, config: &PurgeDivConfig) -> HashSet<u64> {
    let diverged_clusters = get_diverged_clusters_dev(ds, thr, config);
    // If the all the cluster is labelled as "diverged", it is the fault of the consensus...,
    let diverged_clusters: Vec<_> = diverged_clusters
        .into_iter()
        .map(|mut xs| {
            if xs.iter().all(|&x| x) {
                xs.iter_mut().for_each(|x| *x = false);
            }
            xs
        })
        .collect();
    let removed_clusters: usize = diverged_clusters
        .iter()
        .map(|cls| cls.iter().filter(|&&x| x).count())
        .sum();
    debug!("PD\tPurge\t{}", removed_clusters);
    for (id, x) in diverged_clusters.iter().enumerate() {
        for (cl, _) in x.iter().enumerate().filter(|x| *x.1) {
            debug!("PD\tPurge\t{}\t{}", id, cl);
        }
    }
    // Modify cluster number.
    for unit in ds.selected_chunks.iter_mut() {
        let id = unit.id as usize;
        let squished = diverged_clusters[id].iter().filter(|&&x| x).count();
        unit.cluster_num -= squished;
    }
    let removed_nodes: usize = ds
        .encoded_reads
        .iter_mut()
        .map(|read| {
            let orig = read.nodes.len();
            read.nodes
                .retain(|node| !diverged_clusters[node.unit as usize][node.cluster as usize]);
            for node in read.nodes.iter_mut() {
                let cluster_info = &diverged_clusters[node.unit as usize];
                if cluster_info.iter().any(|&b| b) {
                    remove_diverged(node, cluster_info);
                }
            }
            orig - read.nodes.len()
        })
        .sum();
    debug!("PD\tREMOVED\t{}\t{THR}", removed_nodes);
    ds.encoded_reads.retain(|read| !read.nodes.is_empty());
    let seqs: HashMap<_, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        let seq = seqs[&read.id];
        *read = crate::encode::nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    });
    diverged_clusters
        .iter()
        .enumerate()
        .filter_map(|(id, is_div)| is_div.iter().any(|&x| x).then(|| id as u64))
        .collect()
}

// Unit -> Cluster -> If the cluster is very diverged.
pub fn get_diverged_clusters_dev(ds: &DataSet, thr: f64, _: &PurgeDivConfig) -> Vec<Vec<bool>> {
    use crate::estimate_error_rate::estimate_error_rate;
    // let error_rate = ds.error_rate();
    // let fallback = error_rate.total + SD_SIGMA * error_rate.total_sd;
    let fallback = crate::determine_units::calc_sim_thr(ds, 0.5);
    debug!("PD\tFallback\t{fallback}");
    let error_rates = estimate_error_rate(ds, fallback);
    error_rates
        .unit_error_rate
        .iter()
        .map(|es| es.iter().map(|&x| thr < x).collect())
        .collect()
}

// Unit -> Cluster -> If the cluster is very diverged.
// fn get_diverged_clusters(ds: &DataSet, thr: f64, _config: &PurgeDivConfig) -> Vec<Vec<bool>> {
//     let max_unit_id = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
//     let mut slots = vec![vec![]; max_unit_id as usize + 1];
//     for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
//         slots[node.unit as usize].push(node);
//     }
//     let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
//     slots
//         .iter()
//         .enumerate()
//         .map(|(id, nodes)| {
//             assert!(nodes.iter().all(|n| n.unit as usize == id));
//             let ref_unit = match chunks.get(&(id as u64)) {
//                 Some(res) => res,
//                 None => return Vec::new(),
//             };
//             let mut div_rates = vec![vec![]; ref_unit.cluster_num];
//             for node in nodes.iter() {
//                 let (_, aln, _) = node.recover(ref_unit);
//                 let div = aln.iter().filter(|&&x| x != b'|').count() as f64 / aln.len() as f64;
//                 div_rates[node.cluster as usize].push(div);
//             }
//             div_rates
//                 .iter()
//                 .map(|rates| {
//                     let sum: f64 = rates.iter().sum();
//                     thr < sum / rates.len() as f64
//                 })
//                 .collect()
//         })
//         .collect()
// }

fn remove_diverged(node: &mut Node, cluster_info: &[bool]) {
    node.cluster -= cluster_info
        .iter()
        .take(node.cluster as usize)
        .filter(|&&x| x)
        .count() as u64;
    let mut i = 0;
    node.posterior.retain(|_| {
        i += 1;
        !cluster_info[i - 1]
    });
}

// ///Return mapping from one chunk to another.
// ///The `selected_chunks` member would be modified.
// fn purge_diverged_nodes_dev(
//     ds: &mut DataSet,
//     thr: f64,
//     config: &PurgeDivConfig,
// ) -> (HashMap<(u64, u64), (u64, u64)>, HashSet<u64>) {
//     let diverged_clusters = get_diverged_clusters(ds, thr, config);
//     // If the all the cluster is labelled as "diverged", it is the fault of the consensus...,
//     let diverged_clusters: Vec<Vec<bool>> = diverged_clusters
//         .into_iter()
//         .map(|mut xs| {
//             if xs.iter().all(|&x| x) {
//                 xs.iter_mut().for_each(|x| *x = false);
//             }
//             xs
//         })
//         .collect();
//     for (uid, cls) in diverged_clusters.iter().enumerate() {
//         for (cl, _) in cls.iter().enumerate().filter(|&(_, &b)| b) {
//             debug!("PD\t{uid}\t{cl}");
//         }
//     }
//     let has_diverged_cluster: Vec<bool> = diverged_clusters
//         .iter()
//         .map(|xs| xs.iter().any(|&x| x))
//         .collect();
//     let removed_clusters: usize = diverged_clusters
//         .iter()
//         .map(|cls| cls.iter().filter(|&&x| x).count())
//         .sum();
//     debug!("PD\tPurge\t{}", removed_clusters);
//     // Modify cluster number.
//     let mut new_units = Vec::new();
//     let max_id = ds.selected_chunks.iter().map(|u| u.id).max().unwrap();
//     let mut changed = HashSet::new();
//     let mappings = ds
//         .selected_chunks
//         .iter_mut()
//         .filter(|u| has_diverged_cluster[u.id as usize])
//         .flat_map(|u| {
//             changed.insert(u.id);
//             // Tune cluster number.
//             let squished = diverged_clusters[u.id as usize]
//                 .iter()
//                 .filter(|&&x| x)
//                 .count();
//             u.cluster_num -= squished;
//             u.copy_num -= squished;
//             // Create new chunk.
//             let id = max_id + 1 + new_units.len() as u64;
//             changed.insert(id);
//             let mut new_unit = Unit::new(id, u.seq.clone().into(), squished);
//             new_unit.cluster_num = squished;
//             new_units.push(new_unit);
//             // Determine mappings.
//             let mut mapping = Vec::with_capacity(diverged_clusters[u.id as usize].len());
//             let (mut new_unit_cl, mut old_unit_cl) = (0, 0);
//             for (from, is_diverged) in diverged_clusters[u.id as usize].iter().enumerate() {
//                 let from = (u.id, from as u64);
//                 match is_diverged {
//                     true => {
//                         mapping.push((from, (id, new_unit_cl)));
//                         new_unit_cl += 1;
//                     }
//                     false => {
//                         mapping.push((from, (u.id, old_unit_cl)));
//                         old_unit_cl += 1;
//                     }
//                 }
//             }
//             mapping
//         })
//         .collect();
//     ds.selected_chunks.extend(new_units);
//     (mappings, changed)
// }
