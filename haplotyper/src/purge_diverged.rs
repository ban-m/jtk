//! This modeule includes methods to newly define chunks
//! From the result of local clustering.
use definitions::*;
use log::*;

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
    compress: bool,
}

impl PurgeLargeDelConfig {
    pub fn new(indel_size: usize, occupy_fraction: f64, compress: bool) -> Self {
        Self {
            indel_size,
            occupy_fraction,
            compress,
        }
    }
}

pub trait PurgeDivergent {
    fn purge(&mut self, config: &PurgeDivConfig);
    // Return the ids of the chunks.
    fn purge_largeindel(&mut self, config: &PurgeLargeDelConfig) -> HashSet<u64>;
}

use rayon::prelude::*;
const THR: f64 = 0.1;
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
        debug!("PLI\tUnits\t{:?}", purged_cluster);
        purged_cluster
    }
}

const ACCEPT_RATE: f64 = 0.5;
const DEL_WEIGHT: i64 = 2;
const MATCH_WEIGHT: i64 = 1;
fn purge_large_deletion_nodes(ds: &mut DataSet, config: &PurgeLargeDelConfig) -> HashSet<u64> {
    let mut indel_size_distr: HashMap<(u64, u64), Vec<_>> = HashMap::new();
    let copy_num: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.copy_num))
        .collect();
    for read in ds.encoded_reads.iter() {
        for (idx, node) in read.nodes.iter().enumerate() {
            let xs = node.cigar.iter().map(|&op| match op {
                definitions::Op::Match(l) => -(l as i64) * MATCH_WEIGHT,
                definitions::Op::Del(l) => l as i64 * DEL_WEIGHT,
                definitions::Op::Ins(l) => l as i64 * DEL_WEIGHT,
            });
            let del_size = crate::misc::max_region(xs) / DEL_WEIGHT;
            let key = match config.compress {
                true => (node.chunk, 0),
                false => (node.chunk, node.cluster),
            };
            let pointer = (read.id, idx, del_size);
            indel_size_distr.entry(key).or_default().push(pointer);
        }
    }
    // TODO: If there are two more clusters, and all the read in one cluster
    // support large indel, we can .. ?
    crate::misc::update_coverage(ds);
    indel_size_distr.retain(|&(chunk, cluster), distr| {
        let copy_num = copy_num[&chunk];
        let hap_coverage = ds.coverage.unwrap();
        match copy_num {
            0 | 1 => filter_single_copy(distr, chunk, cluster, hap_coverage, config),
            _ => filter_multi_copy(distr, chunk, cluster, hap_coverage, config),
        }
    });
    let mut read_bucket: HashMap<_, Vec<_>> = HashMap::new();
    for (&(chunk, cluster), distr) in indel_size_distr.iter() {
        for &(readid, idx, _) in distr.iter() {
            read_bucket
                .entry(readid)
                .or_default()
                .push((chunk, cluster, idx));
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
        for &(chunk, cluster, idx) in bucket {
            let node = &read.nodes[idx];
            if config.compress {
                assert_eq!(chunk, node.chunk);
            } else {
                assert_eq!((chunk, cluster), (node.chunk, node.cluster));
            }
            read.remove(idx);
        }
    }
    indel_size_distr.keys().map(|x| x.0).collect()
}

fn filter_single_copy(
    distr: &mut Vec<(u64, usize, i64)>,
    chunk: u64,
    cluster: u64,
    _hap_coverage: f64,
    config: &PurgeLargeDelConfig,
) -> bool {
    let indel_size = config.indel_size;
    let accept_size = (indel_size as f64 * ACCEPT_RATE).ceil() as usize;
    let num = distr.iter().filter(|&x| (indel_size as i64) < x.2).count();
    let upper_thr = ((distr.len() as f64) * (1f64 - config.occupy_fraction / 2f64)).ceil() as usize;
    if (1..upper_thr).contains(&num) {
        let sum: i64 = distr
            .iter()
            .map(|x| x.2)
            .filter(|&x| (accept_size as i64) < x)
            .sum();
        let mean = sum / (num as i64).max(1);
        let cov = distr.len();
        debug!("PLI\t{chunk}\t{cluster}\t{mean}\t{num}\t{cov}\tSingle");
        distr.retain(|&(_, _, size)| accept_size < size as usize);
        true
    } else {
        false
    }
}

fn filter_multi_copy(
    distr: &mut Vec<(u64, usize, i64)>,
    chunk: u64,
    cluster: u64,
    hap_coverage: f64,
    config: &PurgeLargeDelConfig,
) -> bool {
    let indel_size = config.indel_size;
    let accept_size = (indel_size as f64 * ACCEPT_RATE).ceil() as usize;
    let lower_thr = (hap_coverage * config.occupy_fraction).floor() as usize;
    let upper_thr = (distr.len() as f64 - lower_thr as f64).max(0f64).ceil() as usize;
    let min_remove = distr.iter().filter(|&x| (indel_size as i64) < x.2).count();
    let max_remove = distr.iter().filter(|&x| (accept_size as i64) < x.2).count();
    if chunk == 1527 {
        debug!("PLI\t{chunk}\t{cluster}\t{min_remove}\t{max_remove}");
    }
    if lower_thr < min_remove && max_remove < upper_thr {
        let sum: i64 = distr
            .iter()
            .map(|x| x.2)
            .filter(|&x| (accept_size as i64) < x)
            .sum();
        let mean = sum / max_remove as i64;
        let cov = distr.len();
        debug!("PLI\t{chunk}\t{cluster}\t{mean}\t{min_remove}\t{max_remove}\t{cov}");
        distr.retain(|&(_, _, size)| accept_size < size as usize);
        true
    } else {
        false
    }
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
    let changed_chunks: HashSet<_> = ds
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
    debug!("PD\tReClustering\t{}", changed_chunks.len());
    use crate::local_clustering::LocalClustering;
    ds.local_clustering_selected(&changed_chunks);
}

fn purge_diverged_nodes(ds: &mut DataSet, thr: f64, config: &PurgeDivConfig) -> HashSet<u64> {
    let diverged_clusters = get_diverged_clusters(ds, thr, config);
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
    for chunk in ds.selected_chunks.iter_mut() {
        let id = chunk.id as usize;
        let squished = diverged_clusters[id].iter().filter(|&&x| x).count();
        chunk.cluster_num -= squished;
    }
    let removed_nodes: usize = ds
        .encoded_reads
        .iter_mut()
        .map(|read| {
            let orig = read.nodes.len();
            read.nodes
                .retain(|node| !diverged_clusters[node.chunk as usize][node.cluster as usize]);
            for node in read.nodes.iter_mut() {
                let cluster_info = &diverged_clusters[node.chunk as usize];
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
        .filter_map(|(id, is_div)| is_div.iter().any(|&x| x).then_some(id as u64))
        .collect()
}

// Unit -> Cluster -> If the cluster is very diverged.
fn get_diverged_clusters(ds: &DataSet, thr: f64, _: &PurgeDivConfig) -> Vec<Vec<bool>> {
    use crate::estimate_error_rate::estimate_error_rate;
    let fallback = crate::determine_chunks::calc_sim_thr(ds, 0.5);
    debug!("PD\tFallback\t{fallback}");
    let error_rates = estimate_error_rate(ds, fallback);
    error_rates
        .chunk_error_rate
        .iter()
        .map(|es| es.iter().map(|&x| thr < x).collect())
        .collect()
}

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
