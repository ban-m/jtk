//! This modeule includes methods to newly define units
//! From the result of local clustering.
use definitions::*;
use std::collections::HashMap;
#[derive(Debug, Clone, Default)]
pub struct PurgeDivConfig {
    threads: usize,
}

impl PurgeDivConfig {
    pub fn new(threads: usize) -> Self {
        Self { threads }
    }
}

pub trait PurgeDivergent {
    fn purge(&mut self, config: &PurgeDivConfig);
}

use rayon::prelude::*;
impl PurgeDivergent for DataSet {
    fn purge(&mut self, config: &PurgeDivConfig) {
        let copy_number: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.copy_num))
            .collect();
        let error_rate = crate::light_stats::error_rate(self);
        let thr = error_rate.total + 4f64 * error_rate.total_sd;
        debug!(
            "PD\tTHR\t{}\t{}\t{}",
            error_rate.total, error_rate.total_sd, thr
        );
        let to_be_removed = purge_diverged_nodes(self, thr, config);
        let removed_nodes: usize = to_be_removed.iter().map(|(_, xs)| xs.len()).sum();
        debug!("PD\tREMOVED\t{}", removed_nodes);
        let seqs: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
        let prev = self.encoded_reads.len();
        self.encoded_reads
            .par_iter_mut()
            .zip(to_be_removed)
            .for_each(|(read, (id, indices))| {
                assert_eq!(read.id, id);
                assert!(indices.is_sorted());
                for &index in indices.iter().rev() {
                    read.remove(index);
                }
            });
        self.encoded_reads.retain(|read| !read.nodes.is_empty());
        self.encoded_reads.par_iter_mut().for_each(|read| {
            let mut nodes = vec![];
            nodes.append(&mut read.nodes);
            let seq = seqs[&read.id];
            *read = crate::encode::nodes_to_encoded_read(read.id, nodes, seq).unwrap();
        });
        debug!("PD\tEncodedRead\t{}\t{}", prev, self.encoded_reads.len());
        // Reclustering...
        {
            // Preserve clustering.
            let preserve: Vec<_> = self
                .encoded_reads
                .iter()
                .map(|r| {
                    let cls: Vec<_> = r.nodes.iter().map(|n| n.cluster).collect();
                    (r.id, cls, r.nodes.len())
                })
                .collect();
            self.encoded_reads
                .iter_mut()
                .flat_map(|r| r.nodes.iter_mut())
                .for_each(|n| n.cluster = 0);
            use crate::multiplicity_estimation::*;
            let multip_config = MultiplicityEstimationConfig::new(config.threads, 230493, None);
            self.estimate_multiplicity(&multip_config);
            // Recover.
            for (r, (id, cls, len)) in self.encoded_reads.iter_mut().zip(preserve) {
                assert_eq!(r.id, id);
                assert_eq!(r.nodes.len(), len);
                r.nodes
                    .iter_mut()
                    .zip(cls)
                    .for_each(|(n, cl)| n.cluster = cl);
            }
        }
        use std::collections::HashSet;
        let changed_units: HashSet<_> = self
            .selected_chunks
            .iter()
            .filter_map(|c| {
                copy_number
                    .get(&c.id)
                    .filter(|&&copy_num| copy_num != c.copy_num)
                    .map(|_| c.id)
            })
            .collect();
        debug!("PD\tReClustering\t{}", changed_units.len());
        crate::local_clustering::local_clustering_selected(self, &changed_units);
    }
}

fn purge_diverged_nodes(
    ds: &mut DataSet,
    thr: f64,
    config: &PurgeDivConfig,
) -> Vec<(u64, Vec<usize>)> {
    let diverged_clusters = get_diverged_clusters(ds, thr, config);
    let has_diverged_cluster: Vec<bool> = diverged_clusters
        .iter()
        .map(|xs| xs.iter().any(|&x| x))
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
    for unit in ds
        .selected_chunks
        .iter_mut()
        .filter(|u| has_diverged_cluster[u.id as usize])
    {
        let squished = diverged_clusters[unit.id as usize]
            .iter()
            .filter(|&&x| x)
            .count();
        unit.cluster_num -= squished;
    }
    // Maybe we need to modify the copy number, thoguth....
    ds.encoded_reads
        .iter_mut()
        .map(|r| {
            let indices: Vec<_> = r
                .nodes
                .iter_mut()
                .enumerate()
                .filter(|(_, node)| has_diverged_cluster[node.unit as usize])
                .filter_map(|(idx, node)| {
                    let cluster_info = &diverged_clusters[node.unit as usize];
                    if cluster_info[node.cluster as usize] {
                        Some(idx)
                    } else {
                        remove_diverged(node, cluster_info);
                        None
                    }
                })
                .collect();
            (r.id, indices)
        })
        .collect()
}

// Unit -> Cluster -> If the cluster is very diverged.
fn get_diverged_clusters(ds: &DataSet, thr: f64, _config: &PurgeDivConfig) -> Vec<Vec<bool>> {
    let max_unit_id = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    let mut slots = vec![vec![]; max_unit_id as usize + 1];
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        slots[node.unit as usize].push(node);
    }
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    slots
        .iter()
        .enumerate()
        .map(|(id, nodes)| {
            assert!(nodes.iter().all(|n| n.unit as usize == id));
            let ref_unit = chunks.get(&(id as u64)).unwrap();
            let mut div_rates = vec![vec![]; ref_unit.cluster_num];
            for node in nodes.iter() {
                let (_, aln, _) = node.recover(ref_unit);
                let div = aln.iter().filter(|&&x| x != b'|').count() as f64 / aln.len() as f64;
                div_rates[node.cluster as usize].push(div);
            }
            div_rates
                .iter()
                .map(|rates| {
                    let sum: f64 = rates.iter().sum();
                    thr < sum / rates.len() as f64
                })
                .collect()
        })
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
