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
    fn purge_dev(&mut self, config: &PurgeDivConfig);
}

use rayon::prelude::*;

const SD_SIGMA: f64 = 6f64;
const REDEFINE: bool = false;
use crate::stats::Stats;
impl PurgeDivergent for DataSet {
    fn purge(&mut self, config: &PurgeDivConfig) {
        let error_rate = self.error_rate();
        let thr = error_rate.total + SD_SIGMA * error_rate.total_sd;
        debug!(
            "PD\tTHR\t{}\t{}\t{}",
            error_rate.total, error_rate.total_sd, thr
        );
        let prev = self.encoded_reads.len();
        if REDEFINE {
            let re_defined_chunks = purge_diverged_nodes_dev(self, thr, config);
            let changed_nodes: usize = self
                .encoded_reads
                .par_iter_mut()
                .map(|read| {
                    read.nodes
                        .iter_mut()
                        .filter(|n| re_defined_chunks.contains_key(&(n.unit, n.cluster)))
                        .map(|n| {
                            let (u, c) = re_defined_chunks[&(n.unit, n.cluster)];
                            n.unit = u;
                            n.cluster = c;
                        })
                        .count()
                })
                .sum();
            debug!("PD\tChanged\t{}", changed_nodes);
        } else {
            let to_be_removed = purge_diverged_nodes(self, thr, config);
            let removed_nodes: usize = to_be_removed.iter().map(|(_, xs)| xs.len()).sum();
            debug!("PD\tREMOVED\t{}", removed_nodes);
            let seqs: HashMap<_, _> = self.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
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
        }
        debug!("PD\tEncodedRead\t{}\t{}", prev, self.encoded_reads.len());
        // Reclustering...
        re_cluster(self, config.threads);
    }
    fn purge_dev(&mut self, config: &PurgeDivConfig) {
        purge_disjoint_cluster(self, config);
        re_cluster(self, config.threads);
    }
}

fn purge_disjoint_cluster(ds: &mut DataSet, _config: &PurgeDivConfig) {
    let mut pileups: HashMap<_, Vec<_>> = HashMap::new();
    ds.encoded_reads
        .iter_mut()
        .flat_map(|r| r.nodes.iter_mut())
        .for_each(|n| {
            pileups.entry(n.unit).or_default().push(n);
        });
    let mut new_units = Vec::new();
    let unit_seqs: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c.seq())).collect();
    let mut max_unit_id = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    for (id, pileup) in pileups.into_iter().filter(|(_, p)| !p.is_empty()) {
        let seq = unit_seqs.get(&id).unwrap();
        let split_units = split_disjoint_cluster(pileup, seq, max_unit_id);
        max_unit_id += split_units.len() as u64;
        new_units.extend(split_units);
    }
}

const LK_THRESHOLD: f64 = -20f64;
fn split_disjoint_cluster(mut nodes: Vec<&mut Node>, seq: &[u8], unit_id: u64) -> Vec<Unit> {
    let num_cluster = nodes[0].posterior.len();
    assert!(nodes.iter().all(|n| n.posterior.len() == num_cluster));
    let mut fu = crate::find_union::FindUnion::new(num_cluster);
    for posterior in nodes.iter().map(|n| n.posterior.as_slice()) {
        let mut high_post = posterior
            .iter()
            .enumerate()
            .filter(|&x| *x.1 > LK_THRESHOLD);
        let mut current = match high_post.next() {
            Some(idx) => idx.0,
            None => continue,
        };
        for (next, _) in high_post {
            fu.unite(next, current);
            current = next;
        }
    }
    let (cluster_mapped_to, cluster_sizes) = {
        let mut map = HashMap::new();
        let mut cluster_size = Vec::new();
        for i in 0..num_cluster {
            let len = map.len();
            let parent = fu.find(i).unwrap();
            if !map.contains_key(&parent) {
                map.insert(i, len);
                cluster_size.push(0);
            }
            cluster_size[i] += 1;
        }
        let cluster_mapped_to: Vec<_> = (0..num_cluster)
            .map(|i| map[&fu.find(i).unwrap()])
            .collect();
        (cluster_mapped_to, cluster_size)
    };
    // There's no split. Retain clusterings.
    if cluster_sizes.len() == 1 {
        Vec::new()
    } else {
        // There ARE splits. Erase clustering information.
        for node in nodes.iter_mut() {
            let mapped_to = cluster_mapped_to[node.cluster as usize];
            if mapped_to != 0 {
                node.unit = unit_id + mapped_to as u64;
            }
            let cluster_size = cluster_sizes[mapped_to];
            node.posterior.clear();
            node.posterior
                .extend(std::iter::repeat((cluster_size as f64).recip()).take(cluster_size));
            node.cluster = 0;
        }
        let seq = String::from_utf8_lossy(seq).to_string();
        cluster_sizes
            .into_iter()
            .enumerate()
            .skip(1)
            .map(|(i, cl)| Unit::new(unit_id + i as u64, seq.clone(), cl))
            .collect()
    }
}

fn re_cluster(ds: &mut DataSet, threads: usize) {
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
        let multip_config = MultiplicityEstimationConfig::new(threads, 230493, None);
        ds.estimate_multiplicity(&multip_config);
        // Recover.
        for (r, (id, cls, len)) in ds.encoded_reads.iter_mut().zip(preserve) {
            assert_eq!(r.id, id);
            assert_eq!(r.nodes.len(), len);
            r.nodes
                .iter_mut()
                .zip(cls)
                .for_each(|(n, cl)| n.cluster = cl);
        }
    }
    use std::collections::HashSet;
    let changed_units: HashSet<_> = ds
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
    crate::local_clustering::local_clustering_selected(ds, &changed_units);
}

fn purge_diverged_nodes(
    ds: &mut DataSet,
    thr: f64,
    config: &PurgeDivConfig,
) -> Vec<(u64, Vec<usize>)> {
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

fn purge_diverged_nodes_dev(
    _ds: &mut DataSet,
    _thr: f64,
    _config: &PurgeDivConfig,
) -> HashMap<(u64, u64), (u64, u64)> {
    todo!()
    // let diverged_clusters = get_diverged_clusters(ds, thr, config);
    // // If the all the cluster is labelled as "diverged", it is the fault of the consensus...,
    // let diverged_clusters: Vec<_> = diverged_clusters
    //     .into_iter()
    //     .map(|mut xs| {
    //         if xs.iter().all(|&x| x) {
    //             xs.iter_mut().for_each(|x| *x = false);
    //         }
    //         xs
    //     })
    //     .collect();
    // let has_diverged_cluster: Vec<bool> = diverged_clusters
    //     .iter()
    //     .map(|xs| xs.iter().any(|&x| x))
    //     .collect();
    // let removed_clusters: usize = diverged_clusters
    //     .iter()
    //     .map(|cls| cls.iter().filter(|&&x| x).count())
    //     .sum();
    // debug!("PD\tPurge\t{}", removed_clusters);
    // for (id, x) in diverged_clusters.iter().enumerate() {
    //     for (cl, _) in x.iter().enumerate().filter(|x| *x.1) {
    //         debug!("PD\tPurge\t{}\t{}", id, cl);
    //     }
    // }
    // // Modify cluster number.
    // for unit in ds
    //     .selected_chunks
    //     .iter_mut()
    //     .filter(|u| has_diverged_cluster[u.id as usize])
    // {
    //     let squished = diverged_clusters[unit.id as usize]
    //         .iter()
    //         .filter(|&&x| x)
    //         .count();
    //     unit.cluster_num -= squished;
    // }
    // // Maybe we need to modify the copy number, thoguth....
    // ds.encoded_reads
    //     .iter_mut()
    //     .map(|r| {
    //         let indices: Vec<_> = r
    //             .nodes
    //             .iter_mut()
    //             .enumerate()
    //             .filter(|(_, node)| has_diverged_cluster[node.unit as usize])
    //             .filter_map(|(idx, node)| {
    //                 let cluster_info = &diverged_clusters[node.unit as usize];
    //                 if cluster_info[node.cluster as usize] {
    //                     Some(idx)
    //                 } else {
    //                     remove_diverged(node, cluster_info);
    //                     None
    //                 }
    //             })
    //             .collect();
    //         (r.id, indices)
    //     })
    //     .collect()
}
