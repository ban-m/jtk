use definitions::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone)]
pub struct ReClusteringConfig {
    threads: usize,
    repeat_num: usize,
    coverage_thr: usize,
}

impl ReClusteringConfig {
    pub fn new(threads: usize, repeat_num: usize, coverage_thr: usize) -> Self {
        Self {
            threads,
            repeat_num,
            coverage_thr,
        }
    }
}

pub trait ReClustering {
    fn re_clustering(self, c: &ReClusteringConfig) -> Self;
}

impl ReClustering for DataSet {
    /// Assemble the dataset, re-clustering units with copy number more than 2.
    fn re_clustering(mut self, c: &ReClusteringConfig) -> DataSet {
        // Remember the initial read clustering.
        let init_clustering: Vec<_> = self
            .encoded_reads
            .iter()
            .map(|r| {
                let id = r.id;
                let cls: Vec<_> = r.nodes.iter().map(|n| (n.unit, n.cluster)).collect();
                (id, cls)
            })
            .collect();
        let cluster_num: HashMap<_, _> = self
            .selected_chunks
            .iter()
            .map(|c| (c.id, c.cluster_num))
            .collect();
        // Polish clusering.
        use crate::em_correction::ClusteringCorrection;
        self = self.correct_clustering_em(c.repeat_num, c.coverage_thr, true);
        // ID of the cluster -> total excess of the clustering.
        let re_cluster = get_units_to_cluster(&self, c);
        // Recover the initial clustering.
        for (read, (id, prev)) in self.encoded_reads.iter_mut().zip(init_clustering) {
            assert_eq!(read.id, id);
            for (node, (u, cl)) in read.nodes.iter_mut().zip(prev) {
                assert_eq!(node.unit, u);
                node.cluster = cl;
            }
        }
        self.selected_chunks
            .iter_mut()
            .for_each(|c| c.cluster_num = cluster_num[&c.id]);
        // Re clustering.
        let target_units: HashSet<_> = re_cluster.keys().copied().collect();
        crate::local_clustering::local_clustering_selected(&mut self, &target_units);
        // let mut pileups: HashMap<u64, Vec<&mut _>> = HashMap::new();
        // self.encoded_reads
        //     .iter_mut()
        //     .flat_map(|r| r.nodes.iter_mut())
        //     .filter(|node| re_cluster.contains_key(&node.unit))
        //     .for_each(|node| pileups.entry(node.unit).or_default().push(node));
        // debug!("There are {} chunks to be re-clustered.", re_cluster.len());
        // // Clustering.
        // use rayon::prelude::*;
        // let coverage = self.coverage.clone();
        // let consensus_and_clusternum: HashMap<_, _> = pileups
        //     .par_iter_mut()
        //     .map(|(&unit_id, units)| {
        //         let cluster_num = re_cluster[&unit_id];
        //         let coverage = coverage.unwrap_or(units.len() as f64 / cluster_num as f64);
        //         let cluster_num = local_clustering(unit_id, units, coverage, cluster_num);
        //         (unit_id, cluster_num)
        //     })
        //     .collect();
        // for unit in self.selected_chunks.iter_mut() {
        //     if let Some(cluster_num) = consensus_and_clusternum.get(&unit.id) {
        //         unit.cluster_num = *cluster_num as usize;
        //     }
        // }
        self
    }
}

// fn local_clustering(unit_id: u64, units: &mut [&mut Node], coverage: f64, cluster_num: u8) -> u8 {
//     use crate::local_clustering::kmeans;
//     use rand::SeedableRng;
//     use rand_xoshiro::Xoroshiro128PlusPlus;
//     let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(unit_id * 23);
//     let seqs: Vec<_> = units.iter().map(|node| node.seq()).collect();
//     let mut config = kmeans::ClusteringConfig::new(100, cluster_num, coverage);
//     let start = std::time::Instant::now();
//     let (asn, _, _) = kmeans::clustering(&seqs, &mut rng, &mut config).unwrap();
//     let end = std::time::Instant::now();
//     let elapsed = (end - start).as_secs();
//     debug!(
//         "RECLUSTER\t{}\t{}\t{}\t{}\t{}",
//         unit_id,
//         elapsed,
//         seqs.len(),
//         cluster_num,
//         config.cluster_num
//     );
//     for (node, asn) in units.iter_mut().zip(asn) {
//         node.cluster = asn as u64;
//     }
//     config.cluster_num
// }

fn get_units_to_cluster(ds: &DataSet, c: &ReClusteringConfig) -> HashMap<u64, u8> {
    let mut re_cluster: HashMap<_, u8> = HashMap::new();
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use crate::assemble::ditch_graph::DitchGraph;
    let config = crate::assemble::AssembleConfig::new(c.threads, 200, false, false);
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    graph.remove_lightweight_edges(2);
    if let Some(cov) = ds.coverage {
        let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
        graph.remove_zero_copy_elements(cov, &lens, 0.3);
        graph.z_edge_selection();
        graph.remove_tips(0.5, 5);
        graph.zip_up_overclustering();
        let (node_copy_num, _) = graph.copy_number_estimation(cov, &lens);
        for (&(node, _), &copy_num) in node_copy_num.iter() {
            *re_cluster.entry(node).or_default() += copy_num as u8;
        }
    }
    // We would re-cluster the unit if the following two conditions hold:
    // 1. If the total copy number is more than the cluster number.
    // 2. If next copy number is more than 2. (If it is two, then cl = 1, which is actually not the bad thing.)
    let chunkid2cluster_num: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|u| (u.id, u.cluster_num))
        .collect();
    re_cluster.retain(
        |chunkid, &mut new_copy_num| match chunkid2cluster_num.get(chunkid) {
            Some(&cluster_num) => 2 < cluster_num && cluster_num < new_copy_num as usize,
            None => false,
        },
    );
    re_cluster
}
