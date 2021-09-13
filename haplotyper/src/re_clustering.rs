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
        debug!("RECLUSTRE\t{}", target_units.len());
        crate::local_clustering::local_clustering_selected(&mut self, &target_units);
        self
    }
}

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
