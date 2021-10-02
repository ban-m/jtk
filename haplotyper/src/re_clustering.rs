use definitions::*;
use std::collections::HashMap;
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
    /// Note that it would not do clustering on chunks with estimated copy-number equal to 2.
    /// As these chunks are already "tried and failed."
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
        // Polish clusering.
        use crate::em_correction::ClusteringCorrection;
        self = self.correct_clustering_em(c.repeat_num, c.coverage_thr, true);
        // The set of the ((unit,cluster), (copy_num, offset)) to be clustered into `cluster_num` clusters.
        let re_cluster = get_units_to_cluster(&self, c);
        // Allocate them in HashMap or recover the initial clustering.
        let mut to_clustered_nodes: HashMap<_, Vec<&mut Node>> =
            re_cluster.keys().map(|&node| (node, Vec::new())).collect();
        for (read, (id, init_clusters)) in self.encoded_reads.iter_mut().zip(init_clustering) {
            assert_eq!(read.id, id);
            for (node, (unit, cl)) in read.nodes.iter_mut().zip(init_clusters) {
                assert_eq!(node.unit, unit);
                match to_clustered_nodes.get_mut(&(node.unit, node.cluster)) {
                    Some(bucket) => bucket.push(node),
                    None => node.cluster = cl,
                }
            }
        }
        let coverage = self.coverage.unwrap_or(30f64);
        // Do clustering, change the clustering information.
        // TODO: Is it OK to define another units for each re-clustered units?
        // Or, is it better to retain the original units,
        // and try to assign the squished cluster number after clustering?
        // The latter is more "flexible", as the information is not lost.
        // However, keep in mind that, by retaining the previous units,
        // we SHOULD take consensus on each (unit,cluter) in the assembly step,
        // newly, from scratch. Is that very dangerous?
        // List of ((unit,cluster), clustered_num).
        use rayon::prelude::*;
        debug!("RECLUSTERE\t{}", re_cluster.len());
        let scores: HashMap<_, _> = to_clustered_nodes
            .par_iter_mut()
            .map(|(key, nodes)| {
                use rand::SeedableRng;
                use rand_xoshiro::Xoshiro256StarStar;
                let (copy_num, offset) = re_cluster[key];
                let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(key.0 * 2312);
                let seqs: Vec<_> = nodes.iter().map(|node| node.seq()).collect();
                use crate::local_clustering::kmeans;
                let config = kmeans::ClusteringConfig::new(100, copy_num, coverage);
                let start = std::time::Instant::now();
                let cov = seqs.len();
                let (asn, consensus, score) = kmeans::clustering(&seqs, &mut rng, &config).unwrap();
                for (node, asn) in nodes.iter_mut().zip(asn) {
                    node.cluster = asn as u64 + offset;
                }
                let end = std::time::Instant::now();
                let elapsed = (end - start).as_secs();
                let len = consensus.len();
                debug!(
                    "RECLUSTER\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                    key.0, key.1, elapsed, len, score, cov
                );
                (key.0, score)
            })
            .collect();
        // Re-naming the newly clustered region.
        // Allocate
        let mut updated_nodes: HashMap<_, Vec<&mut Node>> =
            re_cluster.keys().map(|x| (x.0, Vec::new())).collect();
        for node in self
            .encoded_reads
            .iter_mut()
            .flat_map(|r| r.nodes.iter_mut())
        {
            if let Some(bucket) = updated_nodes.get_mut(&node.unit) {
                bucket.push(node);
            }
        }
        // Renaming.
        for (_, buckets) in updated_nodes {
            let mut new_cluster_id: HashMap<u64, u64> = HashMap::new();
            for node in buckets.iter() {
                if !new_cluster_id.contains_key(&node.cluster) {
                    let new_id = new_cluster_id.len() as u64;
                    new_cluster_id.insert(node.cluster, new_id);
                }
            }
            for mut node in buckets {
                node.cluster = new_cluster_id[&node.cluster];
            }
        }
        // Update the estimated copy-number of the chunks.
        let increased_copy_number: HashMap<_, usize> = {
            let mut increase = HashMap::new();
            for (&(node, _), &(copy_num, _)) in re_cluster.iter() {
                // This cluster would be copy_num clusters, thus, the increased copy number is copy_num -1.
                *increase.entry(node).or_default() += (copy_num - 1) as usize;
            }
            increase
        };
        for chunk in self.selected_chunks.iter_mut() {
            chunk.cluster_num += increased_copy_number.get(&chunk.id).unwrap_or(&0);
            chunk.score += scores.get(&chunk.id).unwrap_or(&0f64);
        }
        // Update the estimated copy-number of the chunks.
        // ID of the cluster -> total excess of the clustering.
        // let re_cluster = get_units_to_cluster(&self, c);
        // // Recover the initial clustering.
        // for (read, (id, prev)) in self.encoded_reads.iter_mut().zip(init_clustering) {
        //     assert_eq!(read.id, id);
        //     for (node, (u, cl)) in read.nodes.iter_mut().zip(prev) {
        //         assert_eq!(node.unit, u);
        //         node.cluster = cl;
        //     }
        // }
        // self.selected_chunks
        //     .iter_mut()
        //     .for_each(|c| c.cluster_num = cluster_num[&c.id]);
        // // Re clustering.
        // let target_units: HashSet<_> = re_cluster.keys().copied().collect();
        // debug!("RECLUSTRE\t{}", target_units.len());
        // crate::local_clustering::local_clustering_selected(&mut self, &target_units);
        // self
        self
    }
}

fn get_units_to_cluster(ds: &DataSet, c: &ReClusteringConfig) -> HashMap<(u64, u64), (u8, u64)> {
    let current_copy_num: HashMap<u64, usize> = ds
        .selected_chunks
        .iter()
        .map(|u| (u.id, u.cluster_num))
        .collect();
    let reads: Vec<_> = ds.encoded_reads.iter().collect();
    use crate::assemble::ditch_graph::DitchGraph;
    let config = crate::assemble::AssembleConfig::new(c.threads, 200, false, false, 6);
    let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
    // TODO: Parametrize here.
    graph.remove_lightweight_edges(2, true);
    let coverage = ds.coverage.unwrap_or(30f64);
    let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
    graph.assign_copy_number(coverage, &lens);
    graph.remove_zero_copy_elements(&lens, 0.3);
    let (node_copy_num, _) = graph.copy_number_estimation(coverage, &lens);
    // If the original copy number is 2 and the re-estimated copy number is 2,
    // then, that region is "strong" diplotig, i.e., we already tried and failed the clustering.
    node_copy_num
        .iter()
        .filter(|&((unit, _), &copy_num)| 1 < copy_num && 2 < current_copy_num[unit])
        .map(|(&(unit, cluster), &copy_num)| {
            let old_cluster_num = current_copy_num[&unit];
            let offset = old_cluster_num as u64 * cluster;
            ((unit, cluster), (copy_num as u8, offset))
        })
        .collect()
}

// fn get_units_to_cluster(
//     ds: &DataSet,
//     current_copy_num: &HashMap<u64, usize>,
//     c: &ReClusteringConfig,
// ) -> HashMap<u64, u8> {
//     let mut re_cluster: HashMap<_, u8> = HashMap::new();
//     let reads: Vec<_> = ds.encoded_reads.iter().collect();
//     use crate::assemble::ditch_graph::DitchGraph;
//     let config = crate::assemble::AssembleConfig::new(c.threads, 200, false, false, 6);
//     let mut graph = DitchGraph::new(&reads, Some(&ds.selected_chunks), &config);
//     // TODO: Parametrize here.
//     graph.remove_lightweight_edges(2);
//     if let Some(cov) = ds.coverage {
//         let lens: Vec<_> = ds.raw_reads.iter().map(|x| x.seq().len()).collect();
//         graph.remove_zero_copy_elements(cov, &lens, 0.3);
//         let (node_copy_num, _) = graph.copy_number_estimation(cov, &lens);
//         for (&(node, _), &copy_num) in node_copy_num.iter() {
//             *re_cluster.entry(node).or_default() += copy_num as u8;
//         }
//     }
//     // We would re-cluster the unit if the following two conditions hold:
//     // 1. If the total copy number is more than the cluster number.
//     // 2. If next copy number is more than 2. (If it is two, then cl = 1, which is actually not the bad thing.)
//     let chunkid2cluster_num: HashMap<_, _> = ds
//         .selected_chunks
//         .iter()
//         .map(|u| (u.id, u.cluster_num))
//         .collect();
//     re_cluster.retain(
//         |chunkid, &mut new_copy_num| match chunkid2cluster_num.get(chunkid) {
//             Some(&cluster_num) => 2 < cluster_num && cluster_num < new_copy_num as usize,
//             None => false,
//         },
//     );
//     re_cluster
// }
