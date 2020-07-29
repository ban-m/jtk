use super::find_union::FindUnion;
use std::collections::HashMap;
use std::collections::HashSet;
mod de_bruijn_graph;
mod error_correction;
use de_bruijn_graph::*;
// use error_correction::local_correction;
use error_correction::CorrectedRead;
pub mod path_clustering;
pub use path_clustering::path_clustering;
#[derive(Debug, Clone, Copy)]
pub struct GlobalClusteringConfig {
    pub k_mer: usize,
    pub min_cluster_size: usize,
    pub mat_score: i32,
    pub mismat_score: i32,
    pub gap_score: i32,
}

impl GlobalClusteringConfig {
    pub fn new(
        k_mer: usize,
        min_cluster_size: usize,
        mat_score: i32,
        mismat_score: i32,
        gap_score: i32,
    ) -> Self {
        Self {
            k_mer,
            min_cluster_size,
            mat_score,
            mismat_score,
            gap_score,
        }
    }
}
pub trait GlobalClustering {
    fn global_clustering(self, c: &GlobalClusteringConfig) -> Self;
}

impl GlobalClustering for definitions::DataSet {
    fn global_clustering(mut self, c: &GlobalClusteringConfig) -> Self {
        if log_enabled!(log::Level::Debug) {
            let length: Vec<_> = self.encoded_reads.iter().map(|r| r.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&length);
            let tot = length.iter().sum::<usize>();
            eprintln!("Read({}){}\n{}", length.len(), tot, hist.format(20, 40));
        }
        // let reads = local_correction(&self, c);
        // debug!("Corrected reads.");
        // if log_enabled!(log::Level::Debug) {
        //     let length: Vec<_> = reads.iter().map(|r| r.nodes.len()).collect();
        //     let hist = histgram_viz::Histgram::new(&length);
        //     let tot = length.iter().sum::<usize>();
        //     eprintln!("Read({}){}\n{}", length.len(), tot, hist.format(20, 40));
        // }
        //let graph = DeBruijnGraph::from_corrected_reads(&reads, c.k_mer);
        let graph = DeBruijnGraph::from_encoded_reads(&self.encoded_reads, c.k_mer);
        if log_enabled!(log::Level::Debug) {
            let mut count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            count.sort();
            let top_15: Vec<_> = count.iter().rev().take(15).collect();
            debug!("Top 15 Occurences:{:?}", top_15);
            let last = **top_15.last().unwrap();
            let count: Vec<_> = count.into_iter().filter(|&x| x < last).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("The rest({}-mer) nodes\n{}", c.k_mer, hist.format(20, 40));
        }
        let mut graph = graph.clean_up_auto();
        //let mut graph = graph;
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("Node({}-mer) occurences\n{}", c.k_mer, hist.format(20, 40));
        }
        if log_enabled!(log::Level::Debug) {
            let mut count: HashMap<_, u32> = HashMap::new();
            for n in graph.nodes.iter() {
                *count.entry(n.edges.len()).or_default() += 1;
            }
            let mut count: Vec<_> = count.into_iter().collect();
            count.sort_by_key(|x| x.0);
            eprintln!("Degree Count\n{:?}", count);
        }
        graph.coloring(c);
        let component_num = graph.nodes.iter().map(|n| n.cluster).max().unwrap() + 1;
        debug!("Resulting in {} clusters.", component_num);
        let mut count: HashMap<_, usize> = HashMap::new();
        let assignments: Vec<_> = self
            .encoded_reads
            .iter()
            .filter_map(|read| {
                let id = read.id;
                graph.assign_encoded_read(&read).map(|cluster| {
                    *count.entry(cluster).or_default() += 1;
                    definitions::Assignment { id, cluster }
                })
            })
            .collect();
        self.assignments = assignments;
        if log_enabled!(log::Level::Debug) {
            let mut count: Vec<_> = count.into_iter().collect();
            count.sort_by_key(|x| x.0);
            debug!("Cluster\tCount");
            for (cl, count) in count {
                debug!("{}\t{}", cl, count);
            }
        }
        self
    }
}

pub fn remove_collapsed_units(mut reads: Vec<CorrectedRead>) -> Vec<CorrectedRead> {
    let unit_counts = {
        let mut counts: HashMap<_, usize> = HashMap::new();
        for read in reads.iter() {
            for node in read.nodes.iter() {
                *counts.entry(node.unit).or_default() += 1;
            }
        }
        counts
    };
    if log_enabled!(log::Level::Debug) {
        let counts: Vec<_> = unit_counts.values().copied().collect();
        let hist = histgram_viz::Histgram::new(&counts);
        eprintln!("Unit Histgram:{}\n{}", counts.len(), hist.format(20, 40));
    }
    let mean = unit_counts.values().copied().sum::<usize>() / unit_counts.len();
    debug!("Removing units having occurence more than {}", mean / 2);
    debug!("And single component.");
    let mut char_count: HashMap<u64, HashMap<u64, usize>> = HashMap::new();
    for read in reads.iter() {
        for node in read.nodes.iter() {
            *char_count
                .entry(node.unit)
                .or_default()
                .entry(node.cluster)
                .or_default() += 1;
        }
    }
    // Determine the units to be discarded.
    let discarded: HashSet<u64> = char_count
        .iter()
        .filter_map(|(unit, counts)| {
            let sum = counts.values().copied().sum::<usize>();
            let num: Vec<_> = counts.values().collect();
            debug!("{:?}->{}", num, counts.len() <= 1 && sum > mean / 2);
            if counts.len() <= 1 && sum > mean / 2 {
                Some(*unit)
            } else {
                None
            }
        })
        .collect();
    debug!(
        "Discarding {} units out of {}.",
        discarded.len(),
        char_count.len()
    );
    let count = reads
        .iter()
        .map(|read| {
            read.nodes
                .iter()
                .filter(|n| discarded.contains(&n.unit))
                .count()
        })
        .sum::<usize>();
    debug!("Dropping {} units in the dataset in total", count);
    // Discard collupsed units.
    let total = reads.iter().map(|r| r.nodes.len()).sum::<usize>();
    reads
        .iter_mut()
        .for_each(|read| read.nodes.retain(|n| !discarded.contains(&n.unit)));
    let total_after = reads.iter().map(|r| r.nodes.len()).sum::<usize>();
    debug!("{}->{}", total, total_after);
    reads
}
