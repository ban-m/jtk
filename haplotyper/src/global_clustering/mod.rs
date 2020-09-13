#![allow(dead_code)]
pub mod error_correction;
use super::unit_correlation;
use de_bruijn_graph::*;
use definitions::DataSet;
use std::collections::{HashMap, HashSet};
struct ReadWrapper<'a>(&'a definitions::EncodedRead);
struct FilteredRead {
    id: u64,
    nodes: Vec<(u64, u64)>,
}

impl FilteredRead {
    fn id(&self) -> u64 {
        self.id
    }
    fn new(x: &definitions::EncodedRead, units: &HashSet<u64>) -> Self {
        let nodes: Vec<_> = x
            .nodes
            .iter()
            .filter_map(|n| {
                if units.contains(&n.unit) {
                    None
                } else {
                    Some((n.unit, n.cluster))
                }
            })
            .collect();
        let id = x.id;
        Self { nodes, id }
    }
}

impl IntoDeBruijnNodes for FilteredRead {
    fn into_de_bruijn_nodes(&self, k: usize) -> Vec<Node> {
        self.nodes
            .windows(k)
            .map(|w| {
                let first = w.first().unwrap().0;
                let last = w.last().unwrap().0;
                let kmer: Vec<_> = if first < last {
                    w.to_vec()
                } else {
                    w.iter().rev().copied().collect()
                };
                Node::new(kmer)
            })
            .collect()
    }
}

impl<'a> ReadWrapper<'a> {
    fn new(x: &'a definitions::EncodedRead) -> Self {
        Self(x)
    }
    fn id(&self) -> u64 {
        self.0.id
    }
}

impl<'a> IntoDeBruijnNodes for ReadWrapper<'a> {
    fn into_de_bruijn_nodes(&self, k: usize) -> Vec<Node> {
        self.0
            .nodes
            .windows(k)
            .map(|w| {
                let first = w.first().unwrap().unit;
                let last = w.last().unwrap().unit;
                let kmer: Vec<_> = if first < last {
                    w.iter().map(|n| (n.unit, n.cluster)).collect()
                } else {
                    w.iter().rev().map(|n| (n.unit, n.cluster)).collect()
                };
                Node::new(kmer)
            })
            .collect()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GlobalClusteringConfig {
    pub k_mer: usize,
    pub min_cluster_size: usize,
    pub mat_score: i32,
    pub mismat_score: i32,
    pub gap_score: i32,
    pub p_value: f64,
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
            p_value: 0.01,
        }
    }
}
pub trait GlobalClustering {
    fn global_clustering(self, c: &GlobalClusteringConfig) -> Self;
}

impl GlobalClustering for definitions::DataSet {
    fn global_clustering(mut self, c: &GlobalClusteringConfig) -> Self {
        if log_enabled!(log::Level::Debug) {
            debug!("------Raw reads------");
            let reads: Vec<_> = self.encoded_reads.iter().map(ReadWrapper::new).collect();
            let graph = DeBruijnGraph::from(&reads, c.k_mer);
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("The {}-mer nodes\n{}", c.k_mer, hist.format(40, 20));
            eprintln!("Max occ:{}", count.iter().max().unwrap());
            let mut count: HashMap<_, u32> = HashMap::new();
            for n in graph.nodes.iter() {
                *count.entry(n.edges.len()).or_default() += 1;
            }
            let mut count: Vec<_> = count.into_iter().collect();
            count.sort_by_key(|x| x.0);
            eprintln!("Degree Count\n{:?}", count);
        }
        let clustered_units = super::unit_correlation::select_uninformative_units(&self, 0.01);
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                if !clustered_units[&node.unit] {
                    node.cluster = 0;
                }
            }
        }
        let reads = error_correction::local_correction(&self, c);
        // let reads: Vec<_> = self.encoded_reads.iter().map(ReadWrapper::new).collect();
        let mut graph = DeBruijnGraph::from(&reads, c.k_mer);
        if log_enabled!(log::Level::Debug) {
            debug!("------Corrected reads------");
            let mut counts: HashMap<Vec<u64>, usize> = HashMap::new();
            for node in graph.nodes.iter() {
                let kmer: Vec<_> = node.kmer.iter().map(|n| n.0).collect();
                *counts.entry(kmer).or_default() += node.occ;
            }
            let counts: Vec<_> = counts.values().copied().collect();
            let hist = histgram_viz::Histgram::new(&counts);
            eprintln!("Node occ:\n{}", hist.format(40, 20));
            eprintln!("Max occ:{}", counts.iter().max().unwrap());
        }
        graph.clean_up_auto();
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("The {}-mer nodes\n{}", c.k_mer, hist.format(40, 20));
            eprintln!("Max occ:{}", count.iter().max().unwrap());
            let mut count: HashMap<_, u32> = HashMap::new();
            for n in graph.nodes.iter() {
                *count.entry(n.edges.len()).or_default() += 1;
            }
            let mut count: Vec<_> = count.into_iter().collect();
            count.sort_by_key(|x| x.0);
            eprintln!("Degree Count\n{:?}", count);
        }
        let component_num = graph.coloring(10);
        debug!("Resulting in {} clusters.", component_num);
        let mut count: HashMap<_, usize> = HashMap::new();
        let assignments: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                let id = read.id();
                let cluster = graph
                    .assign_read(read)
                    .or_else(|| graph.assign_read_by_unit(read));
                cluster.map(|cluster| {
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

pub fn filter_uninformative_units(
    dataset: &DataSet,
    c: &GlobalClusteringConfig,
) -> Vec<error_correction::CorrectedRead> {
    let is_ok = unit_correlation::select_uninformative_units(dataset, c.p_value);
    let removed = is_ok.iter().filter(|x| !x.1).count();
    debug!("Removed {} units from {} units.", removed, is_ok.len());
    if log_enabled!(log::Level::Trace) {
        let mut is_ok: Vec<_> = is_ok.iter().collect();
        is_ok.sort_by_key(|x| x.0);
        for (unit, flag) in is_ok {
            trace!("UNIT\t{}\t{}", unit, flag);
        }
    }
    dataset
        .encoded_reads
        .iter()
        .map(|eread| {
            let id = eread.id;
            let nodes: Vec<_> = eread
                .nodes
                .iter()
                .filter(|n| is_ok[&n.unit])
                .map(|n| (n.unit, n.cluster))
                .map(|(unit, cluster)| error_correction::Unit { unit, cluster })
                .collect();
            error_correction::CorrectedRead { id, nodes }
        })
        .collect()
}
