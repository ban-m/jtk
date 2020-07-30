use de_bruijn_graph::*;
use std::collections::HashMap;
struct ReadWrapper<'a>(&'a definitions::EncodedRead);

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
                let first = {
                    let f = w.first().unwrap();
                    (f.unit, f.cluster)
                };
                let last = {
                    let l = w.last().unwrap();
                    (l.unit, l.cluster)
                };
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
        let reads: Vec<_> = self.encoded_reads.iter().map(ReadWrapper::new).collect();
        let graph = DeBruijnGraph::from(&reads, c.k_mer);
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
        graph.resolve_crossings(&reads);
        graph.resolve_bubbles(&reads);
        if log_enabled!(log::Level::Debug) {
            let mut count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            count.sort();
            let top_15: Vec<_> = count.iter().rev().take(15).collect();
            debug!("Top 15 Occurences:{:?}", top_15);
            let last = **top_15.last().unwrap();
            let count: Vec<_> = count.into_iter().filter(|&x| x < last).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("The rest({}-mer) nodes\n{}", c.k_mer, hist.format(20, 40));
            let mut count: HashMap<_, u32> = HashMap::new();
            for n in graph.nodes.iter() {
                *count.entry(n.edges.len()).or_default() += 1;
            }
            let mut count: Vec<_> = count.into_iter().collect();
            count.sort_by_key(|x| x.0);
            eprintln!("Degree Count\n{:?}", count);
        }
        graph.coloring();
        let component_num = graph.nodes.iter().filter_map(|n| n.cluster).max().unwrap() + 1;
        debug!("Resulting in {} clusters.", component_num);
        let mut count: HashMap<_, usize> = HashMap::new();
        let assignments: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                let id = read.id();
                graph.assign_read(read).map(|cluster| {
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
