#![allow(dead_code)]
use super::find_union::FindUnion;
use definitions;
use std::collections::HashMap;
use std::collections::HashSet;
mod error_correction;
use error_correction::local_correction;
use error_correction::CorrectedRead;
pub mod path_clustering;
pub use path_clustering::path_clustering;
#[derive(Debug, Clone, Copy)]
pub struct GlobalClusteringConfig {
    pub threads: usize,
    pub k_mer: usize,
    pub min_cluster_size: usize,
}
impl GlobalClusteringConfig {
    pub fn new(threads: usize, k_mer: usize, min_cluster_size: usize) -> Self {
        Self {
            threads,
            k_mer,
            min_cluster_size,
        }
    }
}
pub trait GlobalClustering {
    fn global_clustering(self, c: &GlobalClusteringConfig) -> Self;
}

#[derive(Clone)]
pub struct DeBruijnGraph {
    k: usize,
    nodes: Vec<Node>,
    indexer: HashMap<Node, usize>,
}

impl std::fmt::Debug for DeBruijnGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (idx, node) in self.nodes.iter().enumerate() {
            writeln!(f, "{}\t{:?}", idx, node)?;
        }
        write!(f, "K:{}", self.k)
    }
}

#[derive(Clone)]
struct Node {
    occ: usize,
    edges: Vec<Edge>,
    kmer: Vec<(u64, u64)>,
    cluster: usize,
}
impl std::fmt::Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let edges: Vec<_> = self
            .edges
            .iter()
            .map(|e| format!("(->{},{})", e.to, e.weight))
            .collect();
        let kmer: Vec<_> = self
            .kmer
            .iter()
            .map(|(u, c)| format!("{}-{}", u, c))
            .collect();
        write!(
            f,
            "{}\t{}\t{}\t[{}]",
            self.cluster,
            self.occ,
            edges.join(","),
            kmer.join(",")
        )
    }
}

impl Node {
    fn new(w: &[definitions::Node]) -> Self {
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
        let (edges, occ, cluster) = (vec![], 0, 0);
        Self {
            kmer,
            edges,
            occ,
            cluster,
        }
    }
    fn from_corrected(w: &[error_correction::Unit]) -> Self {
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
        let (edges, occ, cluster) = (vec![], 0, 0);
        Self {
            kmer,
            edges,
            occ,
            cluster,
        }
    }

    fn push(&mut self, to: usize) {
        match self.edges.iter_mut().find(|e| e.to == to) {
            Some(x) => {
                x.weight += 1;
            }
            None => self.edges.push(Edge { to, weight: 1 }),
        }
    }
}

#[derive(Debug, Clone)]
struct Edge {
    to: usize,
    weight: u64,
}

use std::hash::Hasher;
impl std::hash::Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        assert!(!self.kmer.is_empty());
        // Normalize and hashing.
        if self.kmer.first().unwrap() < self.kmer.last().unwrap() {
            for (unit, cluster) in self.kmer.iter() {
                unit.hash(state);
                cluster.hash(state);
            }
        } else {
            for (unit, cluster) in self.kmer.iter().rev() {
                unit.hash(state);
                cluster.hash(state);
            }
        }
    }
}

impl std::cmp::PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        assert!(!self.kmer.is_empty());
        assert!(!other.kmer.is_empty());
        if self.kmer.len() != other.kmer.len() {
            return false;
        }
        let is_self_normed = self.kmer.first().unwrap() < self.kmer.last().unwrap();
        let is_other_normed = other.kmer.first().unwrap() < other.kmer.last().unwrap();
        match (is_self_normed, is_other_normed) {
            (false, false) | (true, true) => self.kmer.iter().zip(&other.kmer).all(|(x, y)| x == y),
            (false, true) | (true, false) => {
                self.kmer.iter().rev().zip(&other.kmer).all(|(x, y)| x == y)
            }
        }
    }
}
impl std::cmp::Eq for Node {}

impl DeBruijnGraph {
    pub fn from_encoded_reads(reads: &[definitions::EncodedRead], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            for w in read.nodes.windows(k + 1) {
                // Calc kmer
                let from = Node::new(&w[..k]);
                let to = Node::new(&w[1..]);
                // Check entry.
                let from = if !indexer.contains_key(&from) {
                    indexer.insert(from.clone(), nodes.len());
                    nodes.push(from);
                    nodes.len() - 1
                } else {
                    *indexer.get(&from).unwrap()
                };
                let to = if !indexer.contains_key(&to) {
                    indexer.insert(to.clone(), nodes.len());
                    nodes.push(to);
                    nodes.len() - 1
                } else {
                    *indexer.get(&to).unwrap()
                };
                nodes[from].occ += 1;
                nodes[to].occ += 1;
                nodes[from].push(to);
                nodes[to].push(from);
            }
        }

        Self { k, nodes, indexer }
    }
    fn from_corrected_reads(reads: &[CorrectedRead], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            for w in read.nodes.windows(k + 1) {
                // Calc kmer
                let from = Node::from_corrected(&w[..k]);
                let to = Node::from_corrected(&w[1..]);
                // Check entry.
                let from = if !indexer.contains_key(&from) {
                    indexer.insert(from.clone(), nodes.len());
                    nodes.push(from);
                    nodes.len() - 1
                } else {
                    *indexer.get(&from).unwrap()
                };
                let to = if !indexer.contains_key(&to) {
                    indexer.insert(to.clone(), nodes.len());
                    nodes.push(to);
                    nodes.len() - 1
                } else {
                    *indexer.get(&to).unwrap()
                };
                nodes[from].occ += 1;
                nodes[to].occ += 1;
                nodes[from].push(to);
                nodes[to].push(from);
            }
        }

        Self { k, nodes, indexer }
    }
    fn calc_thr(&self) -> usize {
        let counts: Vec<_> = self.nodes.iter().map(|n| n.occ).collect();
        let mean = counts.iter().sum::<usize>() / counts.len();
        mean / 2
    }
    fn calc_thr_edge(&self) -> u64 {
        let counts = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().fold(0, |x, w| x + w.weight))
            .sum::<u64>();
        let len: usize = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        counts / len as u64 / 3
    }
    fn clean_up_auto(self) -> Self {
        // let thr = self.calc_thr();
        // debug!("Removing nodes with occ less than {}", thr);
        // self.clean_up(thr)
        let thr = self.calc_thr_edge();
        debug!("Removing edges with weight less than {}", thr);
        self.clean_up_2(thr)
    }
    fn clean_up_2(mut self, thr: u64) -> Self {
        // Removing weak and non-solo edge.
        // This is because, solo-edge would be true edge
        // Even if the weight is small.
        // It is not a problem to leave consective false-edge ,
        // as such a cluster would be removed by filtering
        // clusters by their sizes.
        self.nodes
            .iter_mut()
            .filter(|node| node.edges.len() > 2)
            .for_each(|node| {
                // Remove weak edges.
                node.edges.retain(|edge| edge.weight > thr);
            });
        self
    }
    fn clean_up(self, thr: usize) -> Self {
        let mut resulting_node = vec![];
        let Self {
            k,
            nodes,
            mut indexer,
        } = self;
        let mut current_idx = 0;
        let mut map = vec![None; nodes.len()];
        for (idx, node) in nodes.into_iter().enumerate() {
            if node.occ > thr {
                map[idx] = Some(current_idx);
                if let Some(x) = indexer.get_mut(&node) {
                    *x = current_idx
                }
                current_idx += 1;
                resulting_node.push(node);
            }
        }
        resulting_node.iter_mut().for_each(|n| {
            n.edges.retain(|edge| map[edge.to].is_some());
            n.edges
                .iter_mut()
                .for_each(|edge| edge.to = map[edge.to].unwrap());
        });
        Self {
            k,
            nodes: resulting_node,
            indexer,
        }
    }
    fn assign_corrected_read(&self, read: &CorrectedRead) -> Option<usize> {
        let mut count = HashMap::<_, u32>::new();
        for w in read.nodes.windows(self.k) {
            let node = Node::from_corrected(w);
            if let Some(&idx) = self.indexer.get(&node) {
                *count.entry(self.nodes[idx].cluster).or_default() += 1;
            }
        }
        count.into_iter().max_by_key(|x| x.1).map(|x| x.0)
    }
    fn coloring(&mut self, _c: &GlobalClusteringConfig) {
        // Coloring node of the de Bruijn graph.
        // As a first try, I just color de Bruijn graph by its connected components.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut current_component = 1;
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            debug!("Find cluster. Containing {} k-mers.", count);
            if count < 10 {
                continue;
            }
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = current_component;
                }
            }
            current_component += 1;
        }
    }
    fn connections(
        &self,
        i: usize,
        j: usize,
        reads: &[CorrectedRead],
        _c: &GlobalClusteringConfig,
    ) -> usize {
        let set_i: HashSet<_> = self
            .nodes
            .iter()
            .filter(|n| n.cluster == i)
            .flat_map(|n| n.kmer.iter())
            .copied()
            .collect();
        let set_j: HashSet<_> = self
            .nodes
            .iter()
            .filter(|n| n.cluster == j)
            .flat_map(|n| n.kmer.iter())
            .copied()
            .collect();
        let connections: Vec<_> = reads
            .iter()
            .filter(|read| {
                let (mut is_i, mut is_j) = (false, false);
                for node in read.nodes.iter() {
                    let unit = (node.unit, node.cluster);
                    is_i |= set_i.contains(&unit);
                    is_j |= set_j.contains(&unit);
                }
                is_i && is_j
            })
            .filter_map(|read| {
                let last = read.nodes.last().unwrap();
                let last = (last.unit, last.cluster);
                let first = read.nodes.first().unwrap();
                let first = (first.unit, first.cluster);
                if (set_i.contains(&last) && set_j.contains(&first))
                    || set_i.contains(&first) && set_j.contains(&last)
                {
                    Some((last, first))
                } else {
                    None
                }
            })
            .collect();
        debug!("{},{}=>{:?}", i, j, connections);
        connections.len()
    }
    pub fn clustering(&self, thr: usize) -> Vec<HashSet<(u64, u64)>> {
        // Clustering de Bruijn graph.
        // As a first try, I implement very naive conneceted component analysis.
        // To this end, I use naive FindUnion Tree. In other words,
        // as I traverse nodes, I merge the two connected nodes.
        // Currently, I ignore very weak connection, i.e., connection
        // with the weight of less than 1.
        let mut fu = FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|(_, n)| n.occ > thr) {
            for edge in node.edges.iter().filter(|e| e.weight > thr as u64) {
                fu.unite(from, edge.to);
            }
        }
        let mut components = vec![];
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            debug!("Find cluster. Containing {} k-mers.", count);
            let component: HashSet<_> = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .flat_map(|node_idx| self.nodes[node_idx].kmer.iter())
                .copied()
                .collect();
            if component.len() > 1 {
                components.push(component);
            }
        }
        components
    }
    fn clustering_mcl(&self, e: i32, r: i32) -> Vec<HashSet<(u64, u64)>> {
        // Convert de Bruijn graph into an usual node-edge graph.
        let nodes = self.nodes.len();
        let edges: Vec<Vec<(usize, u64)>> = self
            .nodes
            .iter()
            .map(|node| node.edges.iter().map(|e| (e.to, e.weight)).collect())
            .collect();
        let num_edge: usize = edges.iter().map(|x| x.len()).sum::<usize>();
        debug!("Node:{}, edges:{}", nodes, num_edge);
        let graph = mcl::Graph::new(nodes, &edges);
        debug!("Constructed a MCL graph.");
        let mut components = graph.clustering(e, r);
        components.retain(|c| c.len() > 10);
        debug!("Clustered. Deduping..");
        components.sort();
        components.dedup();
        debug!("Convert into HashSet.");
        components
            .into_iter()
            .map(|comp| {
                comp.into_iter()
                    .flat_map(|idx| {
                        self.nodes
                            .iter()
                            .find(|n| match self.indexer.get(n) {
                                Some(i) => *i == idx,
                                None => false,
                            })
                            .map(|n| n.kmer.clone())
                            .unwrap_or_else(Vec::new)
                    })
                    .collect()
            })
            .collect()
    }
}

#[derive(Clone)]
struct PlugGraph {
    nodes: usize,
    edges: Vec<Vec<PlugEdge>>,
}

impl std::fmt::Debug for PlugGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Nodes:{}", self.nodes)?;
        let lines: Vec<_> = self.edges.iter().map(|e| format!("{:?}", e)).collect();
        write!(f, "{}", lines.join("\n"))
    }
}

#[derive(Clone, Default)]
struct PlugEdge {
    to: usize,
    weight: u32,
}

impl std::fmt::Debug for PlugEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(->{}:{})", self.to, self.weight)
    }
}

impl PlugGraph {
    fn from_corrected_reads(dbg: &DeBruijnGraph, reads: &[CorrectedRead]) -> Self {
        let nodes = dbg.nodes.iter().map(|n| n.cluster).max().unwrap() + 1;
        debug!("Nodes: {}", nodes);
        // let mut plugs: Vec<HashSet<(u64, u64)>> = (0..nodes).map(|_| HashSet::new()).collect();
        let plugs: Vec<HashSet<_>> = (0..nodes)
            .map(|c| {
                dbg.nodes
                    .iter()
                    .filter(|n| n.cluster == c)
                    .flat_map(|n| n.kmer.iter())
                    .copied()
                    .collect()
            })
            .collect();
        let edges: Vec<Vec<PlugEdge>> = (0..nodes)
            .map(|from| {
                (0..nodes)
                    .filter(|&to| to != from)
                    .filter_map(|to| {
                        eprintln!("PlugEdds of {}->{}", from, to);
                        let weight = reads
                            .iter()
                            .filter(|read| read.nodes.len() > 2)
                            .filter(|read| {
                                let start = {
                                    let s = &read.nodes.first().unwrap();
                                    (s.unit, s.cluster)
                                };
                                let end = {
                                    let e = &read.nodes.last().unwrap();
                                    (e.unit, e.cluster)
                                };
                                let to = &plugs[to];
                                let from = &plugs[from];
                                let connected = (to.contains(&start) && from.contains(&end))
                                    || (to.contains(&end) && from.contains(&start));
                                if connected {
                                    let nodes: String = read
                                        .nodes
                                        .iter()
                                        .map(|u| {
                                            let u = (u.unit, u.cluster);
                                            match (from.contains(&u), to.contains(&u)) {
                                                (true, true) => 'B',
                                                (true, false) => 'F',
                                                (false, true) => 'T',
                                                (false, false) => 'X',
                                            }
                                        })
                                        .collect();
                                    eprintln!("{}", nodes);
                                    let nodes: Vec<_> =
                                        read.nodes.iter().map(|u| (u.unit, u.cluster)).collect();
                                    eprintln!("{:?}", nodes);
                                }
                                connected
                            })
                            .count() as u32;
                        if weight == 0 {
                            None
                        } else {
                            Some(PlugEdge { to, weight })
                        }
                    })
                    .collect()
            })
            .collect();
        Self { nodes, edges }
    }
    // Mapping: de Bruijn graph component -> Aggregated cluster.
    fn clustering(&self) -> HashMap<usize, usize> {
        debug!("Clustering {:?}", self);
        let mut fu = FindUnion::new(self.nodes);
        for (from, edges) in self.edges.iter().enumerate() {
            for edge in edges.iter().filter(|e| e.weight > 1) {
                fu.unite(from, edge.to).unwrap();
            }
        }
        let mut mapping = HashMap::new();
        let mut current_component = 0;
        for cluster in 0..self.nodes {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes)
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            debug!("Find cluster. Containing {} de Bruijn SCC.", count);
            for node in 0..self.nodes {
                if fu.find(node).unwrap() == cluster {
                    mapping.insert(node, current_component);
                }
            }
            current_component += 1;
        }
        debug!("Resulting in {} clusters", current_component);
        debug!("{:?}", mapping);
        mapping
    }
}

impl GlobalClustering for definitions::DataSet {
    fn global_clustering(mut self, c: &GlobalClusteringConfig) -> Self {
        if log_enabled!(log::Level::Debug) {
            let length: Vec<_> = self.encoded_reads.iter().map(|r| r.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&length);
            let tot = length.iter().sum::<usize>();
            eprintln!("Read({}){}\n{}", length.len(), tot, hist.format(20, 40));
        }
        let reads = local_correction(&self);
        debug!("Corrected reads.");
        if log_enabled!(log::Level::Debug) {
            let length: Vec<_> = reads.iter().map(|r| r.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&length);
            let tot = length.iter().sum::<usize>();
            eprintln!("Read({}){}\n{}", length.len(), tot, hist.format(20, 40));
        }
        let graph = DeBruijnGraph::from_corrected_reads(&reads, c.k_mer);
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("Node({}-mer) occurences\n{}", c.k_mer, hist.format(20, 40));
        }
        let mut graph = graph.clean_up_auto();
        //let mut graph = graph;
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("Node({}-mer) occurences\n{}", c.k_mer, hist.format(20, 40));
        }
        // let mut components: Vec<HashSet<(u64, u64)>> = graph.clustering(0);
        // components.retain(|cmp| cmp.len() > c.min_cluster_size);
        // let components = graph.clustering_mcl(2, 2);
        // let components = merge(components, &self.encoded_reads);
        // debug!("Resulting in {} clusters.", components.len());
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
        // let mapping = PlugGraph::from_corrected_reads(&graph, &reads).clustering();
        // graph
        //     .nodes
        //     .iter_mut()
        //     .for_each(|n| n.cluster = mapping[&n.cluster]);
        let component_num = graph.nodes.iter().map(|n| n.cluster).max().unwrap() + 1;
        debug!("Resulting in {} clusters.", component_num);
        let components: Vec<HashSet<_>> = (0..component_num)
            .map(|c| {
                graph
                    .nodes
                    .iter()
                    .filter(|n| n.cluster == c)
                    .flat_map(|n| n.kmer.iter())
                    .copied()
                    .collect()
            })
            .collect();
        for i in 0..component_num {
            for j in i + 1..component_num {
                let count = components[i].intersection(&components[j]).count();
                debug!("{}\t{}\t{}", i, j, count);
            }
        }
        let mut count: HashMap<_, usize> = HashMap::new();
        let assignments: Vec<_> = reads
            .into_iter()
            .filter_map(|read| {
                let id = read.id;
                graph.assign_corrected_read(&read).map(|cluster| {
                    *count.entry(cluster).or_default() += 1;
                    definitions::Assignment { id, cluster }
                })
            })
            .collect();
        // let assignments: Vec<_> = self
        //     .encoded_reads
        //     .iter()
        //     .map(|read| {
        //         let id = read.id;
        //         let cluster = components
        //             .iter()
        //             .map(|c| sim(c, read))
        //             .enumerate()
        //             .max_by_key(|x| x.1)
        //             .unwrap()
        //             .0;
        //         *count.entry(cluster).or_default() += 1;
        //         definitions::Assignment { id, cluster }
        //     })
        //     .collect();
        self.assignments = assignments;
        if log_enabled!(log::Level::Debug) {
            eprintln!("Cluster\tCount");
            for (cl, count) in count {
                eprintln!("{}\t{}", cl, count);
            }
        }
        self
    }
}

fn merge(
    mut components: Vec<HashSet<(u64, u64)>>,
    reads: &[definitions::EncodedRead],
) -> Vec<HashSet<(u64, u64)>> {
    'outer: loop {
        let len = components.len();
        for i in 0..len {
            for j in (i + 1)..len {
                let int_reads = reads
                    .iter()
                    .filter(|r| sim(&components[i], r) > 0 && sim(&components[j], r) > 0)
                    .count();
                if int_reads > 3 {
                    debug!(
                        "Merge {} and {} as they share {} reads in common.",
                        i, j, int_reads
                    );
                    let c: Vec<_> = components[j].drain().collect();
                    components[i].extend(c);
                    continue 'outer;
                }
            }
        }
        break;
    }
    components.into_iter().filter(|c| !c.is_empty()).collect()
}

fn sim(cluster: &HashSet<(u64, u64)>, read: &definitions::EncodedRead) -> u32 {
    read.nodes
        .iter()
        .filter(|n| cluster.contains(&(n.unit, n.cluster)))
        .count() as u32
}

#[cfg(test)]
mod tests {
    impl Node {
        fn new_raw(w: &[(u64, u64)]) -> Self {
            let first = w.first().unwrap();
            let last = w.last().unwrap();
            let kmer: Vec<_> = if first < last {
                w.iter().copied().collect()
            } else {
                w.iter().rev().copied().collect()
            };
            let (edges, occ, cluster) = (vec![], 0, 0);
            Self {
                kmer,
                edges,
                occ,
                cluster,
            }
        }
    }
    impl DeBruijnGraph {
        fn new(reads: &[Vec<(u64, u64)>], k: usize) -> Self {
            let (mut nodes, mut indexer) = (vec![], HashMap::new());
            for read in reads {
                for w in read.windows(k + 1) {
                    // Calc kmer
                    let from = Node::new_raw(&w[..k]);
                    let to = Node::new_raw(&w[1..]);
                    // Check entry.
                    let from = if !indexer.contains_key(&from) {
                        indexer.insert(from.clone(), nodes.len());
                        nodes.push(from);
                        nodes.len() - 1
                    } else {
                        *indexer.get(&from).unwrap()
                    };
                    let to = if !indexer.contains_key(&to) {
                        indexer.insert(to.clone(), nodes.len());
                        nodes.push(to);
                        nodes.len() - 1
                    } else {
                        *indexer.get(&to).unwrap()
                    };
                    nodes[from].occ += 1;
                    nodes[to].occ += 1;
                    nodes[from].push(to);
                    nodes[to].push(from);
                }
            }
            Self { k, nodes, indexer }
        }
    }
    #[derive(Clone, Copy, Debug)]
    struct TestConfig {
        cl: usize,
        num: usize,
        fail: f64,
        skip: f64,
        max_len: usize,
        min_len: usize,
        unit_len: usize,
    }
    use super::*;
    use rand::Rng;
    fn gen_dataset<R: Rng>(r: &mut R, conf: TestConfig) -> (Vec<Vec<(u64, u64)>>, Vec<usize>) {
        let TestConfig {
            cl,
            num,
            fail,
            skip,
            max_len,
            min_len,
            unit_len,
        } = conf;
        let mut answer = vec![];
        let mut reads = vec![];
        for i in 0..num {
            let cluster = (i % cl) as u64;
            let len = r.gen::<usize>() % (max_len - min_len) + min_len;
            let start = r.gen::<usize>() % (unit_len - len);
            let units: Vec<_> = if r.gen_bool(0.5) {
                (start..=start + len).collect()
            } else {
                let start = start + len;
                (start - len..=start).rev().collect()
            };
            let mut read = vec![];
            for unit in units {
                if r.gen_bool(skip) {
                    continue;
                } else if r.gen_bool(fail) {
                    let cluster = r.gen::<u64>() % cl as u64;
                    read.push((unit as u64, cluster));
                } else {
                    read.push((unit as u64, cluster));
                }
            }
            answer.push(cluster as usize);
            reads.push(read);
        }
        (reads, answer)
    }
    #[test]
    fn construction_test() {
        let read = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.nodes.len(), 3, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.nodes.len(), 3, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.nodes.len(), 5, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 1), (0, 2), (0, 3), (1, 4), (0, 5), (0, 6), (0, 7)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.nodes.len(), 8, "{:?}", graph);
    }
    #[test]
    fn clustering_test() {
        let read = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
        let read = vec![
            vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)],
            vec![(0, 1), (0, 2), (0, 3), (1, 4), (0, 5), (0, 6), (0, 7)],
            vec![(0, 5), (0, 4), (0, 3), (0, 2), (0, 1)],
            vec![(0, 3), (0, 4), (0, 5), (0, 6), (0, 7)],
        ];
        let graph = DeBruijnGraph::new(&read, 3);
        assert_eq!(graph.clustering(0).len(), 1, "{:?}", graph);
    }
    #[test]
    fn clustering_test_2() {
        // Case1
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 2
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 0), (1, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 3
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(0, 1), (1, 1), (2, 2), (3, 1), (4, 1), (5, 1)],
            vec![(4, 1), (3, 1), (2, 2), (1, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        // Case 4
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(0, 0), (1, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
    }
    #[test]
    fn clustering_test_3() {
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1), (1, 0)],
            vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 0), (1, 1), (0, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3).clean_up(1);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
        let reads = vec![
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(5, 0), (4, 0), (3, 0), (2, 0), (1, 0), (0, 0)],
            vec![(1, 0), (2, 0), (3, 0), (4, 0), (5, 0)],
            vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(5, 1), (4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
            vec![(1, 1), (2, 1), (3, 1), (4, 1), (5, 1)],
            vec![(0, 1), (1, 1), (2, 0), (3, 0), (4, 0), (5, 1)],
            vec![(5, 1), (4, 0), (3, 0), (2, 0), (1, 1), (0, 1)],
        ];
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up(1);
        assert_eq!(graph.clustering(0).len(), 2, "\n{:?}", graph);
    }
    #[test]
    fn path_clustering_test_large_noisy() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 2000,
            fail: 0.02,
            skip: 0.02,
            max_len: 40,
            min_len: 5,
            unit_len: 500,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(1);
        components.retain(|c| c.len() > 100);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    #[test]
    fn path_clustering_test_large_hard() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 5000,
            fail: 0.02,
            skip: 0.02,
            max_len: 40,
            min_len: 10,
            unit_len: 800,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(1);
        components.retain(|cl| cl.len() > 100);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    #[test]
    fn path_clustering_test_short() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 200,
            fail: 0.01,
            skip: 0.01,
            max_len: 20,
            min_len: 10,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let mut components = graph.clustering(0);
        components.retain(|c| c.len() > 20);
        assert_eq!(components.len(), 2);
        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    fn path_clustering_test_mcl() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 100,
            fail: 0.05,
            skip: 0.05,
            max_len: 15,
            min_len: 8,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let components = graph.clustering_mcl(10, 2);
        assert_eq!(components.len(), 2, "{:?}", components,);
        {
            for c in components.iter() {
                let mut c: Vec<_> = c.iter().copied().collect();
                c.sort();
                eprintln!("{:?}", c);
            }
        }

        let preds: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .enumerate()
                    .max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(preds.len(), answer.len());
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| {
                eprintln!("{}\t{}", ans, pred);
                ans == pred
            })
            .count();
        eprintln!("{}/{}", correct, reads.len());
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
    fn path_clustering_test_mcl_with() {
        use rand::SeedableRng;
        use rand_xoshiro::Xoshiro256StarStar;
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
        let conf = TestConfig {
            cl: 2,
            num: 200,
            fail: 0.05,
            skip: 0.05,
            max_len: 20,
            min_len: 10,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let components = graph.clustering_mcl(10, 2);
        assert_eq!(components.len(), 2, "{:?}", components,);
        {
            for c in components.iter() {
                let mut c: Vec<_> = c.iter().copied().collect();
                c.sort();
                eprintln!("{:?}", c);
            }
        }
        let init: Vec<_> = reads
            .iter()
            .filter_map(|read| {
                let counts: Vec<_> = components
                    .iter()
                    .map(|c| read.iter().filter(|x| c.contains(x)).count())
                    .collect();
                counts.into_iter().enumerate().max_by_key(|x| x.1)
            })
            .map(|x| x.0)
            .collect();
        assert_eq!(init.len(), answer.len());
        let preds = {
            let units: HashSet<_> = reads.iter().flat_map(|e| e).copied().collect();
            let map: HashMap<(u64, u64), usize> =
                units.into_iter().enumerate().map(|(x, y)| (y, x)).collect();
            let reads: Vec<Vec<_>> = reads
                .iter()
                .map(|read| read.iter().map(|u| map[u]).collect())
                .collect();
            path_clustering::path_clustering(&reads, &init)
        };
        eprintln!("------------------");
        {
            let max = *preds.iter().max().unwrap();
            for cl in 0..=max {
                let mut component: Vec<_> = reads
                    .iter()
                    .zip(preds.iter())
                    .filter(|&(_, &x)| x == cl)
                    .flat_map(|(p, _)| p.iter())
                    .collect::<HashSet<_>>()
                    .into_iter()
                    .collect();
                component.sort();
                eprintln!("{:?}", component);
            }
        }
        let correct = answer
            .iter()
            .zip(preds.iter())
            .filter(|&(ans, pred)| ans == pred)
            .count();
        eprintln!("{}/{}", correct, reads.len());
        for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
            eprintln!("{}\t{}\t{}", ans, assign, read.len());
        }
        let correct = correct.max(reads.len() - correct);
        assert!(correct > reads.len() * 8 / 10);
    }
}
