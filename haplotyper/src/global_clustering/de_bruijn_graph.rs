use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct DeBruijnGraph {
    k: usize,
    pub nodes: Vec<Node>,
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
pub struct Node {
    pub occ: usize,
    pub edges: Vec<Edge>,
    kmer: Vec<(u64, u64)>,
    pub cluster: usize,
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
    #[allow(dead_code)]
    pub fn new(w: &[definitions::Node]) -> Self {
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
    pub fn from_corrected(w: &[super::error_correction::Unit]) -> Self {
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

    pub fn push(&mut self, to: usize) {
        match self.edges.iter_mut().find(|e| e.to == to) {
            Some(x) => {
                x.weight += 1;
            }
            None => self.edges.push(Edge { to, weight: 1 }),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Edge {
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
    #[allow(dead_code)]
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
    pub fn from_corrected_reads(reads: &[super::CorrectedRead], k: usize) -> Self {
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
    fn calc_thr_edge(&self) -> u64 {
        let counts = self
            .nodes
            .iter()
            .map(|n| n.edges.iter().fold(0, |x, w| x + w.weight))
            .sum::<u64>();
        let len: usize = self.nodes.iter().map(|n| n.edges.len()).sum::<usize>();
        counts / len as u64 / 3
    }
    pub fn clean_up_auto(self) -> Self {
        let thr = self.calc_thr_edge();
        debug!("Removing edges with weight less than {}", thr);
        self.clean_up(thr)
    }
    fn clean_up(mut self, thr: u64) -> Self {
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
    pub fn assign_corrected_read(&self, read: &super::CorrectedRead) -> Option<usize> {
        let mut count = HashMap::<_, u32>::new();
        for w in read.nodes.windows(self.k) {
            let node = Node::from_corrected(w);
            if let Some(&idx) = self.indexer.get(&node) {
                *count.entry(self.nodes[idx].cluster).or_default() += 1;
            }
        }
        count.into_iter().max_by_key(|x| x.1).map(|x| x.0)
    }
    pub fn coloring(&mut self, _c: &super::GlobalClusteringConfig) {
        // Coloring node of the de Bruijn graph.
        // As a first try, I just color de Bruijn graph by its connected components.
        let mut fu = super::FindUnion::new(self.nodes.len());
        for (from, node) in self.nodes.iter().enumerate().filter(|x| x.1.occ > 0) {
            for edge in node.edges.iter().filter(|e| e.weight > 0) {
                fu.unite(from, edge.to);
            }
        }
        let mut current_component = 1;
        debug!("ClusterID\tNumberOfKmer");
        let mut ignored = 0;
        for cluster in 0..self.nodes.len() {
            if fu.find(cluster).unwrap() != cluster {
                continue;
            }
            let count = (0..self.nodes.len())
                .filter(|&x| fu.find(x).unwrap() == cluster)
                .count();
            if count < 10 {
                ignored += 1;
                continue;
            }
            debug!("{}\t{}", current_component, count);
            for (idx, node) in self.nodes.iter_mut().enumerate() {
                if fu.find(idx).unwrap() == cluster {
                    node.cluster = current_component;
                }
            }
            current_component += 1;
        }
        debug!("Ignored small component (<10 kmer):{}", ignored);
    }
    #[allow(dead_code)]
    pub fn clustering(&self, thr: usize) -> Vec<HashSet<(u64, u64)>> {
        // Clustering de Bruijn graph.
        // As a first try, I implement very naive conneceted component analysis.
        // To this end, I use naive FindUnion Tree. In other words,
        // as I traverse nodes, I merge the two connected nodes.
        // Currently, I ignore very weak connection, i.e., connection
        // with the weight of less than 1.
        let mut fu = super::FindUnion::new(self.nodes.len());
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
    #[allow(dead_code)]
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
    // #[allow(dead_code)]
    // fn path_clustering_test_mcl() {
    //     use rand::SeedableRng;
    //     use rand_xoshiro::Xoshiro256StarStar;
    //     let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
    //     let conf = TestConfig {
    //         cl: 2,
    //         num: 100,
    //         fail: 0.05,
    //         skip: 0.05,
    //         max_len: 15,
    //         min_len: 8,
    //         unit_len: 50,
    //     };
    //     let (reads, answer) = gen_dataset(&mut rng, conf);
    //     let graph = DeBruijnGraph::new(&reads, 3);
    //     let components = graph.clustering_mcl(10, 2);
    //     assert_eq!(components.len(), 2, "{:?}", components,);
    //     {
    //         for c in components.iter() {
    //             let mut c: Vec<_> = c.iter().copied().collect();
    //             c.sort();
    //             eprintln!("{:?}", c);
    //         }
    //     }

    //     let preds: Vec<_> = reads
    //         .iter()
    //         .filter_map(|read| {
    //             components
    //                 .iter()
    //                 .map(|c| read.iter().filter(|x| c.contains(x)).count())
    //                 .enumerate()
    //                 .max_by_key(|x| x.1)
    //         })
    //         .map(|x| x.0)
    //         .collect();
    //     assert_eq!(preds.len(), answer.len());
    //     let correct = answer
    //         .iter()
    //         .zip(preds.iter())
    //         .filter(|&(ans, pred)| {
    //             eprintln!("{}\t{}", ans, pred);
    //             ans == pred
    //         })
    //         .count();
    //     eprintln!("{}/{}", correct, reads.len());
    //     let correct = correct.max(reads.len() - correct);
    //     assert!(correct > reads.len() * 8 / 10);
    // }
    // #[allow(dead_code)]
    // fn path_clustering_test_mcl_with() {
    //     use rand::SeedableRng;
    //     use rand_xoshiro::Xoshiro256StarStar;
    //     let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(3205);
    //     let conf = TestConfig {
    //         cl: 2,
    //         num: 200,
    //         fail: 0.05,
    //         skip: 0.05,
    //         max_len: 20,
    //         min_len: 10,
    //         unit_len: 50,
    //     };
    //     let (reads, answer) = gen_dataset(&mut rng, conf);
    //     let graph = DeBruijnGraph::new(&reads, 3);
    //     let components = graph.clustering_mcl(10, 2);
    //     assert_eq!(components.len(), 2, "{:?}", components,);
    //     {
    //         for c in components.iter() {
    //             let mut c: Vec<_> = c.iter().copied().collect();
    //             c.sort();
    //             eprintln!("{:?}", c);
    //         }
    //     }
    //     let init: Vec<_> = reads
    //         .iter()
    //         .filter_map(|read| {
    //             let counts: Vec<_> = components
    //                 .iter()
    //                 .map(|c| read.iter().filter(|x| c.contains(x)).count())
    //                 .collect();
    //             counts.into_iter().enumerate().max_by_key(|x| x.1)
    //         })
    //         .map(|x| x.0)
    //         .collect();
    //     assert_eq!(init.len(), answer.len());
    //     let preds = {
    //         let units: HashSet<_> = reads.iter().flat_map(|e| e).copied().collect();
    //         let map: HashMap<(u64, u64), usize> =
    //             units.into_iter().enumerate().map(|(x, y)| (y, x)).collect();
    //         let reads: Vec<Vec<_>> = reads
    //             .iter()
    //             .map(|read| read.iter().map(|u| map[u]).collect())
    //             .collect();
    //         path_clustering::path_clustering(&reads, &init)
    //     };
    //     eprintln!("------------------");
    //     {
    //         let max = *preds.iter().max().unwrap();
    //         for cl in 0..=max {
    //             let mut component: Vec<_> = reads
    //                 .iter()
    //                 .zip(preds.iter())
    //                 .filter(|&(_, &x)| x == cl)
    //                 .flat_map(|(p, _)| p.iter())
    //                 .collect::<HashSet<_>>()
    //                 .into_iter()
    //                 .collect();
    //             component.sort();
    //             eprintln!("{:?}", component);
    //         }
    //     }
    //     let correct = answer
    //         .iter()
    //         .zip(preds.iter())
    //         .filter(|&(ans, pred)| ans == pred)
    //         .count();
    //     eprintln!("{}/{}", correct, reads.len());
    //     for ((read, assign), ans) in reads.iter().zip(preds.iter()).zip(answer.iter()) {
    //         eprintln!("{}\t{}\t{}", ans, assign, read.len());
    //     }
    //     let correct = correct.max(reads.len() - correct);
    //     assert!(correct > reads.len() * 8 / 10);
    // }
}
