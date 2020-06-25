use super::find_union::FindUnion;
use definitions;
use std::collections::HashMap;
use std::collections::HashSet;
#[derive(Debug, Clone, Copy)]
pub struct GlobalClusteringConfig {
    pub threads: usize,
    pub k_mer: usize,
}
impl GlobalClusteringConfig {
    pub fn new(threads: usize, k_mer: usize) -> Self {
        Self { threads, k_mer }
    }
}
pub trait GlobalClustering {
    fn global_clustering(self, c: &GlobalClusteringConfig) -> Self;
}

#[derive(Clone)]
struct DeBruijnGraph {
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
        write!(f, "{}\t{}\t[{}]", self.occ, edges.join(","), kmer.join(","))
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
        let (edges, occ) = (vec![], 0);
        Self { kmer, edges, occ }
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
    fn from_encoded_reads(reads: &[definitions::EncodedRead], k: usize) -> Self {
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
    fn calc_thr(&self) -> usize {
        let counts: Vec<_> = self.nodes.iter().map(|n| n.occ).collect();
        let mean = counts.iter().sum::<usize>() / counts.len();
        let thr = mean / 2;
        return thr;
    }
    fn clean_up_auto(self) -> Self {
        let thr = self.calc_thr();
        debug!("Removing nodes with occ less than {}", thr);
        self.clean_up(thr)
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
            n.edges = n
                .edges
                .iter()
                .filter_map(|edge| {
                    map[edge.to].map(|t| Edge {
                        to: t,
                        weight: edge.weight,
                    })
                })
                .collect();
        });
        Self {
            k,
            nodes: resulting_node,
            indexer,
        }
    }
    fn clustering(&self, thr: usize) -> Vec<HashSet<(u64, u64)>> {
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

impl GlobalClustering for definitions::DataSet {
    fn global_clustering(mut self, c: &GlobalClusteringConfig) -> Self {
        let graph = DeBruijnGraph::from_encoded_reads(&self.encoded_reads, c.k_mer);
        if log_enabled!(log::Level::Debug) {
            let length: Vec<_> = self.encoded_reads.iter().map(|r| r.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&length);
            eprintln!("Read({})\n{}", length.len(), hist.format(20, 40));
        }
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("Node({}-mer) occurences\n{}", c.k_mer, hist.format(20, 40));
        }
        let graph = graph.clean_up_auto();
        if log_enabled!(log::Level::Debug) {
            let count: Vec<_> = graph.nodes.iter().map(|n| n.occ).collect();
            let hist = histgram_viz::Histgram::new(&count);
            eprintln!("Node({}-mer) occurences\n{}", c.k_mer, hist.format(20, 40));
        }
        let components: Vec<HashSet<(u64, u64)>> = graph.clustering(0);
        let components = merge(components, &self.encoded_reads);
        debug!("Resulting in {} clusters.", components.len());
        let mut count: HashMap<_, usize> = HashMap::new();
        let assignments: Vec<_> = self
            .encoded_reads
            .iter()
            .map(|read| {
                let id = read.id;
                let cluster = components
                    .iter()
                    .map(|c| sim(c, read))
                    .enumerate()
                    .max_by_key(|x| x.1)
                    .unwrap()
                    .0;
                *count.entry(cluster).or_default() += 1;
                definitions::Assignment { id, cluster }
            })
            .collect();
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
    loop {
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
                    continue;
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
            let (edges, occ) = (vec![], 0);
            Self { kmer, edges, occ }
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
        let components = graph.clustering(1);
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
            fail: 0.01,
            skip: 0.01,
            max_len: 20,
            min_len: 10,
            unit_len: 800,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let components = graph.clustering(1);
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
            fail: 0.05,
            skip: 0.05,
            max_len: 20,
            min_len: 10,
            unit_len: 50,
        };
        let (reads, answer) = gen_dataset(&mut rng, conf);
        let graph = DeBruijnGraph::new(&reads, 3);
        let graph = graph.clean_up_auto();
        eprintln!("{:?}", graph);
        let components = graph.clustering(0);
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
}
