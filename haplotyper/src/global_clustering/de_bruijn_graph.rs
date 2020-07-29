#![allow(dead_code)]
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct DeBruijnGraph {
    k: usize,
    pub nodes: Vec<Node>,
    indexer: HashMap<Node, usize>,
}

impl std::fmt::Debug for DeBruijnGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "NumOfNodes:{}", self.nodes.len())?;
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
    fn remove_edge(&mut self, to: usize) {
        self.edges.retain(|x| x.to != to);
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
    pub fn from_encoded_reads(reads: &[definitions::EncodedRead], k: usize) -> Self {
        let (mut nodes, mut indexer) = (vec![], HashMap::new());
        for read in reads {
            for (idx, w) in read.nodes.windows(k + 1).enumerate() {
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
                nodes[from].push(to);
                nodes[to].push(from);
                if idx == 0 {
                    nodes[from].occ += 1;
                }
                nodes[to].occ += 1;
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
    pub fn assign_encoded_read(&self, read: &definitions::EncodedRead) -> Option<usize> {
        let mut count = HashMap::<_, u32>::new();
        for w in read.nodes.windows(self.k) {
            let node = Node::new(w);
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
    pub fn resolve_crossings(&mut self, _reads: &[definitions::EncodedRead]) {}
    pub fn resolve_bubbles(&mut self, reads: &[definitions::EncodedRead]) {
        let mut bubble_spec = self.enumerate_bubbles(reads);
        let mut queue_and_parent: std::collections::VecDeque<_> = bubble_spec
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| e.map(|_| (idx, idx)))
            .collect();
        let mut to_be_removed = vec![false; self.nodes.len()];
        while let Some((idx, parent)) = queue_and_parent.pop_back() {
            let bubble = match bubble_spec[idx] {
                Some(x) => x,
                None => continue,
            };
            // search nearest bubble.
            eprintln!("{:?}", bubble);
            let pair_pos = match self.search_nearest_bubble(bubble.root, bubble.shoot, &bubble_spec)
            {
                Ok(res) => res,
                Err(ReachedNodeType::Root) => {
                    // We cut the bubble arbitrary.
                    let branch = bubble.branches.0;
                    self.nodes[branch].remove_edge(bubble.shoot);
                    self.nodes[bubble.shoot].remove_edge(branch);
                    continue;
                }
                Err(ReachedNodeType::BubbleBranch(p)) if parent == p => {
                    // We reached the same parent again.
                    // I think this bubble could not be resolved. Just abandon.
                    continue;
                }
                Err(ReachedNodeType::BubbleBranch(p)) => {
                    // We reached a branch of a bubble. However, thre is
                    // some progress.
                    // Hope we can resolve this bubble in the next time.
                    queue_and_parent.push_front((idx, p));
                    continue;
                }
                Err(ReachedNodeType::ComplexNode) => continue,
            };
            // Never panic.
            let pair = bubble_spec[pair_pos].unwrap();
            // Resolve node.
            bubble_spec[idx] = None;
            bubble_spec[pair_pos] = None;
            let resolved_pairs = self.resolve_bubble(reads, bubble, pair);
            // First, remove all the edges from branches to the root.
            for bbl in &[pair, bubble] {
                let ((b0, b1), shoot) = (bbl.branches, bbl.shoot);
                self.nodes[b0].remove_edge(shoot);
                self.nodes[b1].remove_edge(shoot);
                self.nodes[shoot].remove_edge(b0);
                self.nodes[shoot].remove_edge(b1);
            }
            // Then, connect anchored pairs.
            for &(n0, n1) in resolved_pairs.iter() {
                self.nodes[n0].push(n1);
                self.nodes[n1].push(n0);
            }
            // Then, remove bubbles if the node is branch-free.
            for &(n0, n1) in resolved_pairs.iter() {
                if self.nodes[n0].edges.len() <= 2 {
                    bubble_spec[n0] = None;
                }
                if self.nodes[n1].edges.len() <= 2 {
                    bubble_spec[n1] = None;
                }
            }
            // Then, check bubbles at anchored positions,
            for &(n0, n1) in resolved_pairs.iter() {
                if let Some(b) = bubble_spec.get_mut(n0).unwrap().as_mut() {
                    if b.branches.0 == bubble.shoot || b.branches.0 == pair.shoot {
                        b.branches.0 = n1;
                    } else if b.branches.1 == bubble.shoot || b.branches.1 == pair.shoot {
                        b.branches.1 = n1;
                    } else if b.root == bubble.shoot || b.root == pair.shoot {
                        b.root = n1;
                    }
                }
                if let Some(b) = bubble_spec.get_mut(n1).unwrap().as_mut() {
                    if b.branches.0 == bubble.shoot || b.branches.0 == pair.shoot {
                        b.branches.0 = n0;
                    } else if b.branches.1 == bubble.shoot || b.branches.1 == pair.shoot {
                        b.branches.1 = n0;
                    } else if b.root == bubble.shoot || b.root == pair.shoot {
                        b.root = n0;
                    }
                }
            }
            // Lastly, register all the nodes between two shoots as `removed`
            let (mut prev, mut current) = (bubble.shoot, bubble.root);
            while current != pair.shoot {
                to_be_removed[current] = true;
                let next_nodes = &self.nodes[current].edges;
                assert!(next_nodes.len() == 2);
                let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                prev = current;
                current = next;
            }
            to_be_removed[bubble.shoot] = true;
            to_be_removed[pair.shoot] = true;
        }
        // Removing nodes.
        self.remove_nodes(&to_be_removed);
    }
    fn resolve_bubble(
        &self,
        reads: &[definitions::EncodedRead],
        bubble0: Bubble,
        bubble1: Bubble,
    ) -> Vec<(usize, usize)> {
        eprintln!("Resolving {:?} {:?}", bubble0, bubble1);
        // [00->10, 01->10, 00->11, 01->11];
        let mut connection_counts = [0; 4];
        // TODO:Need faster algorithm here.
        for read in reads.iter() {
            let (mut b00, mut b01) = (false, false);
            let (mut b10, mut b11) = (false, false);
            for &node in read
                .nodes
                .windows(self.k)
                .filter_map(|w| self.indexer.get(&Node::new(&w)))
            {
                b00 |= node == bubble0.branches.0;
                b01 |= node == bubble0.branches.1;
                b10 |= node == bubble1.branches.0;
                b11 |= node == bubble1.branches.1;
            }
            assert!(!b00 || !b01);
            assert!(!b10 || !b11);
            if (b00 || b01) && (b10 || b11) {
                let b0 = if b00 { 0 } else { 1 };
                let b1 = if b10 { 0 } else { 1 };
                connection_counts[(b1 << 1) + b0] += 1;
            }
        }
        // TODO: parametrize here.
        connection_counts.iter_mut().for_each(|x| {
            if *x <= 1 {
                *x = 0;
            }
        });
        // Case0: bubble0.0 -> bubble1.0 and bubble0.1 -> bubble1.1
        let case0 = connection_counts[0] + connection_counts[3];
        // Case1: bubble0.0 -> bubble1.1 and bubble0.0 -> bubble1.0
        let case1 = connection_counts[1] + connection_counts[2];
        let mut pairs = vec![];
        // TODO: parametrize here.
        if case0 > case1 {
            if connection_counts[0] > 0 {
                pairs.push((bubble0.branches.0, bubble1.branches.0));
            }
            if connection_counts[3] > 0 {
                pairs.push((bubble0.branches.1, bubble1.branches.1));
            }
        } else {
            if connection_counts[1] > 0 {
                pairs.push((bubble0.branches.1, bubble1.branches.0));
            }
            if connection_counts[2] > 0 {
                pairs.push((bubble0.branches.0, bubble1.branches.1));
            }
        }
        pairs
    }
    fn enumerate_bubbles(&self, reads: &[definitions::EncodedRead]) -> Vec<Option<Bubble>> {
        // Enumerate bubble.
        let mut edge_counts: Vec<_> = (0..self.nodes.len()).map(|idx| vec![0; idx]).collect();
        for read in reads.iter() {
            for w in read.nodes.windows(self.k + 2) {
                let from = match self.indexer.get(&Node::new(&w[..self.k])) {
                    Some(&res) => res,
                    None => continue,
                };
                let to = match self.indexer.get(&Node::new(&w[2..])) {
                    Some(&res) => res,
                    None => continue,
                };
                edge_counts[from.max(to)][from.min(to)] += 1;
            }
        }
        self.nodes
            .iter()
            .enumerate()
            .map(|(shoot, node)| {
                if node.edges.len() != 3 {
                    return None;
                }
                let (to0, to1, to2) = (node.edges[0].to, node.edges[1].to, node.edges[2].to);
                let bet0and1 = edge_counts[to0.max(to1)][to0.min(to1)];
                let bet0and2 = edge_counts[to0.max(to2)][to0.min(to2)];
                let bet1and2 = edge_counts[to1.max(to2)][to1.min(to2)];
                let (branches, root) = if bet0and1 == 0 && bet0and2 > 0 && bet1and2 > 0 {
                    ((to0, to1), to2)
                } else if bet0and1 > 0 && bet0and2 == 0 && bet1and2 > 0 {
                    ((to0, to2), to1)
                } else if bet0and1 > 0 && bet0and2 > 0 && bet1and2 == 0 {
                    ((to1, to2), to0)
                } else {
                    return None;
                };
                Some(Bubble {
                    branches,
                    shoot,
                    root,
                })
            })
            .collect()
    }
    fn remove_nodes(&mut self, to_be_removed: &[bool]) {
        let mut next_index = vec![];
        {
            let mut index = 0;
            for &b in to_be_removed.iter() {
                next_index.push(index);
                index += !b as usize;
            }
        }
        let mut index = 0;
        self.nodes.retain(|_| {
            index += 1;
            !to_be_removed[index - 1]
        });
        self.nodes.iter_mut().for_each(|n| {
            n.edges.iter_mut().for_each(|x| x.to = next_index[x.to]);
        });
        self.indexer
            .iter_mut()
            .for_each(|(_, x)| *x = next_index[*x]);
    }
    fn search_nearest_bubble(
        &self,
        root: usize,
        shoot: usize,
        bubbles: &[Option<Bubble>],
    ) -> Result<usize, ReachedNodeType> {
        let (mut prev, mut current) = (shoot, root);
        while bubbles[current].is_none() {
            let next_nodes = &self.nodes[current].edges;
            if next_nodes.len() == 1 {
                // Fail to find pair. But this is bubble is terminating.
                return Err(ReachedNodeType::Root);
            } else if next_nodes.len() == 2 {
                // Proceed.
                let next = next_nodes.iter().find(|e| e.to != prev).unwrap().to;
                prev = current;
                current = next;
            } else {
                // Fail to find pair. We've reached very complex bubble.
                return Err(ReachedNodeType::ComplexNode);
            }
        }
        if bubbles[current].unwrap().root == prev {
            // We entered current node from the root node. Well done!
            Ok(current)
        } else {
            // We entered current node from either of branches.
            let (b0, b1) = bubbles[current].unwrap().branches;
            assert!(prev == b0 || prev == b1);
            Err(ReachedNodeType::BubbleBranch(current))
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Bubble {
    // Index of bubble.
    branches: (usize, usize),
    shoot: usize,
    // Root. Where this bubble collapse.
    root: usize,
}

#[derive(Debug, Clone, Copy)]
enum ReachedNodeType {
    Root,
    BubbleBranch(usize),
    ComplexNode,
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
    fn to_encoded_reads(reads: &[Vec<(u64, u64)>]) -> Vec<definitions::EncodedRead> {
        reads
            .iter()
            .enumerate()
            .map(|(idx, units)| {
                let edges = vec![];
                let nodes: Vec<_> = units
                    .iter()
                    .map(|&(unit, cluster)| definitions::Node {
                        position_from_start: 0,
                        unit: unit,
                        cluster: cluster,
                        seq: String::new(),
                        is_forward: true,
                        cigar: vec![],
                    })
                    .collect();
                definitions::EncodedRead {
                    original_length: 0,
                    leading_gap: 0,
                    trailing_gap: 0,
                    id: idx as u64,
                    edges: edges,
                    nodes: nodes,
                }
            })
            .collect()
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
    // Bubble resolving functionality.
    #[test]
    fn bubble_works() {
        let reads = vec![vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        graph.resolve_bubbles(&reads)
    }
    fn sample_from_to(
        template: &[(u64, u64)],
        from: usize,
        to: usize,
        rev: bool,
    ) -> Vec<(u64, u64)> {
        if rev {
            template[from..to].iter().rev().copied().collect()
        } else {
            template[from..to].to_vec()
        }
    }
    #[test]
    fn bubble_case_0() {
        let templates: Vec<_> = vec![vec![0; 8], vec![1, 1, 0, 0, 0, 0, 1, 0]];
        let templates: Vec<Vec<(u64, u64)>> = templates
            .iter()
            .map(|ts| {
                ts.iter()
                    .enumerate()
                    .map(|(idx, &x)| (idx as u64, x as u64))
                    .collect()
            })
            .collect();
        let reads = vec![
            sample_from_to(&templates[0], 0, 8, false),
            sample_from_to(&templates[1], 0, 8, false),
            sample_from_to(&templates[0], 0, 8, true),
            sample_from_to(&templates[1], 0, 8, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 2, "{:?}", graph);
        let pair_idx = graph.search_nearest_bubble(3, 2, &bubbles).unwrap();
        let bubble = bubbles[pair_idx].unwrap();
        assert_eq!(bubble.shoot, 3);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
    }
    #[test]
    fn bubble_case_1() {
        let templates: Vec<_> = vec![
            vec![0; 9],
            vec![0, 1, 0, 0, 0, 0, 1, 1, 0],
            vec![0, 1, 0, 0, 0, 0, 1, 1, 2],
        ];
        let templates: Vec<Vec<(u64, u64)>> = templates
            .iter()
            .map(|ts| {
                ts.iter()
                    .enumerate()
                    .map(|(idx, &x)| (idx as u64, x as u64))
                    .collect()
            })
            .collect();
        let reads = vec![
            sample_from_to(&templates[0], 0, 9, false),
            sample_from_to(&templates[1], 0, 9, false),
            sample_from_to(&templates[2], 0, 9, false),
            sample_from_to(&templates[0], 0, 9, true),
            sample_from_to(&templates[1], 0, 9, true),
            sample_from_to(&templates[2], 0, 9, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 3, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
    }
    #[test]
    fn bubble_case_2() {
        let template_a: Vec<(u64, u64)> = (0..11).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![1, 1, 0, 0, 1, 1, 2, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64 + 3, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 11, false),
            sample_from_to(&template_a, 0, 11, true),
            sample_from_to(&template_b, 0, 11, false),
            sample_from_to(&template_b, 0, 11, true),
            sample_from_to(&template_c, 0, 8, false),
            sample_from_to(&template_c, 0, 8, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
    }
    #[test]
    fn bubble_case_3() {
        let template_a: Vec<(u64, u64)> = (0..13).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 13, false),
            sample_from_to(&template_a, 0, 13, true),
            sample_from_to(&template_b, 0, 13, false),
            sample_from_to(&template_b, 0, 13, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
    }
    #[test]
    fn bubble_case_4() {
        let template_a: Vec<(u64, u64)> = (0..12).map(|idx| (idx, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 12, false),
            sample_from_to(&template_a, 0, 12, true),
            sample_from_to(&template_b, 0, 12, false),
            sample_from_to(&template_b, 0, 12, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 3, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 2, "{:?}", graph);
    }
    #[test]
    fn bubble_case_5() {
        let template_a: Vec<(u64, u64)> = (0..9).map(|idx| (idx + 2, 0)).collect();
        let template_b: Vec<(u64, u64)> = vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let template_c: Vec<(u64, u64)> = vec![2, 2, 1, 1, 0, 0, 0, 0, 0, 1, 1, 2, 1]
            .into_iter()
            .enumerate()
            .map(|(idx, c)| (idx as u64, c))
            .collect();
        let reads = vec![
            sample_from_to(&template_a, 0, 9, false),
            sample_from_to(&template_a, 0, 9, true),
            sample_from_to(&template_b, 0, 13, false),
            sample_from_to(&template_b, 0, 13, true),
            sample_from_to(&template_c, 0, 13, false),
            sample_from_to(&template_c, 0, 13, true),
        ];
        let reads = to_encoded_reads(&reads);
        let mut graph = DeBruijnGraph::from_encoded_reads(&reads, 3);
        eprintln!("{:?}", graph);
        let bubbles = graph.enumerate_bubbles(&reads);
        let num_bubbles = bubbles.iter().filter(|e| e.is_some()).count();
        assert_eq!(num_bubbles, 4, "{:?}", graph);
        graph.resolve_bubbles(&reads);
        assert_eq!(graph.clustering(0).len(), 3, "{:?}", graph);
    }
}
