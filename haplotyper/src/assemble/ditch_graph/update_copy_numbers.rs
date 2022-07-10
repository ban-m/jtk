//! Methods to update copy numbers.
//!
use super::copy_num_by_mst;
use super::DitEdge;
use super::DitchEdge;
use super::DitchGraph;
use super::EdgeBetweenSimplePath;
use super::NodeIndex;
use super::Position;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct CopyNumbers {
    inner: HashMap<NodeIndex, usize>,
}

impl std::convert::From<HashMap<NodeIndex, usize>> for CopyNumbers {
    fn from(inner: HashMap<NodeIndex, usize>) -> Self {
        Self { inner }
    }
}

impl CopyNumbers {
    pub fn new(inner: HashMap<NodeIndex, usize>) -> Self {
        Self { inner }
    }
}

impl std::ops::Index<NodeIndex> for CopyNumbers {
    type Output = usize;
    fn index(&self, index: NodeIndex) -> &Self::Output {
        &self.inner[&index]
    }
}

impl<'a> DitchGraph<'a> {
    /// Estimoate copy number of nodes and edges.
    /// *This function does not modify the graph content*. If you want
    /// to assign copy number to each node, call `assign_copy_number` instead.
    pub fn copy_number_estimation(
        &self,
        cov: f64,
        _lens: &[usize],
    ) -> (CopyNumbers, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (terminals, mut edges) =
            self.convert_connecting_edges(&node_to_pathid, &connecting_edges);
        let mut nodes = self.convert_path_weight(&node_to_pathid);
        edges.iter_mut().for_each(|x| x.4 /= cov);
        nodes.iter_mut().for_each(|x| x.0 /= cov);
        let (node_cp, edge_cp) = crate::assemble::copy_number::estimate_copy_number(&nodes, &edges);
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }

    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) = self.copy_number_estimation(naive_cov, lens);
        self.modify_by_with_index(|(index, node)| {
            let cp = node_copy_number[index];
            node.copy_number = Some(cp);
            for edge in node.edges.iter_mut() {
                edge.copy_number = edge_copy_number.get(&edge.key()).copied();
            }
        });
    }
    /// Estimoate copy number of nodes and edges by a gibbs sampler.
    /// *This function does not modify the graph content*.
    /// If you want to assign copy number to each node, call `assign_copy_number_gbs` instead.
    pub fn copy_number_estimation_gbs(
        &self,
        cov: f64,
        _lens: &[usize],
    ) -> (CopyNumbers, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (terminals, edges) = self.convert_connecting_edges(&node_to_pathid, &connecting_edges);
        let nodes = self.convert_path_weight(&node_to_pathid);
        let nodes: Vec<_> = nodes.iter().map(|x| x.0).collect();
        let (node_cp, edge_cp) =
            crate::assemble::copy_number::estimate_copy_number_gbs(&nodes, &edges, cov);
        if log_enabled!(log::Level::Trace) {
            trace!("COVCP\tType\tCov\tCp");
            for (n, cp) in nodes.iter().zip(node_cp.iter()) {
                trace!("COVCP\tNODE\t{n:.2}\t{cp}");
            }
            for (e, cp) in edges.iter().zip(edge_cp.iter()) {
                trace!("COVCP\tEDGE\t{:.2}\t{cp}", e.4);
            }
        }
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number_gbs(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) = self.copy_number_estimation_gbs(naive_cov, lens);
        self.modify_by_with_index(|(index, node)| {
            let cp = node_copy_number[index];
            node.copy_number = Some(cp);
            for edge in node.edges.iter_mut() {
                edge.copy_number = edge_copy_number.get(&edge.key()).copied();
            }
        });
    }
    /// Estimoate copy number of nodes and edges by MCMC.
    /// *This function does not modify the graph content*.
    /// If you want to assign copy number to each node, call `assign_copy_number_gbs` instead.
    pub fn copy_number_estimation_mcmc(
        &self,
        cov: f64,
        _lens: &[usize],
    ) -> (CopyNumbers, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (terminals, edges) = self.convert_connecting_edges(&node_to_pathid, &connecting_edges);
        let nodes = self.convert_path_weight(&node_to_pathid);
        let (node_cp, edge_cp) =
            crate::assemble::copy_number::estimate_copy_number_mcmc(&nodes, &edges, cov);
        if log_enabled!(log::Level::Trace) {
            trace!("COVCP\tType\tCov\tCp");
            for ((cov, len), cp) in nodes.iter().zip(node_cp.iter()) {
                trace!("COVCP\tNODE\t{}\t{}\t{}", cov, len, cp,);
            }
            for ((f, fp, t, tp, cov), cp) in edges.iter().zip(edge_cp.iter()) {
                trace!("COVCP\tEDGE\t{f}\t{fp}\t{t}\t{tp}\t{cov}\t{cp}");
            }
        }
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }
    // TODO:Fasten this function.
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number_mcmc(&mut self, naive_cov: f64, lens: &[usize]) {
        let (node_copy_number, edge_copy_number) =
            self.copy_number_estimation_mcmc(naive_cov, lens);
        self.modify_by_with_index(|(index, node)| {
            let cp = node_copy_number[index];
            node.copy_number = Some(cp);
            for edge in node.edges.iter_mut() {
                edge.copy_number = edge_copy_number.get(&edge.key()).copied();
            }
        });
    }
    /// Estimoate copy number of nodes and edges by MCMC.
    /// *This function does not modify the graph content*.
    /// If you want to assign copy number to each node, call `assign_copy_number_gbs` instead.
    pub fn copy_number_estimation_mst<R: rand::Rng>(
        &self,
        hap_cov: f64,
        rng: &mut R,
    ) -> (CopyNumbers, HashMap<DitEdge, usize>) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let (terminals, edges) = self.convert_connecting_edges(&node_to_pathid, &connecting_edges);
        let nodes = self.convert_path_weight(&node_to_pathid);
        let mut fat_edges = vec![];
        let mut fat_self_loops = vec![];
        for (i, &(target, len)) in nodes.iter().enumerate() {
            let occ = target.ceil() as usize;
            fat_edges.push(copy_num_by_mst::FatEdge::new(2 * i, 2 * i + 1, occ, len));
        }
        for &(from, fdir, to, todir, target) in edges.iter() {
            let (from, to) = (2 * from + fdir as usize, 2 * to + todir as usize);
            let edge = copy_num_by_mst::FatEdge::new(from, to, target.ceil() as usize, 1);
            if from / 2 != to / 2 {
                fat_edges.push(edge);
            } else {
                fat_self_loops.push(edge);
            }
        }
        debug!(
            "MST\t{}\t{}\t{}",
            hap_cov,
            fat_edges.len(),
            fat_self_loops.len()
        );
        let mut graph = copy_num_by_mst::Graph::new(hap_cov, fat_edges, fat_self_loops);
        let config = copy_num_by_mst::MSTConfig::default();
        graph.update_copy_numbers(rng, &config);
        let mut node_cps = vec![];
        let mut edge_cps = HashMap::new();
        for edge in graph.edges() {
            if edge.from / 2 == edge.to / 2 {
                node_cps.push((edge.from / 2, edge.copy_number.max(0) as usize));
            } else {
                edge_cps.insert((edge.from, edge.to), edge.copy_number.max(0) as usize);
            }
        }
        for edge in graph.self_loops() {
            edge_cps.insert((edge.from, edge.to), edge.copy_number.max(0) as usize);
        }
        node_cps.sort_by_key(|x| x.0);
        let node_cp: Vec<_> = node_cps.into_iter().map(|x| x.1).collect();
        let edge_cp: Vec<_> = edges
            .iter()
            .map(|&(from, fdir, to, todir, _)| {
                let from = 2 * from + fdir as usize;
                let to = 2 * to + todir as usize;
                edge_cps[&(from.min(to), from.max(to))]
            })
            .collect();
        assert_eq!(nodes.len(), node_cp.len());
        assert_eq!(edges.len(), edge_cp.len());
        self.gather_answer(&edges, &node_cp, &edge_cp, &node_to_pathid, &terminals)
    }
    /// (Re-)estimate copy number on each node and edge.
    pub fn assign_copy_number_mst<R: rand::Rng>(&mut self, naive_cov: f64, rng: &mut R) {
        let (node_copy_number, edge_copy_number) = self.copy_number_estimation_mst(naive_cov, rng);
        self.modify_by_with_index(|(index, node)| {
            let cp = node_copy_number[index];
            node.copy_number = Some(cp);
            for edge in node.edges.iter_mut() {
                edge.copy_number = edge_copy_number.get(&edge.key()).copied();
            }
        });
    }

    // Partition edges whether or not it is in simple path/not.
    // Only (from, _, to, _) with from <= to edges would be in these vectors.
    fn partition_edges_by_simple_path(&self) -> (Vec<&DitchEdge>, Vec<&DitchEdge>) {
        let degrees: HashMap<NodeIndex, (u32, u32)> =
            self.nodes()
                .map(|(index, node)| {
                    let (head, tail) = node.edges.iter().fold((0, 0), |(head, tail), edge| {
                        match edge.from_position {
                            Position::Head => (head + 1, tail),
                            Position::Tail => (head, tail + 1),
                        }
                    });
                    (index, (head, tail))
                })
                .collect();
        self.nodes()
            .flat_map(|(_, node)| node.edges.iter())
            .filter(|edge| edge.from <= edge.to)
            .partition(|edge| {
                let from_deg = match edge.from_position {
                    Position::Head => degrees[&edge.from].0,
                    Position::Tail => degrees[&edge.from].1,
                };
                let to_deg = match edge.to_position {
                    Position::Head => degrees[&edge.to].0,
                    Position::Tail => degrees[&edge.to].1,
                };
                from_deg == 1 && to_deg == 1
            })
    }

    // Return Node->serialized simple-path id hashmapping and
    // edges between simple paths.
    pub fn reduce_simple_path(&self) -> (HashMap<NodeIndex, usize>, Vec<&DitchEdge>) {
        let (edges_in_simple_path, edges_between_simple_path) =
            self.partition_edges_by_simple_path();
        use crate::find_union::FindUnion;
        let mut fu = FindUnion::new(self.nodes.len());
        for edge in edges_in_simple_path.iter() {
            fu.unite(edge.from.0, edge.to.0);
        }
        let cluster_index: HashMap<_, _> = self
            .nodes()
            .map(|(index, _)| index.0)
            .filter(|&n| fu.find(n) == Some(n))
            .enumerate()
            .map(|(idx, node)| (node, idx))
            .collect();
        {
            let mut plug_num = vec![HashSet::new(); cluster_index.len()];
            for edge in edges_between_simple_path.iter() {
                let from_index = cluster_index[&fu.find(edge.from.0).unwrap()];
                plug_num[from_index].insert((edge.from, edge.from_position));
                let to_index = cluster_index[&fu.find(edge.to.0).unwrap()];
                plug_num[to_index].insert((edge.to, edge.to_position));
            }
            for plugs in plug_num {
                assert!(plugs.len() <= 2, "{:?}", plugs);
            }
        }
        let node_to_cluster: HashMap<_, _> = self
            .nodes()
            .map(|(index, _)| (index, cluster_index[&fu.find(index.0).unwrap()]))
            .collect();
        (node_to_cluster, edges_between_simple_path)
    }
    // Aggregate nodes, aggregates edges.
    fn convert_connecting_edges(
        &self,
        node_to_pathid: &HashMap<NodeIndex, usize>,
        edges: &[&DitchEdge],
    ) -> (Vec<Vec<(NodeIndex, Position)>>, Vec<EdgeBetweenSimplePath>) {
        let path_num: usize = *node_to_pathid.values().max().unwrap_or(&0) + 1;
        let mut terminals: Vec<_> = vec![vec![]; path_num];
        let edges: Vec<_> = edges
            .iter()
            .map(|edge| {
                let from_index = node_to_pathid[&edge.from];
                let to_index = node_to_pathid[&edge.to];
                // Decide the "position" of these nodes.
                // There should not be more than two plugs on each simple path.
                let from_node = (edge.from, edge.from_position);
                let fp = match terminals[from_index].iter().position(|n| n == &from_node) {
                    Some(idx) => idx == 1,
                    None => {
                        terminals[from_index].push(from_node);
                        terminals[from_index].len() == 2
                    }
                };
                assert!(
                    terminals[from_index].len() <= 2,
                    "{:?}",
                    terminals[from_index]
                );
                let to_node = (edge.to, edge.to_position);
                let tp = match terminals[to_index].iter().position(|n| n == &to_node) {
                    Some(idx) => idx == 1,
                    None => {
                        terminals[to_index].push(to_node);
                        terminals[to_index].len() == 2
                    }
                };
                assert!(terminals[to_index].len() <= 2);
                (from_index, fp, to_index, tp, edge.occ as f64)
            })
            .collect();
        (terminals, edges)
    }
    fn convert_path_weight(&self, node_to_pathid: &HashMap<NodeIndex, usize>) -> Vec<(f64, usize)> {
        let path_num: usize = *node_to_pathid.values().max().unwrap_or(&0) + 1;
        let mut node_weights = vec![(0f64, 0usize); path_num];
        for (index, node) in self.nodes() {
            let simple_path_index = node_to_pathid[&index];
            node_weights[simple_path_index].1 += 1;
            node_weights[simple_path_index].0 += node.occ as f64;
        }
        node_weights
            .iter()
            .map(|&(sum, len)| (sum / len as f64, len))
            .collect()
    }
    fn gather_answer(
        &self,
        edges: &[(usize, bool, usize, bool, f64)],
        node_cp: &[usize],
        edge_cp: &[usize],
        node_to_pathid: &HashMap<NodeIndex, usize>,
        terminals: &[Vec<(NodeIndex, Position)>],
    ) -> (CopyNumbers, HashMap<DitEdge, usize>) {
        let node_copy_number: HashMap<_, _> = node_to_pathid
            .iter()
            .map(|(&idx, &pathid)| (idx, node_cp[pathid]))
            .collect();
        let mut edge_copy_number = HashMap::new();
        for edge in self.nodes().flat_map(|(_, node)| node.edges.iter()) {
            let from = node_to_pathid[&edge.from];
            let to = node_to_pathid[&edge.to];
            if from == to {
                edge_copy_number.insert(edge.key(), node_cp[from]);
            }
        }
        for (&edge_cp, &(from, fplus, to, tplus, _)) in edge_cp.iter().zip(edges.iter()) {
            let from = terminals[from]
                .get(fplus as usize)
                .map(|&(idx, pos)| (idx, pos))
                .unwrap();
            let to = terminals[to]
                .get(tplus as usize)
                .map(|&(idx, pos)| (idx, pos))
                .unwrap();
            edge_copy_number.insert((from, to), edge_cp);
            edge_copy_number.insert((to, from), edge_cp);
        }
        let node_copy_number = CopyNumbers::new(node_copy_number);
        (node_copy_number, edge_copy_number)
    }
}
