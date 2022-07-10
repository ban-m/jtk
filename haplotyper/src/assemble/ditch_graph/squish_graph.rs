//! Squishing graph and other modules...
//!
//!
use super::DitchEdge;
use super::DitchGraph;
use super::NodeIndex;
use std::collections::{HashMap, HashSet};
type TempNode = (usize, [Vec<(usize, usize)>; 2]);
use rayon::prelude::*;
impl<'a> DitchGraph<'a> {
    /// Squish small net-like-structure like below:
    /// [Long contig]---[Small contig]---[Long contig]
    ///               X                X
    /// [Long contig]---[Small contig]---[Long contig]
    /// By squishing, only one side of the [small contig] would be retained.
    pub fn squish_small_net(&mut self, len: usize) {
        let (node_to_pathid, connecting_edges) = self.reduce_simple_path();
        let nodes = Self::temp_graph(&node_to_pathid, &connecting_edges);
        // vector of found gruops
        let mut suspicious_nodes: Vec<Vec<usize>> = nodes
            .par_iter()
            .enumerate()
            .filter(|(_, n)| len < n.0)
            .flat_map(|(i, _)| {
                let mut suspic = Vec::new();
                if let Some(simple_paths) = Self::is_net(&nodes, i, 0, len) {
                    suspic.push(simple_paths);
                }
                if let Some(simple_paths) = Self::is_net(&nodes, i, 1, len) {
                    suspic.push(simple_paths);
                }
                suspic
            })
            .collect();
        // Dedup.
        suspicious_nodes.sort();
        suspicious_nodes.dedup();
        let to_remove: HashSet<_> = suspicious_nodes
            .iter()
            .flat_map(|nodes| nodes[1..].iter())
            .collect();
        let to_remove: HashSet<_> = self
            .nodes()
            .filter(|&(nodeidx, _)| to_remove.contains(&node_to_pathid[&nodeidx]))
            .map(|(i, _)| i)
            .collect();
        // Removing this nodes...
        self.remove_nodes(&to_remove);
    }
    fn is_net(nodes: &[TempNode], from: usize, from_slot: usize, len: usize) -> Option<Vec<usize>> {
        // check if this is branching...
        if nodes[from].1[from_slot].len() <= 1 {
            return None;
        }
        // Check if connedted nodes are all short.
        if nodes[from].1[from_slot]
            .iter()
            .any(|&(to, _)| len < nodes[to].0)
        {
            return None;
        }
        let (first_child_node, first_child_pos) = nodes[from].1[from_slot][0];
        let sibs = &nodes[first_child_node].1[first_child_pos];
        // At least there are one additional sibling.
        assert!(sibs.contains(&(from, from_slot)));
        if sibs.len() <= 1 {
            return None;
        }
        // Check if all the children reflects the same siblings.
        if nodes[from].1[from_slot]
            .iter()
            .any(|&(node, pos)| &nodes[node].1[pos] != sibs)
        {
            return None;
        }
        let rev_pos = (first_child_pos == 0) as usize;
        assert_ne!(rev_pos, first_child_pos);
        let destination_after_path = &nodes[first_child_node].1[rev_pos];
        // Check if all the destinations are long.
        if destination_after_path
            .iter()
            .any(|&(dst, _)| nodes[dst].0 <= len)
        {
            return None;
        }
        // Check if all the children converged into the same destination.
        if nodes[from].1[from_slot].iter().any(|&(node, pos)| {
            let rev_pos = (pos == 0) as usize;
            &nodes[node].1[rev_pos] != destination_after_path
        }) {
            return None;
        }
        let nodes: Vec<_> = nodes[from].1[from_slot].iter().map(|x| x.0).collect();
        Some(nodes)
    }
    fn temp_graph(
        nodes_to_pathid: &HashMap<NodeIndex, usize>,
        edges: &[&DitchEdge],
    ) -> Vec<TempNode> {
        let path_num = *nodes_to_pathid.values().max().unwrap() + 1;
        let mut nodes = vec![(0, [Vec::new(), Vec::new()]); path_num];
        for &path_id in nodes_to_pathid.values() {
            nodes[path_id].0 += 1;
        }
        let mut terminals = vec![Vec::with_capacity(2); path_num];
        for edge in edges.iter() {
            let from = nodes_to_pathid[&edge.from];
            let from_slot = terminals[from]
                .iter()
                .position(|&e| e == (edge.from, edge.from_position));
            let from_slot = match from_slot {
                Some(idx) => idx,
                None => {
                    assert!(terminals[from].len() < 2);
                    terminals[from].push((edge.from, edge.from_position));
                    terminals[from].len() - 1
                }
            };
            let to = nodes_to_pathid[&edge.to];
            let to_slot = terminals[to]
                .iter()
                .position(|&e| e == (edge.to, edge.to_position));
            let to_slot = match to_slot {
                Some(idx) => idx,
                None => {
                    assert!(terminals[to].len() < 2);
                    terminals[to].push((edge.to, edge.to_position));
                    terminals[to].len() - 1
                }
            };
            nodes[from].1[from_slot].push((to, to_slot));
            nodes[to].1[to_slot].push((from, from_slot));
        }
        nodes.iter_mut().for_each(|x| {
            x.1[0].sort_unstable();
            x.1[1].sort_unstable();
        });
        nodes
    }
}
