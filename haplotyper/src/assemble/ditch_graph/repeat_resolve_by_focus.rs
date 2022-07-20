use std::collections::HashMap;

const ERROR_PROB: f64 = 0.05;
use super::super::AssembleConfig;
use super::DitchGraph;
use super::DitchNode;
use super::Node;
use super::NodeIndex;
use super::Position;
use definitions::*;
/// A focus from a node of a ditch graph.
/// Here, we tag a node with its position, i.e., Position.
#[derive(Debug, Clone)]
pub struct Focus {
    pub from: NodeIndex,
    pub from_node: Node,
    pub from_position: Position,
    pub to: NodeIndex,
    pub to_node: Node,
    pub to_position: Position,
    pub dist: usize,
    /// Log Likelihood ratio between the alt hypothesis/null hypothesis.
    pub log_likelihood_ratio: f64,
    pub counts: Vec<usize>,
}

impl std::fmt::Display for Focus {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let &Self {
            from_node,
            from_position,
            to_node,
            to_position,
            dist,
            log_likelihood_ratio,
            ref counts,
            ..
        } = self;
        let (fu, fc) = from_node;
        let (tu, tc) = to_node;
        write!(
            f,
            "{fu},{fc},{from_position},{tu},{tc},{to_position},{dist},{:.2},{:?}",
            log_likelihood_ratio, counts
        )
    }
}

impl Focus {
    /// LogLilihood ratio
    fn llr(&self) -> f64 {
        self.log_likelihood_ratio
    }
    fn new(
        (from, from_node, from_position): (NodeIndex, Node, Position),
        (to, to_node, to_position): (NodeIndex, Node, Position),
        dist: usize,
        lk: f64,
        counts: Vec<usize>,
    ) -> Self {
        Self {
            from,
            from_node,
            from_position,
            to,
            to_node,
            to_position,
            dist,
            log_likelihood_ratio: lk,
            counts,
        }
    }
}

use std::collections::HashSet;
impl<'b, 'a: 'b> DitchGraph<'a> {
    /// Check the dynamic of foci.
    fn survey_foci(&'b mut self, foci: &[Focus]) -> usize {
        let mut solved = 0;
        let mut affected_nodes: HashSet<NodeIndex> = HashSet::new();
        for focus in foci {
            if self.is_deleted(focus.from) || self.is_deleted(focus.to) {
                continue;
            }
            if affected_nodes.contains(&focus.from) || affected_nodes.contains(&focus.to) {
                continue;
            }
            // Check if the targets are present in the graph.
            let from_copy_num = self.node(focus.from).and_then(|n| n.copy_number);
            let to_copy_num = self.node(focus.to).and_then(|node| node.copy_number);
            if !matches!(from_copy_num, Some(1)) || !matches!(to_copy_num, Some(1)) {
                continue;
            }
            // Check if this is branching.
            let (key, pos) = (focus.from, focus.from_position);
            if self.count_edges(key, pos) != 1 {
                continue;
            }
            // Check this branch has two-in-degree.
            let mut edges = self.edges_from(key, pos);
            let edge = edges.pop().unwrap();
            let sibs = self.count_edges(edge.to, edge.to_position);
            if sibs <= 1 {
                continue;
            }
            debug!("FOCUS\tTRY\t{}", focus);
            solved += self.survey_focus(focus, &mut affected_nodes).is_some() as usize;
        }
        solved
    }
    // If resolveing successed, return Some(())
    // othewise, it encountered a stronger focus, and return None.
    fn survey_focus(&'b mut self, focus: &Focus, affected: &mut HashSet<NodeIndex>) -> Option<()> {
        let path_to_focus = self.dfs_to_the_target(focus)?;
        if path_to_focus.iter().any(|(n, _)| affected.contains(n)) {
            return None;
        }
        // If this path is no-branching, it is no use.
        if !self.is_path_branching(focus, &path_to_focus) {
            return None;
        }
        for (i, &(index, pos)) in path_to_focus.iter().enumerate() {
            let node = self.node(index).unwrap().node;
            debug!("FOCUS\tTrace\t{}\t{:?}\t{}", i, node, pos);
        }
        let new_nodes = self.duplicate_along(focus, &path_to_focus);
        affected.extend(new_nodes);
        affected.extend(path_to_focus.iter().map(|x| x.0));
        self.remove_along(focus, &path_to_focus);
        Some(())
    }
    fn is_path_branching(&self, focus: &Focus, path: &[(NodeIndex, Position)]) -> bool {
        let len = path.len();
        let mut is_branching = false;
        is_branching |= 1 < self.count_edges(focus.from, focus.from_position);
        let mut check_both = path.iter().take(len.saturating_sub(1));
        is_branching |= check_both.any(|&(node, _)| {
            let edges = self.node(node).unwrap().edges.iter();
            let (head, tail) = edges.fold((0, 0), |(h, t), e| match e.from_position {
                Position::Head => (h + 1, t),
                Position::Tail => (h, t + 1),
            });
            1 < head || 1 < tail
        });
        if let Some(&(node, position)) = path.last() {
            is_branching |= 1 < self.count_edges(node, position);
        }
        is_branching
    }
    // Return nodes newly allocated.
    fn duplicate_along(
        &'b mut self,
        focus: &Focus,
        path_spaning: &[(NodeIndex, Position)],
    ) -> Vec<NodeIndex> {
        let (mut c_index, mut c_pos) = (focus.from, focus.from_position);
        let mut prev_id = focus.from;
        let mut newly_allocated = vec![];
        for (i, &(to_id, to_pos)) in path_spaning.iter().enumerate() {
            // Duplicate if nessesary.
            let new_idx = match i + 1 == path_spaning.len() {
                true => to_id,
                false => {
                    let occ = self
                        .node_mut(to_id)
                        .map(|n| n.decrement_copy_number())
                        .unwrap_or(0);
                    let new_idx = self.duplicate(to_id);
                    if let Some(n) = self.node_mut(new_idx) {
                        n.occ = occ;
                        n.copy_number = Some(1);
                    }
                    newly_allocated.push(new_idx);
                    new_idx
                }
            };
            let edge_occ = self.decrement_edge_copy_number((c_index, c_pos), (to_id, to_pos));
            let mut edge = self
                .node(c_index)
                .unwrap()
                .edges
                .iter()
                .find(|e| (e.from_position, e.to, e.to_position) == (c_pos, to_id, to_pos))
                .unwrap()
                .clone();
            edge.occ = edge_occ;
            edge.copy_number = Some(1);
            edge.from = prev_id;
            edge.to = new_idx;
            self.add_edge(edge);
            c_index = to_id;
            c_pos = !to_pos;
            prev_id = new_idx;
        }
        newly_allocated
    }
    // fn duplicate_along(&'b mut self, focus: &Focus, path_spanning: &[(NodeIndex, Position)]) {
    //     let mut new_id_runs = vec![focus.from];
    //     for &(to_id, _) in path_spanning.iter() {
    //         if to_id != focus.to {
    //             let occ = self
    //                 .node_mut(to_id)
    //                 .map(|n| n.decrement_copy_number())
    //                 .unwrap_or(0);
    //             let new_idx = self.duplicate(to_id);
    //             let new_node = self.node_mut(new_idx).unwrap();
    //             new_node.occ = occ;
    //             new_node.copy_number = Some(1);
    //             new_id_runs.push(new_idx);
    //         } else {
    //             new_id_runs.push(to_id);
    //         }
    //     }
    //     assert_eq!(new_id_runs.len(), path_spanning.len() + 1);
    //     let (mut c_index, mut c_pos) = (focus.from, focus.from_position);
    //     for (window, &(old_idx, old_pos)) in new_id_runs.windows(2).zip(path_spanning.iter()) {
    //         let occ = self.decrement_edge_copy_number((c_index, c_pos), (old_idx, old_pos));
    //         let mut edge = self
    //             .node(c_index)
    //             .unwrap()
    //             .edges
    //             .iter()
    //             .find(|e| (e.from_position, e.to, e.to_position) == (c_pos, old_idx, old_pos))
    //             .unwrap()
    //             .clone();
    //         edge.occ = occ;
    //         edge.copy_number = Some(1);
    //         edge.from = window[0];
    //         edge.to = window[1];
    //         self.add_edge(edge);
    //         c_index = old_idx;
    //         c_pos = !old_pos;
    //     }
    // }
    fn remove_along(&'b mut self, focus: &Focus, path_spaning: &[(NodeIndex, Position)]) {
        // Removing zero-copy edges along the path.
        let (mut prev, mut prev_pos) = (focus.from, focus.from_position);
        for &(node, pos) in path_spaning {
            let prev_node = self.node(prev).unwrap();
            let copy_num = prev_node
                .edges
                .iter()
                .find(|e| e.from_position == prev_pos && e.to == node && e.to_position == pos)
                .and_then(|e| e.copy_number);
            if copy_num == Some(0) {
                self.remove_edge_between((prev, prev_pos), (node, pos));
            }
            prev = node;
            prev_pos = !pos;
        }
        // Removing other branching edges along the path.
        let remove: Vec<_> = {
            let mut remove_edges = self.edges_from(focus.from, focus.from_position);
            remove_edges.retain(|edge| edge.copy_number == Some(0));
            remove_edges.iter().map(|e| e.key()).collect()
        };
        for (f, t) in remove {
            self.remove_edge_between(f, t);
        }
        for &(node, _) in path_spaning {
            let remove: Vec<_> = self
                .node(node)
                .unwrap()
                .edges
                .iter()
                .filter(|edge| edge.copy_number == Some(0))
                .map(|e| ((e.from, e.from_position), (e.to, e.to_position)))
                .collect();
            for (f, t) in remove {
                self.remove_edge_between(f, t);
            }
        }
        for &(node, _) in path_spaning {
            // node might be absent, because it is included in a loop!
            if !self.is_deleted(node) {
                let edge_num = self.node(node).unwrap().edges.len();
                if edge_num == 0 {
                    self.delete(node);
                }
            }
        }
    }
    // Span the region by `path_spaning`.
    // fn span_region(&'b mut self, focus: &Focus, path_spaning: &[(NodeIndex, Position)]) {
    //     assert_eq!(*path_spaning.last().unwrap(), (focus.to, focus.to_position));
    //     // 1. Duplicate nodes along with the path, except the last node.
    //     let mut duplicated_nodes = vec![];
    //     for &(node, _) in path_spaning.iter().take(len - 1) {
    //         let new_idx = self.duplicate(node);
    //         duplicated_nodes.push(new_idx);
    //     }
    //     // 2. Remove all the edges in the newly allocated nodes.
    //     for &index in duplicated_nodes.iter() {
    //         self.node_mut(index).unwrap().edges.clear()
    //     }
    //     // 3. Add edges appropriately.
    //     let (mut from, mut from_pos) = (focus.from, focus.from_position);
    //     for (&(old, old_pos), &new) in path_spaning.iter().zip(duplicated_nodes.iter()) {
    //         let from_node = self.node(from).unwrap();
    //         let mut new_edge = from_node
    //             .edges
    //             .iter()
    //             .find(|e| e.from_position == from_pos && e.to == old && e.to_position == old_pos)
    //             .unwrap()
    //             .clone();
    //         new_edge.to = new;
    //         self.add_edge(new_edge);
    //         from = old;
    //         from_pos = !old_pos;
    //     }
    //     // 4. removing zero-copy edges along the path.
    //     let (mut prev, mut prev_pos) = (focus.from, focus.from_position);
    //     for &(node, pos) in path_spaning {
    //         let prev_node = self.node(prev).unwrap();
    //         let edge = prev_node
    //             .edges
    //             .iter()
    //             .find(|e| e.from_position == prev_pos && e.to == node && e.to_position == pos);
    //         let copy_num = edge.and_then(|e| e.copy_number);
    //         if copy_num == Some(0) {
    //             self.remove_edge_between((prev, prev_pos), (node, pos));
    //         }
    //         prev = node;
    //         prev_pos = !pos;
    //     }
    //     // Removing other branching edges along the path.
    //     let remove: Vec<_> = {
    //         let mut remove_edges = self.edges_from(from, from_pos);
    //         remove_edges.retain(|edge| edge.copy_number == Some(0));
    //         remove_edges.iter().map(|e| e.key()).collect()
    //     };
    //     for (f, t) in remove {
    //         self.remove_edge_between(f, t);
    //     }
    //     for &(node, _) in path_spaning {
    //         let remove: Vec<_> = self
    //             .node(node)
    //             .unwrap()
    //             .edges
    //             .iter()
    //             .filter(|edge| edge.copy_number == Some(0))
    //             .map(|e| ((e.from, e.from_position), (e.to, e.to_position)))
    //             .collect();
    //         for (f, t) in remove {
    //             self.remove_edge_between(f, t);
    //         }
    //     }
    //     for &(node, _) in path_spaning {
    //         // node might be absent, because it is included in a loop!
    //         if !self.is_deleted(node) {
    //             let edge_num = self.node(node).unwrap().edges.len();
    //             if edge_num == 0 {
    //                 self.delete(node);
    //             }
    //         }
    //     }
    // }
    // From (start,pos) position to the focus target(start=focus.node, position=focus.position).
    // Note that, the path would contain the (target,target pos) at the end of the path and
    // not contain the (start,pos) position itself.
    // TODO: Maybe this function failed to find the target focus,
    // because of the removed nodes.
    // If it happened, eigther this focus or previous
    // branching is false. Maybe just dropping this focus would be better?
    // Note that this function does not include any proxying nodes.
    fn dfs_to_the_target(&self, focus: &Focus) -> Option<Vec<(NodeIndex, Position)>> {
        // Breadth first search until get to the target.
        // Just !pos to make code tidy.
        let mut nodes_at = vec![vec![]; focus.dist + 1];
        let mut parents = vec![vec![]; focus.dist + 1];
        nodes_at[0].push((focus.from, !focus.from_position));
        // (distance, index, copy number of edge)
        parents[0].push((0, 0, 0));
        for dist in 0..focus.dist + 1 {
            // To avoid the borrow checker.
            // TODO: Do not need filter...
            let mut nodes_at_temp = vec![];
            for (idx, &(index, pos)) in nodes_at[dist].iter().enumerate() {
                let edges = self
                    .edges_from(index, !pos)
                    .into_iter()
                    .map(|e| (e, dist + 1));
                for (edge, d) in edges.filter(|&(_, d)| d <= focus.dist) {
                    let cp = edge.copy_number.unwrap_or(0);
                    nodes_at_temp.push((d, (edge.to, edge.to_position)));
                    parents[d].push((dist, idx, cp));
                }
            }
            for (d, elm) in nodes_at_temp {
                nodes_at[d].push(elm);
            }
        }
        assert_eq!(nodes_at.len(), parents.len());
        nodes_at
            .iter()
            .zip(parents.iter())
            .for_each(|(n, p)| assert_eq!(n.len(), p.len()));
        // Back-track.
        let mut target = (focus.to, focus.to_position);
        let mut dist = focus.dist;
        let mut back_track = vec![];
        while 0 < dist {
            back_track.push(target);
            let hit_at = nodes_at[dist]
                .iter()
                .zip(parents[dist].iter())
                .filter(|(&x, _)| x == target)
                .max_by_key(|(_, (_, _, cp))| cp);
            (target, dist) = match hit_at {
                Some((_, &(pdist, pidx, _))) => (nodes_at[pdist][pidx], pdist),
                None => {
                    debug!("FOCUS\tFailedToFindTarget\t{}", focus);
                    return None;
                }
            };
        }
        back_track.reverse();
        Some(back_track)
    }
    // Spell along the given path,
    // reducing the copy number of the nodes, the occs of the nodes, and the
    // occs of the edges.
    // fn spell_along_path(&'b mut self, focus: &Focus, path: &[(NodeIndex, Position)]) {
    //     let (mut c_index, mut c_pos) = (focus.from, focus.from_position);
    //     for &(to_id, to_pos) in path.iter() {
    //         if c_index != focus.from {
    //             self.node_mut(c_index)
    //                 .map(|n| n.decrement_copy_number())
    //                 .unwrap_or(0);
    //         }
    //         self.decrement_edge_copy_number((c_index, c_pos), (to_id, to_pos));
    //         // Reduce the occ and copy number of the node, and occ of the edge.
    //         c_index = to_id;
    //         c_pos = !to_pos;
    //     }
    // }

    /// Resolve repeats by "foci" algorithm.
    /// The algorithm consists three step.
    /// 1. For each branching node u, find another node v where reads outgoing from u are
    /// "gathering". In other words, we check the likelihood ratio between "gathering"-model
    /// and null-model and check whether or not the likelihood ration is larger than `thr`.
    /// 2. Select an arbitrary path u -> v, make an new edge connecting u and v,
    /// and remove all other edges from u and to v. By doing this, u and v would be in the same simple path(repeat resolved).
    /// 3. Recursively remove nodes and their edges, whose in-degree are zero.
    /// TODO: This code has bug when the graph contains loops?
    /// TODO: This code has bug, which breaks up a continuous graph ...?
    pub fn resolve_repeats(
        &'b mut self,
        reads: &[&EncodedRead],
        config: &AssembleConfig,
        thr: f64,
    ) {
        debug!("FOCI\tRESOLVE\t{:.3}\t{}", thr, config.min_span_reads);
        warn!("Please fix lk computation!");
        let mut count = 1;
        while count != 0 {
            let mut foci = self.get_foci(reads, config);
            let prev = foci.len();
            foci.retain(|val| thr < val.llr());
            foci.retain(|val| val.llr() != std::f64::INFINITY);
            if foci.is_empty() {
                break;
            }
            debug!("FOCI\tNUM\t{}\t{}", prev, foci.len());
            foci.sort_by(|x, y| match y.llr().partial_cmp(&x.llr()).unwrap() {
                std::cmp::Ordering::Equal => (y.dist, y.from).cmp(&(x.dist, x.from)),
                x => x,
            });
            count = self.survey_foci(&foci);
            debug!("FOCI\tTryAndSuccess\t{}\t{}", foci.len(), count);
        }
    }

    /// Return a hash map containing all foci, thresholded by `config` parameter.
    /// For each (node, position), keep the strongest focus.
    fn get_foci(&self, reads: &[&EncodedRead], config: &AssembleConfig) -> Vec<Focus> {
        let mut foci = vec![];
        for (index, node) in self
            .nodes()
            .filter(|(_, n)| matches!(n.copy_number, Some(1)))
        {
            // Find
            // --Node(Copy=1)--|
            //                 |--Node(Copy>1)--
            // --Node----------|
            //  (We are here now)
            for pos in [Position::Head, Position::Tail] {
                if self.count_edges(index, pos) != 1 {
                    continue;
                }
                let edge = node.edges.iter().find(|e| e.from_position == pos).unwrap();
                let siblings = self.count_edges(edge.to, edge.to_position);
                assert!(1 <= siblings);
                if siblings == 1 {
                    continue;
                }
                // If this edge does not flow into multi-copy contig, continue.
                match self.node(edge.to).unwrap().copy_number {
                    None | Some(1) | Some(0) => continue,
                    _ => {}
                }
                if let Some(focus) = self.examine_focus(index, node, pos, reads, config) {
                    foci.push(focus);
                }
            }
        }
        foci
    }
    // Return the focus if given node has a focus. Otherwise, return None.
    fn examine_focus(
        &self,
        node_index: NodeIndex,
        node: &DitchNode,
        pos: Position,
        reads: &[&EncodedRead],
        config: &AssembleConfig,
    ) -> Option<Focus> {
        // TODO:Maybe we can faster this process by using reverse indexing.
        let reads: Vec<_> = reads
            .iter()
            .filter(|r| r.contains(node.node))
            .copied()
            .collect();
        let dist_nodes =
            self.enumerate_candidate_nodes(&reads, config.min_span_reads, node_index, pos);
        let mut focus: Option<Focus> = None;
        for (dist, node_indices) in dist_nodes.iter().enumerate().filter(|(_, ns)| 1 < ns.len()) {
            let nodes: Vec<_> = node_indices
                .iter()
                .map(|(n, _)| self.node(*n).unwrap())
                .collect();
            if nodes.iter().map(|n| n.occ).sum::<usize>() <= config.min_span_reads {
                continue;
            }
            let dist = dist + 1;
            let occs = Self::count_occurences((node, pos), dist, &nodes, &reads);
            let nodes_with_lk = Self::max_lk_node(&occs, &nodes, node_indices);
            if let Some((lk_ratio, to, (to_idx, to_pos))) = nodes_with_lk {
                if focus.as_ref().map(|f| f.llr()).unwrap_or(0f64) < lk_ratio {
                    let to = (to_idx, to.node, to_pos);
                    let from = (node_index, node.node, pos);
                    let f = Focus::new(from, to, dist, lk_ratio, occs);
                    focus.replace(f);
                }
            }
        }
        focus
    }
    fn count_occurences(
        (node, pos): (&DitchNode, Position),
        dist: usize,
        nodes: &[&DitchNode],
        reads: &[&EncodedRead],
    ) -> Vec<usize> {
        let mut occs: Vec<_> = vec![0; nodes.len()];
        for read in reads.iter() {
            let start = read
                .nodes
                .iter()
                .position(|n| (n.unit, n.cluster) == node.node)
                .unwrap();
            let check_node = match (read.nodes[start].is_forward, pos == Position::Tail) {
                (true, true) | (false, false) if start + dist < read.nodes.len() => {
                    &read.nodes[start + dist]
                }
                (false, true) | (true, false) if dist <= start => &read.nodes[start - dist],
                _ => continue,
            };
            let check_node = (check_node.unit, check_node.cluster);
            if let Some(hit_node) = nodes.iter().position(|n| n.node == check_node) {
                occs[hit_node] += 1;
            }
        }
        occs
    }
    #[allow(dead_code)]
    fn max_lk_node_dev<'c>(
        occs: &[usize],
        nodes: &[&'c DitchNode<'a>],
        node_indices: &[(NodeIndex, Position)],
    ) -> Option<(f64, &'c DitchNode<'a>, (NodeIndex, Position))> {
        let cov_total: f64 = nodes.iter().map(|n| n.occ as f64).sum();
        let occ_total: f64 = occs.iter().map(|&occ| occ as f64).sum();
        nodes
            .iter()
            .zip(node_indices.iter())
            .zip(occs.iter())
            .filter(|((n, _), &occ)| 0 < occ && matches!(n.copy_number, Some(1)))
            .map(|((n, idx), &occ)| {
                let occ = occ as f64;
                let null_lk = Self::likelihood_with(occs, occ, n.occ as f64 / cov_total);
                let alt_lk = Self::likelihood_with(occs, occ, occ / occ_total);
                (alt_lk - null_lk, *n, *idx)
            })
            .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
    }
    fn likelihood_with(occs: &[usize], occ: f64, prob: f64) -> f64 {
        assert!(prob <= 1f64);
        assert!(0f64 <= prob);
        // modify so that .ln() would not fail
        let prob = prob.max(0.0001).min(0.999);
        let lk: f64 = occs
            .iter()
            .map(|&occ| occ as f64 * (1f64 - prob).ln())
            .sum();
        lk + occ * prob.ln() - occ * (1f64 - prob).ln()
    }
    #[allow(dead_code)]
    fn max_lk_node<'c>(
        occs: &[usize],
        nodes: &[&'c DitchNode<'a>],
        node_indices: &[(NodeIndex, Position)],
    ) -> Option<(f64, &'c DitchNode<'a>, (NodeIndex, Position))> {
        let ith_ln = {
            let mut probs: Vec<_> = nodes.iter().map(|n| n.occ as f64).collect();
            let sum: f64 = probs.iter().sum();
            probs.iter_mut().for_each(|x| *x = (*x / sum).ln());
            assert!(probs.iter().all(|x| !x.is_nan()), "{:?}", probs);
            probs
        };
        let null_prob: f64 = occs
            .iter()
            .zip(ith_ln.iter())
            .map(|(&occ, ln)| match occ {
                0 => 0f64,
                _ => occ as f64 * ln,
            })
            .sum();
        assert!(!null_prob.is_nan(), "{:?}\t{:?}", ith_ln, occs);
        assert!(null_prob < std::f64::INFINITY, "{:?}\t{:?}", ith_ln, occs);
        let (correct_lk, error_lk) = lk_pairs(nodes.len());
        assert_eq!(nodes.len(), node_indices.len());
        nodes
            .iter()
            .zip(node_indices.iter())
            .enumerate()
            .filter(|(_, (n, _))| n.copy_number == Some(1))
            .map(|(k, (n, idx))| {
                let occs_with_idx = occs.iter().map(|&x| x as f64).enumerate();
                let lk: f64 = occs_with_idx
                    .map(|(i, occ)| occ * (if i == k { correct_lk } else { error_lk }))
                    .sum();
                assert!(!lk.is_nan());
                (lk - null_prob, *n, *idx)
            })
            .max_by(|x, y| (x.0).partial_cmp(&(y.0)).unwrap())
    }
    fn enumerate_candidate_nodes(
        &self,
        reads: &[&EncodedRead],
        min_span: usize,
        node_index: NodeIndex,
        pos: Position,
    ) -> Vec<Vec<(NodeIndex, Position)>> {
        let node = self.node(node_index).unwrap().node;
        let max_reach = Self::max_reach(reads, min_span, node, pos);
        let mut nodes_at: Vec<Vec<(NodeIndex, Position)>> = vec![];
        let mut active_nodes = vec![(node_index, pos)];
        for _ in 0..max_reach {
            let mut found_nodes: Vec<(NodeIndex, Position)> = vec![];
            for &(index, pos) in active_nodes.iter() {
                for edge in self.edges_from(index, pos) {
                    found_nodes.push((edge.to, edge.to_position));
                }
            }
            found_nodes.sort_unstable();
            found_nodes.dedup();
            active_nodes.clear();
            active_nodes.extend(found_nodes.iter().map(|&(i, p)| (i, !p)));
            nodes_at.push(found_nodes);
        }
        // TODO: I remove duplicated nodes here, but it may cause some problems I do not know currently.
        for nodes in nodes_at.iter_mut() {
            let duplicated: Vec<_> = nodes
                .iter()
                .enumerate()
                .map(|(i, &(n, _))| {
                    nodes[i + 1..]
                        .iter()
                        .any(|&(m, _)| self.node(n).unwrap().node == self.node(m).unwrap().node)
                })
                .collect();
            let mut idx = 0;
            nodes.retain(|_| {
                idx += 1;
                !duplicated[idx - 1]
            });
        }
        nodes_at
    }
    fn max_reach(reads: &[&EncodedRead], min_span: usize, node: Node, pos: Position) -> usize {
        let mut dist_node: HashMap<_, _> = HashMap::new();
        for read in reads.iter() {
            let start = read
                .nodes
                .iter()
                .position(|n| (n.unit, n.cluster) == node)
                .unwrap();
            match (read.nodes[start].is_forward, pos) {
                (true, Position::Tail) | (false, Position::Head) => {
                    // ---> (Here to start) -- <---> -- <---> .... Or
                    // <--- (Here to start) -- <---> -- <---> ....
                    for (dist, node) in read.nodes.iter().skip(start).enumerate() {
                        let node = (node.unit, node.cluster);
                        *dist_node.entry((node, dist)).or_default() += 1;
                    }
                }
                (true, Position::Head) | (false, Position::Tail) => {
                    let read = read.nodes.iter().take(start + 1).rev();
                    for (dist, node) in read.enumerate() {
                        let node = (node.unit, node.cluster);
                        *dist_node.entry((node, dist)).or_default() += 1;
                    }
                }
            }
        }
        dist_node
            .iter()
            .filter_map(|((_, d), &count)| (min_span <= count).then(|| *d))
            .max()
            .unwrap_or(0)
    }
}

fn lk_pairs(len: usize) -> (f64, f64) {
    let choice_num = len as f64;
    let correct_lk = ((1f64 - ERROR_PROB).powi(2) + ERROR_PROB / choice_num).ln();
    let error_lk = {
        let correct_to_error = (1.0 - ERROR_PROB) * ERROR_PROB * (choice_num - 1.0) / choice_num;
        let error_to_error =
            ERROR_PROB / choice_num * ERROR_PROB * (choice_num - 1f64) / choice_num;
        (correct_to_error + error_to_error).ln()
    };
    assert!(!correct_lk.is_nan());
    assert!(!error_lk.is_nan());
    (correct_lk, error_lk)
}
