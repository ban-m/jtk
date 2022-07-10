// selecting a branche from branches by using phaing sets.
// In other words, whenever faced with branches,
// the algorithm searches further node in the same phaseset.
// If the search successed within `len` time BFS,
// it select an arbitrary path from the current branching point to the found node,
// removing other unnessesary edges.
// pub fn resolve_by_phaseset(&mut self, phaseset: &[HashSet<Node>], len: usize) -> Option<()> {
//     let mut cand_nodes: Vec<_> = vec![];
//     for node in self.nodes().filter(|n| n.copy_number == Some(1)) {
//         cand_nodes.push((node.node, Position::Tail));
//         cand_nodes.push((node.node, Position::Head));
//     }
//     // Make phase-set into node->phaseblock.
//     let phaseset: HashMap<(u64, u64), usize> = phaseset
//         .iter()
//         .enumerate()
//         .flat_map(|(block, nodes)| nodes.iter().map(|&node| (node, block)).collect::<Vec<_>>())
//         .collect();
//     for (from, from_pos) in cand_nodes {
//         // Check if this is branching.
//         let have_sibs = self
//             .get_edges(from, from_pos)
//             .all(|edge| 1 < self.get_edges(edge.to, edge.to_position).count());
//         let num_edges = self.get_edges(from, from_pos).count();
//         if num_edges < 2 || !have_sibs {
//             continue;
//         }
//         // resolving branches.
//         let from = (from, from_pos);
//         let path_phasing = self.bfs_to_the_same_phase(from, len, &phaseset)?;
//         let to = *path_phasing.last().unwrap();
//         let edge = self.spell_along_path(&path_phasing, from);
//         debug!(
//             "PHASEPATH\tSpan\t{:?}\t{:?}\t{}",
//             from.0,
//             to.0,
//             edge.seq.len()
//         );
//         self.span_region(from, to, &path_phasing, edge);
//     }
//     Some(())
// }

// Fron (start,pos) position to the node with the same phase block.
// Note that, the path would contain the (target,target pos) at the end of the path and
// not contain the (start,pos) position itself.
// fn bfs_to_the_same_phase(
//     &self,
//     (start, pos): (Node, Position),
//     len: usize,
//     phaseset: &HashMap<Node, usize>,
// ) -> Option<Vec<(Node, Position)>> {
//     // Breadth first search until get to the target.
//     // Just !pos to make code tidy.
//     let target_phase = Some(phaseset.get(&start)?);
//     let mut nodes_at = vec![vec![(start, !pos)]];
//     let mut parents = vec![vec![]];
//     'outer: for dist in 0..len {
//         let mut next_nodes = vec![];
//         let mut parent = vec![];
//         assert_eq!(dist + 1, nodes_at.len());
//         for (idx, &(node, pos)) in nodes_at[dist].iter().enumerate() {
//             if phaseset.get(&node) == target_phase && node != start {
//                 break 'outer;
//             }
//             // Move to the end of the node.
//             let pos = !pos;
//             for edge in self.get_edges(node, pos) {
//                 next_nodes.push((edge.to, edge.to_position));
//                 parent.push(idx);
//             }
//         }
//         assert_eq!(next_nodes.len(), parent.len());
//         nodes_at.push(next_nodes);
//         parents.push(parent);
//     }
//     // Back-track.
//     assert_eq!(nodes_at.len(), parents.len());
//     let mut dist = nodes_at.len() - 1;
//     let mut idx = match nodes_at[dist]
//         .iter()
//         .position(|(n, _)| phaseset.get(n) == target_phase)
//     {
//         Some(idx) => idx,
//         None => {
//             debug!("PHASEPATH\tFailedToFindTarget\t{:?}", start);
//             return None;
//         }
//     };
//     let mut back_track = vec![];
//     while 0 < dist {
//         back_track.push(nodes_at[dist][idx]);
//         idx = parents[dist][idx];
//         dist -= 1;
//     }
//     back_track.reverse();
//     for (i, (node, pos)) in back_track.iter().enumerate() {
//         debug!("PHASEPATH\tTrace\t{}\t{:?}\t{}", i, node, pos);
//     }
//     Some(back_track)
// }
