//! Iterators on a ditch graph.
//! The access is currently serial. In other words, `nth(n)` has O(n) complexity.
//! Do not use usual iterator directory, i.e., `self.nodes.iter()`, as it would include nodes already deleted.
use super::DitchEdge;
use super::DitchGraph;
use super::DitchNode;
use super::NodeIndex;
impl<'b> DitchGraph<'b> {
    pub fn nodes(&self) -> NodeIter {
        NodeIter::new(self)
    }
    /// Modify each node by `f`.
    /// Because implememting `iter_mut` requires unsafe codes,
    /// we restrict users to only modify the elements by this `black-box` modifier.
    /// This function is not parallelized.
    pub fn modify_by<F>(&mut self, f: F)
    where
        F: FnMut(&mut DitchNode<'b>),
    {
        self.nodes.iter_mut().filter(|n| n.is_deleted).for_each(f);
    }
    pub fn edges(&self) -> impl std::iter::Iterator<Item = &DitchEdge> {
        self.nodes().flat_map(|(_, n)| n.edges.iter())
    }
    pub fn modify_by_with_index<F>(&mut self, f: F)
    where
        F: FnMut((NodeIndex, &mut DitchNode<'b>)),
    {
        self.nodes
            .iter_mut()
            .enumerate()
            .map(|(i, n)| (NodeIndex(i), n))
            .filter(|n| !n.1.is_deleted)
            .for_each(f)
    }
}

/// Iterating over nodes and its index (unique identifier)
#[derive(Debug, Clone)]
pub struct NodeIter<'a, 'b> {
    graph: &'a DitchGraph<'b>,
    position: usize,
}

impl<'a, 'b> NodeIter<'a, 'b> {
    fn new(graph: &'a DitchGraph<'b>) -> Self {
        Self { graph, position: 0 }
    }
}

impl<'a, 'b> std::iter::Iterator for NodeIter<'a, 'b> {
    type Item = (NodeIndex, &'a DitchNode<'b>);
    fn next(&mut self) -> Option<Self::Item> {
        while let Some(node) = self.graph.nodes.get(self.position) {
            self.position += 1;
            if !node.is_deleted {
                return Some((NodeIndex(self.position - 1), node));
            }
        }
        None
    }
}
