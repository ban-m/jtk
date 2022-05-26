//!
//! `Indexable` trait
//! Abstraction of types that can be used as an index of vector
//!
pub use petgraph::graph::{EdgeIndex, NodeIndex};

/// Generalized indexable types
/// that have one-to-one correspondence with `usize`
pub trait Indexable: Copy + Eq + std::hash::Hash + std::fmt::Debug {
    /// convert `usize` to the `Indexable` type
    fn new(x: usize) -> Self;
    /// convert-back the `Indexable` type to `usize`
    fn index(&self) -> usize;
}

impl Indexable for usize {
    #[inline]
    fn new(x: usize) -> Self {
        x
    }
    #[inline]
    fn index(&self) -> usize {
        *self
    }
}

impl Indexable for NodeIndex {
    #[inline]
    fn new(x: usize) -> Self {
        NodeIndex::new(x)
    }
    #[inline]
    fn index(&self) -> usize {
        NodeIndex::index(*self)
    }
}

impl Indexable for EdgeIndex {
    #[inline]
    fn new(x: usize) -> Self {
        EdgeIndex::new(x)
    }
    #[inline]
    fn index(&self) -> usize {
        EdgeIndex::index(*self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexable() {
        let x: usize = Indexable::new(11);
        assert_eq!(x.index(), 11);
        let x: NodeIndex = Indexable::new(9);
        assert_eq!(x.index(), 9);
        let x: EdgeIndex = Indexable::new(9);
        assert_eq!(x.index(), 9);
    }
}
