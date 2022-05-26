//!
//! measureing difference between two vectors
//!
use super::{Indexable, Storage, Vector};
use crate::prob::Prob;
use itertools::{chain, izip, Itertools};

impl<S, Ix> Vector<S, Ix>
where
    S: Storage<Item = f64>,
    Ix: Indexable,
{
    ///
    /// diff of two vectors whose element is f64.
    ///
    pub fn diff_f64<T>(&self, other: &Vector<T, Ix>) -> f64
    where
        T: Storage<Item = f64>,
    {
        assert_eq!(self.len(), other.len());
        let mut diff = 0f64;
        for i in 0..self.len() {
            diff += (self[Ix::new(i)] - other[Ix::new(i)]).abs();
        }
        diff
    }
}

impl<S, Ix> Vector<S, Ix>
where
    S: Storage<Item = Prob>,
    Ix: Indexable,
{
    ///
    /// return the maximum probability and its index of the vector.
    ///
    pub fn max(&self) -> Option<(Ix, Prob)> {
        self.iter_all().max_by_key(|(_, p)| *p)
    }
    ///
    /// return the minimum probability and its index of the vector.
    ///
    pub fn min(&self) -> Option<(Ix, Prob)> {
        self.iter_all().min_by_key(|(_, p)| *p)
    }
    ///
    /// get the sorted iterator on element (in a descending order)
    ///
    pub fn iter_sorted_desc<'a>(&'a self) -> impl Iterator<Item = (Ix, Prob)> + 'a {
        self.iter_all()
            .sorted_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap())
    }
    ///
    /// check if elements in vector whose probability satisfies `p > min_prob`
    /// has same index?
    ///
    pub fn is_same_ranking<T>(&self, other: &Vector<T, Ix>, min_prob: Prob) -> bool
    where
        T: Storage<Item = Prob>,
    {
        for ((ia, _pa), (ib, _pb)) in izip!(
            self.iter_sorted_desc().filter(|(_, p)| *p >= min_prob),
            other.iter_sorted_desc().filter(|(_, p)| *p >= min_prob)
        ) {
            // println!("{:?} {} {:?} {}", ia, pa, ib, pb);
            if ia != ib {
                return false;
            }
        }
        true
    }
    ///
    /// Log difference score of two vectors whose item is prob.
    /// it checks only on active values
    ///
    pub fn diff_on_active<T>(&self, other: &Vector<T, Ix>) -> f64
    where
        T: Storage<Item = Prob>,
    {
        let mut diff = 0f64;
        for (_, d) in self.iter_diff_on_active(other) {
            diff += d
        }
        diff
    }
    ///
    /// Log difference score of two vectors whose item is prob.
    ///
    pub fn diff<T>(&self, other: &Vector<T, Ix>) -> f64
    where
        T: Storage<Item = Prob>,
    {
        let mut diff = 0f64;
        for (_, d) in self.iter_diff(other) {
            diff += d
        }
        diff
    }
    ///
    /// Calculate `\sum_i |p[i]-p[i]|` for all active i
    ///
    pub fn iter_diff_on_active<'a, T>(
        &'a self,
        other: &'a Vector<T, Ix>,
    ) -> impl Iterator<Item = (Ix, f64)> + 'a
    where
        T: Storage<Item = Prob>,
    {
        assert_eq!(self.len(), other.len());
        chain!(self.iter().map(|(i, _)| i), other.iter().map(|(i, _)| i))
            .unique()
            .map(move |i| {
                let diff = (self[i].to_value() - other[i].to_value()).abs();
                (i, diff)
            })
    }
    ///
    /// Calculate `\sum_i |p[i]-p[i]|` for all i
    ///
    pub fn iter_diff<'a, T>(
        &'a self,
        other: &'a Vector<T, Ix>,
    ) -> impl Iterator<Item = (Ix, f64)> + 'a
    where
        T: Storage<Item = Prob>,
    {
        assert_eq!(self.len(), other.len());
        (0..self.len()).map(move |i| {
            let index = Ix::new(i);
            let diff = (self[index].to_value() - other[index].to_value()).abs();
            (index, diff)
        })
    }
    ///
    /// visualize difference of two vectors
    ///
    pub fn show_diff<T>(&self, other: &Vector<T, Ix>)
    where
        T: Storage<Item = Prob>,
    {
        for (i, diff) in self
            .iter_diff(other)
            .sorted_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap())
            .take(100)
        {
            println!("{}: {}", i.index(), diff);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prob::p;
    use crate::vector::{DenseStorage, SparseStorage};

    #[test]
    fn vector_diff() {
        let mut v1: Vector<DenseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v1[0] = p(0.3);
        let mut v2: Vector<DenseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.5);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);

        let mut v1: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v1[0] = p(0.3);
        let mut v2: Vector<SparseStorage<Prob>, usize> = Vector::new(5, p(0.0));
        v2[1] = p(0.5);
        let d = v1.diff(&v2);
        assert_abs_diff_eq!(d, 0.8);

        v1.show_diff(&v2);
    }
    #[test]
    fn vector_min_max_sorted() {
        // dense
        let mut v: Vector<DenseStorage<Prob>, usize> =
            Vector::from_slice(&[p(0.3), p(0.0), p(1.0), p(0.8), p(0.0)], p(0.0));
        assert_eq!(v.max(), Some((2, p(1.0))));
        assert_eq!(v.min(), Some((1, p(0.0))));
        let sorted: Vec<_> = v.iter_sorted_desc().collect();
        assert_eq!(
            sorted,
            vec![
                (2, p(1.0)),
                (3, p(0.8)),
                (0, p(0.3)),
                (1, p(0.0)),
                (4, p(0.0))
            ]
        );

        // sparse
        let mut v: Vector<SparseStorage<Prob>, usize> =
            Vector::from_slice(&[p(0.3), p(0.0), p(1.0), p(0.8), p(0.0)], p(0.0));
        assert_eq!(v.max(), Some((2, p(1.0))));
        assert_eq!(v.min(), Some((1, p(0.0))));
        let sorted: Vec<_> = v.iter_sorted_desc().collect();
        assert_eq!(
            sorted,
            vec![
                (2, p(1.0)),
                (3, p(0.8)),
                (0, p(0.3)),
                (1, p(0.0)),
                (4, p(0.0))
            ]
        );
    }
    #[test]
    fn vector_compare() {
        let mut v: Vector<SparseStorage<Prob>, usize> =
            Vector::from_slice(&[p(0.3), p(0.011), p(1.0), p(0.8), p(0.0), p(0.11)], p(0.0));
        let mut w: Vector<DenseStorage<Prob>, usize> =
            Vector::from_slice(&[p(0.3), p(0.022), p(1.0), p(0.8), p(0.0), p(0.01)], p(0.0));
        assert!(v.is_same_ranking(&w, p(0.2)));
        assert!(!v.is_same_ranking(&w, p(0.0)));

        let mut w: Vector<SparseStorage<Prob>, usize> =
            Vector::from_slice(&[p(0.3), p(0.022), p(1.0), p(0.2), p(0.0), p(0.12)], p(0.0));
        assert!(!v.is_same_ranking(&w, p(0.2)));
    }
}
