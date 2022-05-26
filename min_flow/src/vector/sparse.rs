//!
//! Sparse storage that uses `ArrayVec`
//!
use super::{DenseStorage, Storage};
use arrayvec::ArrayVec;

/// SparseStorage max index size parameter
pub const SIZE: usize = 400;

/// Sparse storage powered by `ArrayVec`
/// A type of items is represented as `T`.
///
/// This sparse storage virtually has the given size
/// `0 <= index < size`
///
/// Hypothesis is that it has at most `SIZE` elements
/// which has a non-default value.
/// `0 <= id < SIZE = 200`
///
#[derive(Debug, Clone)]
pub struct SparseStorage<T> {
    /// virtual size of this storage
    size: usize,
    /// default value
    /// The value of unused index wiil be filled with
    /// this `default_value`.
    default_value: T,
    /// ArrayVec of elements `(index, value)`
    elements: ArrayVec<(usize, T), SIZE>,
}

impl<T> Storage for SparseStorage<T>
where
    T: Copy + PartialEq + Default,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> SparseStorage<T> {
        SparseStorage {
            size,
            default_value,
            elements: ArrayVec::<(usize, T), SIZE>::new(),
        }
    }
    #[inline]
    fn size(&self) -> usize {
        self.size
    }
    #[inline]
    fn n_ids(&self) -> usize {
        self.elements.len()
    }
    #[inline]
    fn get_by_id(&self, id: usize) -> (usize, T) {
        assert!(id < self.n_ids());
        self.elements[id]
    }
    fn get(&self, index: usize) -> &T {
        assert!(index < self.size());

        // search for existing entries
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return &self.elements[i].1;
            }
        }

        // not found
        &self.default_value
    }
    fn get_mut(&mut self, index: usize) -> &mut T {
        assert!(index < self.size());

        // search for existing entries
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return &mut self.elements[i].1;
            }
        }

        // add a new entry and return the reference to it
        self.elements.push((index, self.default_value));
        let n = self.elements.len();
        return &mut self.elements[n - 1].1;
    }
    fn set(&mut self, index: usize, value: T) {
        assert!(index < self.size());
        if value == self.default_value {
            // delete existing element if exists
            for i in 0..self.elements.len() {
                if self.elements[i].0 == index {
                    self.elements.remove(i);
                    break;
                }
            }
        } else {
            // store the value using get_mut
            *self.get_mut(index) = value;
        }
    }
    fn has(&self, index: usize) -> bool {
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return true;
            }
        }
        false
    }
    fn try_get(&self, index: usize) -> Option<&T> {
        for i in 0..self.elements.len() {
            if self.elements[i].0 == index {
                return Some(&self.elements[i].1);
            }
        }
        None
    }
    fn mutate<F: FnMut(usize, &mut T)>(&mut self, mut f: F) {
        let default_value = self.default_value;
        self.elements.retain(|(i, x)| {
            f(*i, x);
            // retain if x is not default_value
            *x != default_value
        })
    }
    fn default_value(&self) -> T {
        self.default_value
    }
    fn set_default_value(&mut self, default_value: T) {
        self.default_value = default_value;
    }
    fn to_dense(&self) -> DenseStorage<T> {
        let mut s: DenseStorage<T> = DenseStorage::new(self.size(), self.default_value);
        for (index, value) in self.iter() {
            *s.get_mut(index) = value;
        }
        s
    }
    fn to_sparse(&self, _: T) -> Self {
        self.clone()
    }
    fn is_dense() -> bool {
        false
    }
}

// Custom Partial Eq for Sparse storage
impl<T> PartialEq for SparseStorage<T>
where
    T: Copy + PartialEq + Default,
{
    fn eq(&self, other: &Self) -> bool {
        self.size == other.size
            && self.default_value == other.default_value
            && self.iter().all(|(i, x)| *other.get(i) == x)
            && other.iter().all(|(i, x)| *self.get(i) == x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::Vector;
    use super::*;

    #[test]
    fn sparse_storage() {
        // u32
        let mut s: SparseStorage<u32> = SparseStorage::new(10, 0);
        *s.get_mut(5) = 111;
        *s.get_mut(3) = 22;
        assert_eq!(*s.get(0), 0);
        assert_eq!(*s.get(3), 22);
        assert_eq!(*s.get(5), 111);
        assert_eq!(s.size(), 10);
        assert_eq!(s.n_ids(), 2);
        assert_eq!(s.get_by_id(0), (5, 111));
        assert_eq!(s.get_by_id(1), (3, 22));

        // clone
        let mut s2 = s.clone();
        assert_eq!(*s2.get(3), 22);
        *s2.get_mut(3) = 21;
        assert_eq!(*s2.get(3), 21);
        assert_eq!(*s.get(3), 22);

        // f64
        let mut s: SparseStorage<f64> = SparseStorage::new(10, 0.0);
        *s.get_mut(5) = 12.11;
        *s.get_mut(3) = 10.0;
        assert_eq!(*s.get(0), 0.0);
        assert_eq!(*s.get(3), 10.0);
        assert_eq!(*s.get(5), 12.11);
    }
    #[test]
    #[should_panic]
    fn sparse_storage_outside() {
        let mut s: SparseStorage<u32> = SparseStorage::new(3, 0);
        *s.get_mut(3) = 22;
    }
    #[test]
    #[should_panic]
    fn sparse_storage_outside_id() {
        let s: SparseStorage<u32> = SparseStorage::new(3, 0);
        s.get_by_id(3);
    }
    #[test]
    fn sparse_storage_iter() {
        let mut s: SparseStorage<u32> = SparseStorage::new(4, 0);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        let v: Vec<(usize, u32)> = s.iter().collect();
        assert_eq!(v, vec![(0, 111), (2, 10)]);
    }
    #[test]
    fn sparse_storage_conversion() {
        let mut s: SparseStorage<u32> = SparseStorage::new(4, 5);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        println!("{:?}", s);
        let s2 = s.to_dense();
        println!("{:?}", s2);
        assert_eq!(*s2.get(0), *s.get(0));
        assert_eq!(*s2.get(1), *s.get(1));
        assert_eq!(*s2.get(2), *s.get(2));
        assert_eq!(*s2.get(3), *s.get(3));
    }
    #[test]
    fn sparse_storage_mutate() {
        let mut s: SparseStorage<u32> = SparseStorage::new(4, 0);
        *s.get_mut(0) = 111;
        *s.get_mut(1) = 222;
        *s.get_mut(2) = 10;
        println!("{:?}", s);
        s.mutate(|i, x| {
            if i % 2 == 0 {
                *x = *x + 10;
            } else {
                *x = 0;
            }
        });
        println!("{:?}", s);
        assert_eq!(*s.get(0), 121);
        assert_eq!(*s.get(1), 0);
        assert_eq!(*s.get(2), 20);
        assert_eq!(*s.get(3), 0);
    }
    #[test]
    fn sparse_storage_dense_check() {
        assert!(!SparseStorage::<u32>::is_dense());
    }
    #[test]
    fn sparse_storage_partial_eq() {
        let mut s1: SparseStorage<u32> = SparseStorage::new(10, 0);
        *s1.get_mut(5) = 111;
        *s1.get_mut(3) = 22;
        let mut s2: SparseStorage<u32> = SparseStorage::new(10, 0);
        *s2.get_mut(3) = 22;
        *s2.get_mut(5) = 111;
        println!("{:?}", s1);
        println!("{:?}", s2);
        assert!(s1 == s2);
        let mut s3: SparseStorage<u32> = SparseStorage::new(10, 0);
        *s2.get_mut(3) = 22;
        *s2.get_mut(2) = 33;
        println!("{}", s1 == s3);
        assert!(!(s1 == s3));
    }
    #[test]
    fn sparse_storage_vector() {
        let mut v: Vector<SparseStorage<u32>> = Vector::new(5, 0);
        v[0] = 100;
        v[3] = 222;
        assert_eq!(v[0], 100);
        assert_eq!(v[1], 0);
        let w: Vec<(usize, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(w, vec![(0, 100), (3, 222)]);
    }
    #[test]
    fn sparse_storage_vector_add_mul() {
        let mut v1: Vector<SparseStorage<u32>> = Vector::new(4, 0);
        v1[0] = 120;
        v1[3] = 111;
        let mut v2: Vector<SparseStorage<u32>> = Vector::new(4, 0);
        v2[0] = 1;
        v2[2] = 111;
        v2[3] = 1;
        println!("{:?}", v1);
        println!("{:?}", v2);
        // v1 + v2
        let added = &v1 + &v2;
        let muled = &v1 * &v2;
        println!("{:?}", added);
        assert_eq!(added[0], 120 + 1);
        assert_eq!(added[1], 0 + 0);
        assert_eq!(added[2], 0 + 111);
        assert_eq!(added[3], 111 + 1);
        println!("{:?}", muled);
        assert_eq!(muled[0], 120 * 1);
        assert_eq!(muled[1], 0 * 0);
        assert_eq!(muled[2], 0 * 111);
        assert_eq!(muled[3], 111 * 1);
    }
    #[test]
    fn arrayvec_retain_test() {
        let b: Vec<usize> = vec![0, 4, 8, 12, 16];
        // (1) arrayvec
        let mut a = ArrayVec::<usize, 10>::new();
        for i in 0..10 {
            a.push(i);
        }
        a.retain(|x| {
            if *x % 2 == 0 {
                *x *= 2;
                true
            } else {
                false
            }
        });
        let av: Vec<usize> = a.into_iter().collect();
        assert_eq!(av, b);

        // (2) std::vec
        let mut v = Vec::new();
        for i in 0..10 {
            v.push(i);
        }
        v.retain_mut(|x| {
            if *x % 2 == 0 {
                *x *= 2;
                true
            } else {
                false
            }
        });
        assert_eq!(v, b);
    }
}
