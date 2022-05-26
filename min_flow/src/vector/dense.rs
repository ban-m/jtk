//!
//! Dense storage that uses `std::Vec`
//!
use super::{SparseStorage, Storage};

/// Dense storage powered by `std::Vec`
///
/// In `DenseStorage`, internal id equals to index.
#[derive(PartialEq, Debug, Clone)]
pub struct DenseStorage<T>(Vec<T>);

impl<T> Storage for DenseStorage<T>
where
    T: Copy + PartialEq + Default,
{
    type Item = T;
    fn new(size: usize, default_value: T) -> DenseStorage<T> {
        DenseStorage(vec![default_value; size])
    }
    #[inline]
    fn size(&self) -> usize {
        self.0.len()
    }
    #[inline]
    fn n_ids(&self) -> usize {
        self.0.len()
    }
    #[inline]
    fn get_by_id(&self, id: usize) -> (usize, T) {
        let index = id;
        let value = self.0[id];
        (index, value)
    }
    #[inline]
    fn get(&self, index: usize) -> &T {
        &self.0[index]
    }
    #[inline]
    fn get_mut(&mut self, index: usize) -> &mut T {
        &mut self.0[index]
    }
    #[inline]
    fn set(&mut self, index: usize, value: T) {
        self.0[index] = value;
    }
    #[inline]
    fn has(&self, _: usize) -> bool {
        // dense storage has all element separatedly
        true
    }
    #[inline]
    fn try_get(&self, index: usize) -> Option<&T> {
        // dense storage has all element separatedly
        Some(&self.0[index])
    }
    fn mutate<F: FnMut(usize, &mut T)>(&mut self, mut f: F) {
        let mut i = 0;
        self.0.retain_mut(|x| {
            // do some modification
            f(i, x);
            i += 1;
            true
        })
    }
    #[inline]
    fn default_value(&self) -> T {
        // DenseStorage do not have default values, so the returned value has no meaning
        T::default()
    }
    #[inline]
    fn set_default_value(&mut self, _: T) {
        // do nothing, because DenseStorage do not have default values
    }
    #[inline]
    fn to_dense(&self) -> Self {
        self.clone()
    }
    fn to_sparse(&self, default_value: T) -> SparseStorage<T> {
        // TODO use the most frequent values as the default_value
        let mut s: SparseStorage<T> = SparseStorage::new(self.size(), default_value);
        for (index, value) in self.iter() {
            if value != default_value {
                *s.get_mut(index) = value;
            }
        }
        s
    }
    // TODO add try_to_sparse
    fn is_dense() -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::super::Vector;
    use super::*;

    #[test]
    fn dense_storage() {
        // u32
        let mut s: DenseStorage<u32> = DenseStorage::new(10, 0);
        *s.get_mut(5) = 111;
        *s.get_mut(3) = 22;
        assert_eq!(*s.get(0), 0);
        assert_eq!(*s.get(3), 22);
        assert_eq!(*s.get(5), 111);
        assert_eq!(s.size(), 10);
        assert_eq!(s.n_ids(), 10);
        assert_eq!(s.get_by_id(3), (3, 22));

        // clone
        let mut s2 = s.clone();
        assert_eq!(*s2.get(3), 22);
        *s2.get_mut(3) = 21;
        assert_eq!(*s2.get(3), 21);
        assert_eq!(*s.get(3), 22);

        // f64
        let mut s: DenseStorage<f64> = DenseStorage::new(10, 0.0);
        *s.get_mut(5) = 12.11;
        *s.get_mut(3) = 10.0;
        assert_eq!(*s.get(0), 0.0);
        assert_eq!(*s.get(3), 10.0);
        assert_eq!(*s.get(5), 12.11);
    }
    #[test]
    #[should_panic]
    fn dense_storage_outside() {
        let mut s: DenseStorage<u32> = DenseStorage::new(3, 0);
        *s.get_mut(3) = 22;
    }
    #[test]
    fn dense_storage_iter() {
        let mut s: DenseStorage<u32> = DenseStorage::new(4, 5);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        let v: Vec<(usize, u32)> = s.iter().collect();
        assert_eq!(v, vec![(0, 111), (1, 5), (2, 10), (3, 5)]);
    }
    #[test]
    fn dense_storage_conversion() {
        // non-zero default values
        let mut s: DenseStorage<u32> = DenseStorage::new(4, 5);
        *s.get_mut(0) = 111;
        *s.get_mut(2) = 10;
        println!("{:?}", s);
        let s2 = s.to_sparse(5);
        println!("{:?}", s2);
        assert_eq!(*s2.get(0), *s.get(0));
        assert_eq!(*s2.get(1), *s.get(1));
        assert_eq!(*s2.get(2), *s.get(2));
        assert_eq!(*s2.get(3), *s.get(3));
    }
    #[test]
    fn dense_storage_mutate() {
        let mut s: DenseStorage<u32> = DenseStorage::new(4, 5);
        *s.get_mut(0) = 111;
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
    fn dense_storage_dense_check() {
        assert!(DenseStorage::<u32>::is_dense());
    }
    #[test]
    fn dense_storage_partial_eq() {
        let mut s1: DenseStorage<u32> = DenseStorage::new(10, 0);
        *s1.get_mut(5) = 111;
        *s1.get_mut(3) = 22;
        let mut s2: DenseStorage<u32> = DenseStorage::new(10, 0);
        *s2.get_mut(3) = 22;
        *s2.get_mut(5) = 111;
        println!("{:?}", s1);
        println!("{:?}", s2);
        assert!(s1 == s2);
    }
}
