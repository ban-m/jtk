//!
//! `Vector` Wrapper of fixed size table
//!
//!
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign};
pub mod bench;
pub mod dense;
pub mod diff;
pub mod graph;
pub mod index;
pub mod parallel;
pub mod sparse;
pub mod test;
pub mod unit;
pub use dense::DenseStorage;
pub use graph::{EdgeVec, NodeVec};
pub use index::Indexable;
pub use sparse::SparseStorage;
use std::marker::PhantomData;
use unit::{UnitAdd, UnitMul};

/// Backend storage of `Vector`
/// an abstruction of a vec with fixed size that is readable/writable
/// by index.
///
/// * `new`
///     create storage with fixed size and filled with the default value
/// * `size`
///     get the fixed size
/// * `get`
///     get the reference to the value in the index
/// * `get_mut`
///     get the mutable reference to the value in the index
/// * `iter`
///     get an iterator of (index, value) whose value is not default.
///
/// Internaly, mapping of index (0 <= index < self.size()) to value
/// is stored with internal element id. (0 <= id < self.n_ids())
/// (index) -> (id) -> (value)
///
pub trait Storage: Clone + Sized + PartialEq {
    /// Item type that this storage stores.
    ///
    type Item: Copy + PartialEq + Default;
    ///
    /// Create a new storage with fixed size and filled with the default value
    fn new(size: usize, default_value: Self::Item) -> Self;
    ///
    /// Get the size of this storage
    fn size(&self) -> usize;
    ///
    /// Get the number of elements.
    /// `0 <= id < self.n_ids()` should be satisfied.
    fn n_ids(&self) -> usize;
    ///
    /// Get the reference to the value at the internal id
    fn get_by_id(&self, id: usize) -> (usize, Self::Item);
    ///
    /// Get the reference to the value at the given index
    fn get(&self, index: usize) -> &Self::Item;
    ///
    /// Get the mutable reference to the given index
    fn get_mut(&mut self, index: usize) -> &mut Self::Item;
    ///
    /// Store the value to the given index
    fn set(&mut self, index: usize, value: Self::Item);
    ///
    /// check if this storage has an element for the given index?
    fn has(&self, index: usize) -> bool;
    ///
    /// if the value of the given index is stored in the storage,
    /// then return Some(value)
    fn try_get(&self, index: usize) -> Option<&Self::Item>;
    ///
    /// mutate all elements
    fn mutate<F: FnMut(usize, &mut Self::Item)>(&mut self, f: F);
    ///
    /// Get the current default value
    fn default_value(&self) -> Self::Item;
    ///
    /// set the current default value
    fn set_default_value(&mut self, default_value: Self::Item);
    ///
    /// get an iterator of (usize, Self::Item) on the storage
    fn iter<'a>(&'a self) -> StorageIterator<'a, Self> {
        StorageIterator {
            id: 0,
            storage: self,
        }
    }
    ///
    /// Convert to the DenseStorage with same contents
    fn to_dense(&self) -> DenseStorage<Self::Item>;
    ///
    /// Convert to the SparseStorage with same contents
    /// with specifying `default_value` in SparseStorage.
    fn to_sparse(&self, default_value: Self::Item) -> SparseStorage<Self::Item>;
    ///
    /// Check if this is dense storage or not
    fn is_dense() -> bool;
}

///
/// Iterator struct of Storage
pub struct StorageIterator<'a, S: Storage> {
    /// current internal id
    id: usize,
    /// reference of the storage
    storage: &'a S,
}

impl<'a, S: Storage> Iterator for StorageIterator<'a, S> {
    type Item = (usize, S::Item);
    fn next(&mut self) -> Option<Self::Item> {
        if self.id < self.storage.n_ids() {
            let id = self.id;
            let item = self.storage.get_by_id(id);
            self.id += 1;
            Some(item)
        } else {
            None
        }
    }
}

/// `Vector` struct
///
/// It generalized of
///
/// 1. item type `Storage::Item`
/// 2. backend storage `S: Storage`
/// 3. index type `Ix: Indexable`
///
/// TODO
/// PartialEq and Eq should be revised.
/// For sparsestorage, PartialEq should be ignore the ordering of elements.
///
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct Vector<S: Storage, Ix: Indexable = usize> {
    /// Backend storage of the Vector
    storage: S,
    /// Hidden marker of index type
    ty: PhantomData<Ix>,
}

impl<S: Storage, Ix: Indexable> Vector<S, Ix> {
    /// Create a new Vector, with fixed size and filled by default_value.
    pub fn new(size: usize, default_value: S::Item) -> Vector<S, Ix> {
        Vector {
            storage: S::new(size, default_value),
            ty: PhantomData,
        }
    }
    /// Create a new Vector, from the elements
    pub fn from_vec(size: usize, default_value: S::Item, vec: &[(Ix, S::Item)]) -> Vector<S, Ix> {
        let mut v = Vector::new(size, default_value);
        for (index, value) in vec.iter() {
            v[*index] = *value;
        }
        v
    }
    /// Create a new Vector, from slice
    pub fn from_slice(vec: &[S::Item], default_value: S::Item) -> Vector<S, Ix> {
        let mut v = Vector::new(vec.len(), default_value);
        for (index, value) in vec.iter().enumerate() {
            if *value != default_value {
                v[Ix::new(index)] = *value;
            }
        }
        v
    }
    /// Get an (virtual) size of the storage
    pub fn len(&self) -> usize {
        self.storage.size()
    }
    /// Get an (sparse) iterator on (index, item).
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = (Ix, S::Item)> {
        self.storage.iter().map(|(i, v)| (Ix::new(i), v))
    }
    /// Get an exhaustive iterator on (index, item).
    pub fn iter_all<'a>(&'a self) -> impl 'a + Iterator<Item = (Ix, S::Item)> {
        (0..self.len()).map(move |i| (Ix::new(i), *self.storage.get(i)))
    }
    /// Convert to the DenseStorage-backed vector.
    pub fn to_dense(&self) -> Vector<DenseStorage<S::Item>, Ix> {
        Vector {
            storage: self.storage.to_dense(),
            ty: PhantomData,
        }
    }
    /// Convert to the SparseStorage-backed vector.
    pub fn to_sparse(&self, default_value: S::Item) -> Vector<SparseStorage<S::Item>, Ix> {
        Vector {
            storage: self.storage.to_sparse(default_value),
            ty: PhantomData,
        }
    }
    /// Convert to the SparseStorage-backed vector
    /// storing values of the specified indexes.
    pub fn to_sparse_by_indexes(
        &self,
        default_value: S::Item,
        indexes: &[Ix],
    ) -> Vector<SparseStorage<S::Item>, Ix> {
        let mut v = Vector {
            storage: SparseStorage::new(self.storage.size(), default_value),
            ty: PhantomData,
        };
        for &index in indexes.iter() {
            v[index] = self[index];
        }
        v
    }
    /// internal storage is dense or sparse?
    ///
    pub fn is_dense(&self) -> bool {
        S::is_dense()
    }
    /// Convert to normal vector
    pub fn to_vec(&self) -> Vec<S::Item> {
        (0..self.len()).map(|i| self[Ix::new(i)]).collect()
    }
    /// Change index
    pub fn switch_index<Jx: Indexable>(self) -> Vector<S, Jx> {
        Vector {
            storage: self.storage,
            ty: PhantomData,
        }
    }
}

/// private associated functions
/// to use math ops between two vector
impl<S: Storage, Ix: Indexable> Vector<S, Ix> {}

impl<S, Ix> Vector<S, Ix>
where
    S: Storage,
    S::Item: std::iter::Sum,
    Ix: Indexable,
{
    /// Calculate the sum of the elements
    ///
    /// TODO
    /// This function can be slow when with sparse storage
    pub fn sum(&self) -> S::Item {
        (0..self.len()).map(|i| self[Ix::new(i)]).sum()
    }
}

/// Implement index access, vec[i]
impl<S: Storage, Ix: Indexable> Index<Ix> for Vector<S, Ix> {
    type Output = S::Item;
    fn index(&self, index: Ix) -> &Self::Output {
        self.storage.get(index.index())
    }
}

/// Implement index write access, vec[i] = 10
impl<S: Storage, Ix: Indexable> IndexMut<Ix> for Vector<S, Ix> {
    fn index_mut(&mut self, index: Ix) -> &mut Self::Output {
        self.storage.get_mut(index.index())
    }
}

/// Implement addition `+` between two vecs
/// if the item of vec supports addition
impl<'a, 'b, Sa, Sb, Ix> Add<&'a Vector<Sa, Ix>> for &'b Vector<Sb, Ix>
where
    Sa: Storage,
    Sb: Storage<Item = Sa::Item>,
    Sa::Item: Add<Output = Sa::Item>,
    Ix: Indexable,
{
    type Output = Vector<Sb, Ix>;
    fn add(self, other: &'a Vector<Sa, Ix>) -> Self::Output {
        assert_eq!(self.len(), other.len());
        let mut ret: Vector<Sb, Ix> = Vector::new(
            self.len(),
            self.storage.default_value() + other.storage.default_value(),
        );
        // fill for the used indexes of self
        for (index, value) in self.iter() {
            ret[index] = other[index] + value;
        }
        // fill for the used indexes of other
        // if self is dense, both loop travarses the all indexes
        if !self.is_dense() {
            for (index, value) in other.iter() {
                ret[index] = self[index] + value;
            }
        }
        ret
    }
}

/// add constant to vector
/// `Vector<S> + S::Item = Vector<S>`
///
impl<'a, 'b, S, Ix> Add<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Add<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn add(mut self, other: S::Item) -> Self::Output {
        // add to inactive indexes
        self.storage
            .set_default_value(self.storage.default_value() + other);

        // add other to active indexes of self
        self.storage.mutate(|_, x| *x = *x + other);
        self
    }
}

/// Implement addition with assignment `+=` between two vecs
/// if the item of vec supports addition
/// This does not cause re-allocation
///
/// `A += B`
///
/// * (A: Dense, B: Dense) -> Dense
///     travarsing all index
/// * (A: Dense, B: Sparse) -> Dense
///     * if B.default_value.is_unit_add()
///         only active elements in B should be considered
///     * otherwise
///         travarsing all index
/// * (A: Sparse, B: Dense) -> Sparse
///     first convert B into sparse
///     then fallback to (A: Sparse, B: Sparse)
/// * (A: Sparse, B: Sparse) -> Sparse
///     merge both active indexes
///
impl<'a, Sa, S, Ix> AddAssign<&'a Vector<Sa, Ix>> for Vector<S, Ix>
where
    S: Storage,
    Sa: Storage<Item = S::Item>,
    S::Item: Add<Output = S::Item> + UnitAdd,
    Ix: Indexable,
{
    fn add_assign(&mut self, other: &'a Vector<Sa, Ix>) {
        assert_eq!(self.len(), other.len());

        match (
            self.is_dense(),
            other.is_dense(),
            other.storage.default_value().is_unit_add(),
        ) {
            // self is dense
            (true, false, true) => {
                // other is sparse with unit default value
                // then it can only add only active indexes in other
                for (index, value) in other.iter() {
                    self[index] = self[index] + value;
                }
            }
            (true, _, _) => {
                // otherwise, all values in self will be modified
                self.storage
                    .mutate(|i, x| *x = *x + (*other.storage.get(i)));
            }
            (false, true, _) => {
                // first convert other into sparse with unit default_value
                let other_sparse = other.to_sparse(S::Item::unit_add());
                // fallback to (sparse, sparse) add_assign
                self.add_assign(&other_sparse);
            }
            (false, false, _) => {
                // both is sparse
                let default_value = other.storage.default_value();
                // add the active indexes of other into self
                for (index, value) in other.iter() {
                    self[index] = self[index] + value;
                }
                // add to inactive indexes
                // this will affect the index-access to non-stored element in self.
                self.storage
                    .set_default_value(self.storage.default_value() + default_value);
                // add the default_value of other into the self-only elements
                self.storage.mutate(|i, x| {
                    if !other.storage.has(i) {
                        *x = *x + default_value;
                    }
                });
            }
        }
    }
}

/// Implement multiplication `*` between two vecs
/// if the item of vec supports multiplication
impl<'a, 'b, Sa, Sb, Ix> Mul<&'a Vector<Sa, Ix>> for &'b Vector<Sb, Ix>
where
    Sa: Storage,
    Sb: Storage<Item = Sa::Item>,
    Sa::Item: Mul<Output = Sa::Item>,
    Ix: Indexable,
{
    type Output = Vector<Sb, Ix>;
    fn mul(self, other: &'a Vector<Sa, Ix>) -> Self::Output {
        assert_eq!(self.len(), other.len());
        let mut ret: Vector<Sb, Ix> = Vector::new(
            self.len(),
            self.storage.default_value() * other.storage.default_value(),
        );
        // fill for the used indexes of self
        for (index, value) in self.iter() {
            ret[index] = other[index] * value;
        }
        if !self.is_dense() {
            // fill for the used indexes of other
            for (index, value) in other.iter() {
                ret[index] = self[index] * value;
            }
        }
        ret
    }
}

/// Implement multiplication with assignment `*=` between two vecs
/// if the item of vec supports multiplication
/// This does not cause re-allocation
impl<'a, Sa, S, Ix> MulAssign<&'a Vector<Sa, Ix>> for Vector<S, Ix>
where
    S: Storage,
    Sa: Storage<Item = S::Item>,
    S::Item: Mul<Output = S::Item> + UnitMul,
    Ix: Indexable,
{
    fn mul_assign(&mut self, other: &'a Vector<Sa, Ix>) {
        assert_eq!(self.len(), other.len());
        if self.is_dense() {
            // self is dense
            if !other.is_dense() && other.storage.default_value().is_unit_mul() {
                // other is sparse with unit default value
                // then it can only add only active indexes in other
                println!("fast");
                for (index, value) in other.iter() {
                    self[index] = self[index] * value;
                }
            } else {
                // otherwise, all values in self will be modified
                self.storage
                    .mutate(|i, x| *x = *x * (*other.storage.get(i)));
            }
        } else {
            // self is sparse
            if other.is_dense() {
                // if other is dense, first convert other into sparse with unit default_value
                // TODO should use zero_mul?
                // count zero/unit and use it as a default_value
                let other_sparse = other.to_sparse(S::Item::zero_mul());
                // fallback to (sparse, sparse) add_assign
                self.mul_assign(&other_sparse);
            } else {
                // both is sparse
                if self.storage.default_value().is_zero_mul() {
                    println!("sparse(zero) x sparse");
                    // only active nodes of self will be remains.
                    self.storage
                        .mutate(|i, x| *x = *x * (*other.storage.get(i)));
                } else if other.storage.default_value().is_zero_mul() {
                    println!("sparse x sparse(zero)");
                    // purge active nodes of self, which is not active in other (=zero)
                    self.storage.mutate(|i, x| {
                        if !other.storage.has(i) {
                            *x = S::Item::zero_mul();
                        }
                    });
                    // set for active nodes of other.
                    for id in 0..other.storage.n_ids() {
                        let (index, other_value) = other.storage.get_by_id(id);
                        *self.storage.get_mut(index) = other_value * *self.storage.get(index);
                    }
                    // only active nodes of other will be remains.
                    self.storage.set_default_value(S::Item::zero_mul());
                } else {
                    println!("sparse x sparse");
                    let default_value = other.storage.default_value();
                    // add the active indexes of other into self
                    for (index, value) in other.iter() {
                        self[index] = self[index] * value;
                    }
                    // add to inactive indexes
                    // this will affect the index-access to non-stored element in self.
                    self.storage
                        .set_default_value(self.storage.default_value() * default_value);
                    // add the default_value of other into the self-only elements
                    self.storage.mutate(|i, x| {
                        if !other.storage.has(i) {
                            *x = *x * default_value;
                        }
                    });
                }
            }
        }
    }
}

/// multiply a constant to vector
/// `Vector<S> * S::Item = Vector<S>`
///
impl<'a, 'b, S, Ix> Mul<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Mul<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn mul(mut self, other: S::Item) -> Self::Output {
        // add to inactive indexes
        self.storage
            .set_default_value(self.storage.default_value() * other);

        // add other to active indexes of self
        self.storage.mutate(|_, x| *x = *x * other);
        self
    }
}

/// divide-by a constant to vector
/// `Vector<S> / S::Item = Vector<S>`
///
impl<'a, 'b, S, Ix> Div<S::Item> for Vector<S, Ix>
where
    S: Storage,
    S::Item: Div<Output = S::Item> + Copy,
    Ix: Indexable,
{
    type Output = Vector<S, Ix>;
    fn div(mut self, other: S::Item) -> Self::Output {
        // add to inactive indexes
        self.storage
            .set_default_value(self.storage.default_value() / other);

        // add other to active indexes of self
        self.storage.mutate(|_, x| *x = *x / other);
        self
    }
}

/// sum of vectors
///
/// ## TODO
///
/// * should implement sum of zero element vector (to avoid errors in rayon calculation).
/// but the size is unknown...
///
impl<S, Ix> std::iter::Sum for Vector<S, Ix>
where
    S: Storage,
    S::Item: Add<Output = S::Item> + UnitAdd + Copy,
    Ix: Indexable,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|mut a, b| {
            a += &b;
            a
        })
        .unwrap()
    }
}

//
// Display
//
impl<S, Ix> std::fmt::Display for Vector<S, Ix>
where
    S: Storage,
    S::Item: std::fmt::Display,
    Ix: Indexable,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        for i in 0..self.len() {
            write!(f, "{}, ", self[Ix::new(i)])?;
        }
        write!(f, "]")?;
        Ok(())
    }
}

/// for approx `assert_abs_diff_eq`
use approx::AbsDiffEq;
impl<S, Ix> AbsDiffEq for Vector<S, Ix>
where
    S: Storage,
    S::Item: AbsDiffEq,
    <<S as Storage>::Item as AbsDiffEq>::Epsilon: Copy,
    Ix: Indexable,
{
    type Epsilon = <<S as Storage>::Item as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        S::Item::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        /*
         * TODO this may cause incorrect result
        let p1 = self
            .iter()
            .all(|(index, value)| S::Item::abs_diff_eq(&other[index], &value, epsilon));
        let p2 = other
            .iter()
            .all(|(index, value)| S::Item::abs_diff_eq(&self[index], &value, epsilon));
        p1 && p2
        */
        (0..self.len())
            .all(|i| S::Item::abs_diff_eq(&self[Ix::new(i)], &other[Ix::new(i)], epsilon))
    }
}
