// //!
// //! Vector benchmark
// //!

// #[cfg(test)]
// mod tests {
//     use super::super::{DenseStorage, SparseStorage, Storage, Vector};
//     use test::Bencher;
//     #[bench]
//     fn bench_dense_vector_add(b: &mut Bencher) {
//         let n = 100000;
//         let mut v: Vector<DenseStorage<u32>> = Vector::new(n, 0);
//         let mut w: Vector<DenseStorage<u32>> = Vector::new(n, 0);
//         b.iter(|| {
//             let u = &v + &w;
//         });
//     }
//     #[bench]
//     fn bench_sparse_vector_add(b: &mut Bencher) {
//         let n = 100000;
//         let mut v: Vector<SparseStorage<u32>> = Vector::new(n, 0);
//         let mut w: Vector<SparseStorage<u32>> = Vector::new(n, 0);
//         b.iter(|| {
//             let u = &v + &w;
//         });
//     }
// }
