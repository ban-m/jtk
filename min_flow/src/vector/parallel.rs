//!
//! Vector parallel summation tests
//!

// #[cfg(test)]
// mod tests {
//     use crate::vector::{DenseStorage, SparseStorage, Storage, Vector};
//     use rayon::prelude::*;

//     #[test]
//     fn vector_parallel_sum() {
//         let inputs: Vec<usize> = vec![1; 100];
//         let mapped = inputs.par_iter().map(|&input| {
//             let v: Vector<DenseStorage<usize>> = Vector::new(1000, input);
//             v
//         });

//         // XXX
//         // this causes exception (but it compiles)
//         // let sum: Vector<DenseStorage<usize>> = mapped.sum()

//         // use this method
//         let sum: Vector<DenseStorage<usize>> = mapped.reduce(
//             || Vector::new(1000, 0),
//             |mut a, b| {
//                 a += &b;
//                 a
//             },
//         );
//         println!("{}", sum);
//     }
// }
