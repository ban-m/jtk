//! Principal component analysis. Specifically, it reduces the dimension into K.
use nalgebra;
use nalgebra::{DMatrix, DVector};
pub fn pca(data: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let len = data[0].len();
    if len < k {
        panic!(
            "This vector's dimension({}) is already less than {}.",
            len, k
        );
    }
    data.iter().for_each(|xs| assert_eq!(xs.len(), len));
    let data: Vec<DVector<_>> = data
        .iter()
        .map(|xs| DVector::from_column_slice(xs.as_slice()))
        .collect();
    let mean = data.iter().fold(DVector::zeros(len), |x, y| x + y) / data.len() as f64;
    let covariance = data
        .iter()
        .map(|x| (x - mean.clone()) * (x - mean.clone()).transpose())
        .fold(DMatrix::zeros(len, len), |x, y| x + y)
        / data.len() as f64;
    let eigens = covariance.clone().symmetric_eigen();
    let mut eigen_and_eigenvec: Vec<_> = eigens
        .eigenvectors
        .column_iter()
        .zip(eigens.eigenvalues.iter())
        .collect();
    eigen_and_eigenvec.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
    eigen_and_eigenvec.reverse();
    for (_, val) in eigen_and_eigenvec.iter() {
        eprintln!("{}", val);
    }
    let pca_vectors = &eigen_and_eigenvec[..k];
    data.iter()
        .map(|x| pca_vectors.iter().map(|(v, _)| x.dot(&v)).collect())
        .collect()
}
