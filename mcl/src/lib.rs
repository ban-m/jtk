//! MCL algorithm implementation.
use log::*;
use nalgebra::DMatrix;
// use serde::{Deserialize, Serialize};
#[derive(Debug, Clone)]
pub struct Graph {
    matrix: DMatrix<f64>,
}

impl Graph {
    /// Take the number of nodes and the adjacency list,
    /// and return the corresponding (normalized) graph.
    /// Note that the weight should be recorded with the edge.
    pub fn new(nodes: usize, edges: &[Vec<(usize, u64)>]) -> Self {
        // Adj matrix
        let mut raw_vector = vec![0.; nodes * nodes];
        for (from, edge) in edges.iter().enumerate() {
            for &(to, weight) in edge {
                raw_vector[from * nodes + to] += weight as f64;
            }
        }
        // Add self loop.
        for i in 0..nodes {
            raw_vector[i * nodes + i] += 1.;
        }
        // Normalize.
        raw_vector.chunks_exact_mut(nodes).for_each(|column| {
            let sum = column.iter().sum::<f64>();
            column.iter_mut().for_each(|x| *x /= sum);
        });
        Self {
            matrix: DMatrix::from_vec(nodes, nodes, raw_vector),
        }
    }
    /// Clustering. Note that the resulting clustering
    /// might be overlapping to each other.
    pub fn clustering(&self, e: i32, r: i32) -> Vec<Vec<usize>> {
        let mut cluster = self.matrix.clone();
        loop {
            let prev = cluster.clone();
            // Exp
            for _ in 0..(e - 1) {
                cluster *= prev.clone();
            }
            // Inflation
            cluster.column_iter_mut().for_each(|mut column| {
                let denom = column.iter().map(|x| x.powi(r)).sum::<f64>();
                column.iter_mut().for_each(|x| *x = x.powi(r) / denom);
                assert!((1. - column.sum()).abs() < 0.001);
            });
            // for row in cluster.row_iter() {
            //     let line: Vec<_> = row.iter().map(|x| format!("{:.2}", x)).collect();
            //     println!("{}", line.join("  "));
            // }
            let diff = (cluster.clone() - prev).norm();
            debug!("Diff:{}", diff);
            if diff < 0.001 {
                break;
            }
        }
        cluster
            .row_iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|&(_, &x)| x > 0.001)
                    .map(|x| x.0)
                    .collect::<Vec<_>>()
            })
            .filter(|row| !row.is_empty())
            .collect()
    }
}
