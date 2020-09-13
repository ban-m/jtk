use definitions::DataSet;
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;
use serde::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiplicityEstimationConfig {
    max_multiplicity: usize,
    seed: u64,
    path: Option<String>,
}

impl MultiplicityEstimationConfig {
    pub fn new(max_multiplicity: usize, seed: u64, path: Option<&str>) -> Self {
        Self {
            max_multiplicity,
            seed,
            path: path.map(|x| x.to_string()),
        }
    }
}

pub trait MultiplicityEstimation {
    fn estimate_multiplicity(self, config: &MultiplicityEstimationConfig) -> Self;
}

impl MultiplicityEstimation for DataSet {
    fn estimate_multiplicity(mut self, config: &MultiplicityEstimationConfig) -> Self {
        for read in self.encoded_reads.iter_mut() {
            for node in read.nodes.iter_mut() {
                node.cluster = 0;
            }
        }
        self.assignments = self
            .encoded_reads
            .iter()
            .map(|r| definitions::Assignment::new(r.id, 0))
            .collect();
        use super::Assemble;
        let assemble_config = super::AssembleConfig::default();
        let graphs = self.assemble_as_graph(&assemble_config);
        debug!("GRAPH\tID\tCoverage\tMean\tLen");
        let estimated_cluster_num: HashMap<u64, usize> = graphs
            .iter()
            .map(|graph| estimate_graph_multiplicity(&self, graph, config))
            .fold(HashMap::new(), |mut x, result| {
                for (unit, cluster) in result {
                    x.insert(unit, cluster);
                }
                x
            });
        for unit in self.selected_chunks.iter_mut() {
            if let Some(&cl_num) = estimated_cluster_num.get(&unit.id) {
                unit.cluster_num = cl_num;
            }
        }
        self
    }
}

fn estimate_graph_multiplicity(
    ds: &DataSet,
    graph: &super::assemble::Graph,
    c: &MultiplicityEstimationConfig,
) -> Vec<(u64, usize)> {
    let covs: Vec<_> = graph
        .nodes
        .iter()
        .map(|node| {
            let len = node.segments.len();
            let unit: HashSet<_> = node.segments.iter().map(|t| t.unit).collect();
            let coverage = ds
                .encoded_reads
                .iter()
                .map(|r| r.nodes.iter().filter(|n| unit.contains(&n.unit)).count())
                .sum::<usize>();
            let mean = (coverage / len) as u64;
            debug!("GRAPH\t{}\t{}\t{}\t{}", node.id, coverage, mean, len);
            mean
        })
        .collect();
    let (model, aic): (Model, f64) = (1..c.max_multiplicity)
        .map(|k| {
            let seed = k as u64 + c.seed;
            let (model, lk) = clustering(&covs, k, seed);
            let aic = -2. * lk + 2. * k as f64;
            (model, aic)
        })
        .min_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
        .unwrap();
    debug!("AIC\t{}", aic);
    let assignments: Vec<_> = covs.iter().map(|&d| model.assign(d)).collect();
    let min_coverage = assignments
        .iter()
        .map(|&x| model.lambdas[x])
        .min_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap();
    let repeat_num: Vec<_> = model
        .lambdas
        .iter()
        .map(|x| ((x / min_coverage) + 0.5).floor() as usize)
        .collect();
    debug!("LAMBDAS:{:?}", model.lambdas);
    debug!("PREDCT:{:?}", repeat_num);
    let mut result = vec![];
    for (&cl, contig) in assignments.iter().zip(graph.nodes.iter()) {
        let repeat_num = repeat_num[cl];
        debug!("REPEATNUM\t{}\t{}\t{}", contig.id, repeat_num, cl);
        for node in contig.segments.iter() {
            result.push((node.unit, repeat_num));
        }
    }
    result
}

struct Model {
    cluster: usize,
    fractions: Vec<f64>,
    lambdas: Vec<f64>,
}
const SMALL: f64 = 0.00000000000000001;
impl Model {
    fn new(data: &[u64], weight: &[Vec<f64>], k: usize) -> Self {
        let sum: Vec<_> = (0..k)
            .map(|cl| weight.iter().map(|ws| ws[cl]).sum::<f64>() + SMALL)
            .collect();
        let fractions: Vec<_> = sum.iter().map(|w| w / data.len() as f64).collect();
        let lambdas: Vec<_> = sum
            .iter()
            .enumerate()
            .map(|(cl, sum)| {
                weight
                    .iter()
                    .zip(data)
                    .map(|(ws, &x)| x as f64 * ws[cl])
                    .sum::<f64>()
                    / sum
            })
            .collect();
        let cluster = k;
        Self {
            cluster,
            fractions,
            lambdas,
        }
    }
    fn lk(&self, data: &[u64]) -> f64 {
        data.iter().map(|&d| self.lk_data(d)).sum::<f64>()
    }
    fn lk_data(&self, data: u64) -> f64 {
        let lks: Vec<_> = (0..self.cluster)
            .map(|cl| {
                self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
                    - self.lambdas[cl]
                    - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
            })
            .collect();
        logsumexp(&lks)
    }
    fn new_weight(&self, data: u64) -> Vec<f64> {
        let lks: Vec<_> = (0..self.cluster)
            .map(|cl| {
                self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
                    - self.lambdas[cl]
                    - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
            })
            .collect();
        let lk = logsumexp(&lks);
        assert!((1. - lks.iter().map(|x| (x - lk).exp()).sum::<f64>()).abs() < 0.0001);
        lks.iter().map(|x| (x - lk).exp()).collect()
    }
    fn assign(&self, data: u64) -> usize {
        let (cl, _) = (0..self.cluster)
            .map(|cl| {
                self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
                    - self.lambdas[cl]
                    - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
            })
            .enumerate()
            .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
            .unwrap();
        cl
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return 0.;
    }
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

fn clustering(data: &[u64], k: usize, seed: u64) -> (Model, f64) {
    let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(seed);
    let (weight, lk) = (0..10)
        .map(|_| {
            let mut weight: Vec<_> = (0..data.len())
                .map(|_| {
                    let mut ws = vec![0.; k];
                    ws[rng.gen::<usize>() % k] = 1.;
                    ws
                })
                .collect();
            let mut lk = std::f64::NEG_INFINITY;
            loop {
                let model = Model::new(data, &weight, k);
                let new_lk = model.lk(data);
                let diff = new_lk - lk;
                if diff < 0.00001 {
                    break;
                }
                lk = new_lk;
                weight = data.iter().map(|&d| model.new_weight(d)).collect();
            }
            (weight, lk)
        })
        .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
        .unwrap();
    let model = Model::new(data, &weight, k);
    (model, lk)
}

//     let mut coverages: Vec<u64> = vec![0; self.selected_chunks.len()];
//     for read in self.encoded_reads.iter() {
//         for node in read.nodes.iter() {
//             coverages[node.unit as usize] += 1;
//         }
//     }
//     // Construct the node grpah.
//     // [[edges from 0], [edges from 1], ... ,[edges from N]].concat();
//     let mut edges: Vec<_> = vec![0.; self.selected_chunks.len().pow(2u32)];
//     let dim = self.selected_chunks.len();
//     for read in self.encoded_reads.iter() {
//         for edge in read.edges.iter() {
//             edges[edge.from as usize * dim + edge.to as usize] += 1.;
//             edges[edge.to as usize * dim + edge.from as usize] += 1.;
//         }
//     }
//     assert_eq!(edges.len() % dim, 0);
//     edges.chunks_mut(dim).for_each(|es| {
//         let sum = es.iter().sum::<f64>();
//         es.iter_mut().for_each(|e| *e /= sum);
//     });
//     let multiplicities = estimate_multiplicity_by_em(&coverages, &edges, dim, config);
//     assert_eq!(coverages.len(), multiplicities.len());
//     self.selected_chunks
//         .iter_mut()
//         .zip(multiplicities)
//         .for_each(|(c, m)| c.cluster_num = m);
//     self
// }
// }

// #[derive(Clone)]
// struct Matrix {
//     column: usize,
//     row: usize,
//     // [X_00, X_01, X_02, ..., X_0C, X_10, ...., X_RC]
//     data: Vec<f64>,
// }

// impl std::fmt::Debug for Matrix {
//     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//         let lines: Vec<_> = self
//             .data
//             .chunks(self.row)
//             .enumerate()
//             .map(|(i, row_vec)| {
//                 let line: Vec<_> = row_vec.iter().map(|x| format!("{:.3}", x)).collect();
//                 format!("{}\t{}", i, line.join(","))
//             })
//             .collect();
//         write!(f, "{}", lines.join("\n"))
//     }
// }

// impl Matrix {
//     pub fn new(data: &[f64], row: usize, column: usize) -> Self {
//         assert_eq!(data.len(), row * column);
//         Self {
//             column,
//             row,
//             data: data.to_vec(),
//         }
//     }
//     pub fn transpose(&self) -> Self {
//         let mut data = vec![0.; self.row * self.column];
//         for r in 0..self.row {
//             for c in 0..self.column {
//                 data[r * self.column + c] = self.data[c * self.row + r];
//             }
//         }
//         Self {
//             column: self.row,
//             row: self.column,
//             data,
//         }
//     }
//     #[allow(dead_code)]
//     pub fn elementwise_mul(&self, other: &Matrix) -> Self {
//         assert_eq!(self.dim(), other.dim());
//         let data: Vec<_> = self
//             .data
//             .iter()
//             .zip(other.data.iter())
//             .map(|(x, y)| x * y)
//             .collect();
//         Self {
//             data,
//             column: self.column,
//             row: self.row,
//         }
//     }
//     pub fn rowwise_normalize(&mut self) {
//         self.data.chunks_mut(self.row).for_each(|xs| {
//             let sum = xs.iter().sum::<f64>();
//             xs.iter_mut().for_each(|x| *x /= sum);
//         });
//     }
//     pub fn diff(&self, other: &Self) -> f64 {
//         assert_eq!(self.column, other.column);
//         assert_eq!(self.row, other.row);
//         self.data
//             .iter()
//             .zip(other.data.iter())
//             .map(|(x, y)| (x - y).powi(2))
//             .sum::<f64>()
//     }
//     pub fn get(&self, i: usize, j: usize) -> f64 {
//         self.data[i * self.row + j]
//     }
//     pub fn rows(&self) -> std::slice::Chunks<f64> {
//         self.data.chunks(self.row)
//     }
//     // Return (column dimension, row dimension)
//     #[allow(dead_code)]
//     fn dim(&self) -> (usize, usize) {
//         (self.column, self.row)
//     }
// }

// impl std::ops::Mul<&Matrix> for &Matrix {
//     type Output = Matrix;
//     fn mul(self, rhs: &Matrix) -> Self::Output {
//         assert_eq!(self.row, rhs.column);
//         let mut data = vec![0.; self.column * rhs.row];
//         for c in 0..self.column {
//             for r in 0..rhs.row {
//                 let inner_prod = (0..self.row)
//                     .map(|k| self.data[c * self.row + k] * rhs.data[k * rhs.row + r])
//                     .sum::<f64>();
//                 data[c * rhs.row + r] = inner_prod;
//             }
//         }
//         Self::Output {
//             data,
//             column: self.column,
//             row: rhs.row,
//         }
//     }
// }

// pub fn prob_of_poiss(obs: u64, lambda: f64) -> f64 {
//     ((obs as f64) * lambda.ln() - (0..obs).map(|t| ((t + 1) as f64).ln()).sum::<f64>() - lambda)
//         .exp()
// }

// pub fn estimate_multiplicity_by_em(
//     coverages: &[u64],
//     edges: &[f64],
//     dim: usize,
//     config: &MultiplicityEstimationConfig,
// ) -> Vec<usize> {
//     // Initialize
//     // initial distribution of units.
//     let total_coverage = coverages.iter().sum::<u64>();
//     let density_on_units: Vec<_> = coverages
//         .iter()
//         .map(|&x| x as f64 / total_coverage as f64)
//         .collect();
//     let graph = Matrix::new(edges, dim, dim);
//     let graph_transpose = graph.transpose();
//     let state_num = config.max_multiplicity;
//     // Estimate parameters for Poisson distribution.
//     let mut lambdas: Vec<f64> = {
//         let mut counts = vec![0.; state_num];
//         let mut sums = vec![0.; state_num];
//         for (idx, &x) in coverages.iter().enumerate() {
//             counts[idx % state_num] += 1.;
//             sums[idx % state_num] += x as f64;
//         }
//         sums.iter().zip(counts).map(|(x, y)| x / y).collect()
//     };
//     // Initial distribution of state.
//     let mut density_on_states: Vec<_> =
//         (0..state_num).map(|_| (state_num as f64).recip()).collect();
//     // Transition matrix of states.
//     let state_transition: Vec<_> = vec![(state_num as f64).recip(); state_num * state_num];
//     // trace!("Density of nodes:{:?}", density_on_units);
//     // trace!("Nodes transition:\n{:?}", graph);
//     let mut state_transition = Matrix::new(&state_transition, state_num, state_num);
//     let weights = loop {
//         trace!("Current lambdas:{:?}", lambdas);
//         trace!("Current fractions:{:?}", density_on_states);
//         trace!("Current transition:\n{:?}", state_transition);
//         let obs_prob: Vec<f64> = (0..dim)
//             .flat_map(|unit| {
//                 (0..state_num)
//                     .map(|state| prob_of_poiss(coverages[unit], lambdas[state]))
//                     .collect::<Vec<_>>()
//             })
//             .collect();
//         let obs_prob = Matrix::new(&obs_prob, state_num, dim);
//         trace!("Current obs prob\n{:?}", obs_prob);
//         // 1. E-step: Estimate gamma_{datum,state} and kappa_{state, state}
//         let (gammas, mut kappas) = {
//             let mut gammas = &Matrix::new(&density_on_units, 1, dim)
//                 * &Matrix::new(&density_on_states, state_num, 1);
//             assert_eq!((gammas.row, gammas.column), (state_num, dim));
//             gammas.rowwise_normalize();
//             let mut diff = 1f64;
//             // Multiple until convergent.
//             while diff.abs() > 0.1 {
//                 let mut next_gammas = &(&graph_transpose * &gammas) * &state_transition;
//                 // let mut next_gammas =
//                 //     (&(&graph_transpose * &gammas) * &state_transition).elementwise_mul(&obs_prob);
//                 next_gammas.rowwise_normalize();
//                 diff = next_gammas.diff(&gammas);
//                 gammas = next_gammas;
//             }
//             // gammas = gammas.elementwise_mul(&obs_prob);
//             // gammas.rowwise_normalize();
//             let mut kappas = vec![0.; state_num * state_num];
//             for i in 0..state_num {
//                 for j in 0..state_num {
//                     let weight = (0..dim)
//                         .map(|node1| {
//                             (0..dim)
//                                 .map(|node2| {
//                                     let p1 = prob_of_poiss(coverages[node1], lambdas[i]);
//                                     let p2 = prob_of_poiss(coverages[node2], lambdas[j]);
//                                     p1 * p2
//                                         * graph.get(node1, node2)
//                                         * gammas.get(node1, i)
//                                         * state_transition.get(i, j)
//                                 })
//                                 .sum::<f64>()
//                         })
//                         .sum::<f64>();
//                     kappas[i * state_num + j] = weight;
//                 }
//             }
//             let mut kappas = Matrix::new(&kappas, state_num, state_num);
//             kappas.rowwise_normalize();
//             (gammas, kappas)
//         };
//         trace!("Probs\n{:?}", gammas);
//         trace!("Transs\n{:?}", kappas);
//         // 2. M-step: Update units_transtion and density on states.
//         lambdas = (0..state_num)
//             .map(|state| {
//                 let (count, sum) =
//                     coverages
//                         .iter()
//                         .enumerate()
//                         .fold((0., 0.), |(count, sum), (node, &x)| {
//                             let weight = gammas.get(node, state);
//                             (count + weight, sum + weight * x as f64)
//                         });
//                 sum / count
//             })
//             .collect();
//         density_on_states = (0..state_num)
//             .map(|state| (0..dim).map(|node| gammas.get(node, state)).sum::<f64>())
//             .map(|x| x / dim as f64)
//             .collect();
//         kappas.rowwise_normalize();
//         let new_state_transition = kappas;
//         let diff = new_state_transition.diff(&state_transition);
//         state_transition = new_state_transition;
//         if diff < 0.00001 {
//             break gammas;
//         }
//     };
//     weights
//         .rows()
//         .map(|row| {
//             let (argmax, _) = row
//                 .iter()
//                 .enumerate()
//                 .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
//                 .unwrap();
//             argmax
//         })
//         .collect()
// }

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn matrix_test() {
//         let m1 = Matrix::new(&vec![1., 2., 3.], 3, 1);
//         let m2 = Matrix::new(&vec![1., 2., 3.], 1, 3);
//         let prod = &m1 * &m2;
//         assert!((prod.get(0, 0) - 14.).abs() < 0.000001);
//         let diad = &m2 * &m1;
//         assert_eq!(diad.dim(), (3, 3));
//         for i in 0..3 {
//             for j in 0..3 {
//                 assert!((diad.get(i, j) - ((i + 1) * (j + 1)) as f64).abs() < 0.00001);
//             }
//         }
//         let m1 = Matrix::new(&vec![1., 2., 3., 4., 5., 6.], 3, 2);
//         let m2 = m1.transpose();
//         let prod = &m1 * &m2;
//         assert_eq!(prod.dim(), (2, 2));
//         assert!((prod.get(0, 0) - 14.).abs() < 0.0001);
//     }
//     use super::*;
//     use rand::{Rng, SeedableRng};
//     use rand_distr::{Distribution, Poisson};
//     use rand::{Rng, SeedableRng};
//     use rand_xoshiro::Xoshiro256PlusPlus;
//     fn create_mock_data_1(
//         num_unit: usize,
//         coverage: f64,
//     ) -> (Vec<u64>, Vec<f64>, usize, MultiplicityEstimationConfig) {
//         let pois = Poisson::new(coverage).unwrap();
//         let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(239);
//         let coverages: Vec<_> = (0..num_unit).map(|_| pois.sample(&mut rng)).collect();
//         let edges: Vec<_> = (0..num_unit)
//             .flat_map(|i| {
//                 let mut es = vec![0.; num_unit];
//                 if i > 0 {
//                     es[i - 1] += 1.;
//                 }
//                 if i + 1 < num_unit {
//                     es[i + 1] += 1.;
//                 }
//                 let sum = es.iter().sum::<f64>();
//                 es.iter_mut().for_each(|x| *x /= sum);
//                 es
//             })
//             .collect();
//         let config = MultiplicityEstimationConfig::new(2);
//         (coverages, edges, num_unit, config)
//     }
//     #[test]
//     fn test_mock_data_1() {
//         let num_unit = 5;
//         let coverage = 20.;
//         let (coverages, edges, dim, config) = create_mock_data_1(num_unit, coverage);
//         eprintln!("start");
//         let result = estimate_multiplicity_by_em(&coverages, &edges, dim, &config);
//         println!("{:?}", result);
//         assert!(result.iter().all(|&x| x == result[0]));
//     }
//     fn create_mock_data_2() -> (
//         Vec<u64>,
//         Vec<f64>,
//         usize,
//         MultiplicityEstimationConfig,
//         Vec<usize>,
//     ) {
//         let coverage = 20.;
//         let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(239);
//         let answer = vec![2, 2, 2, 1, 1, 1, 2, 2, 2];
//         let num_unit = answer.len();
//         let coverages: Vec<_> = answer
//             .iter()
//             .map(|&s| Poisson::new(coverage * s as f64).unwrap().sample(&mut rng))
//             .collect();
//         let edges: Vec<_> = (0..num_unit)
//             .flat_map(|i| {
//                 let mut es = vec![0.; num_unit];
//                 if i > 0 {
//                     es[i - 1] += 1.;
//                 }
//                 if i + 1 < num_unit {
//                     es[i + 1] += 1.;
//                 }
//                 let sum = es.iter().sum::<f64>();
//                 es.iter_mut().for_each(|x| *x /= sum);
//                 es
//             })
//             .collect();
//         let config = MultiplicityEstimationConfig::new(2);
//         (coverages, edges, num_unit, config, answer)
//     }
//     fn create_mock_data_3() -> (
//         Vec<u64>,
//         Vec<f64>,
//         usize,
//         MultiplicityEstimationConfig,
//         Vec<usize>,
//     ) {
//         let coverage = 20.;
//         let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(239);
//         let answer = vec![2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2];
//         let num_unit = answer.len();
//         let coverages: Vec<_> = answer
//             .iter()
//             .map(|&s| Poisson::new(coverage * s as f64).unwrap().sample(&mut rng))
//             .collect();
//         let mut edges: Vec<_> = vec![];
//         //                0   1   2   3   4   5   6   7   8   9   10
//         edges.extend(vec![0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.]); // 0
//         edges.extend(vec![2., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.]); // 1
//         edges.extend(vec![0., 2., 0., 1., 0., 1., 0., 0., 0., 0., 0.]); // 2
//         edges.extend(vec![0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0.]); // 3
//         edges.extend(vec![0., 0., 0., 1., 0., 0., 0., 1., 0., 0., 0.]); // 4
//         edges.extend(vec![0., 0., 1., 0., 0., 0., 1., 0., 0., 0., 0.]); // 5
//         edges.extend(vec![0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0.]); // 6
//         edges.extend(vec![0., 0., 0., 0., 1., 0., 1., 0., 2., 0., 0.]); // 7
//         edges.extend(vec![0., 0., 0., 0., 0., 0., 0., 2., 0., 2., 0.]); // 8
//         edges.extend(vec![0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 2.]); // 9
//         edges.extend(vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0.]); // 10
//         edges.chunks_mut(num_unit).for_each(|es| {
//             let sum = es.iter().sum::<f64>();
//             es.iter_mut().for_each(|x| *x /= sum)
//         });
//         let config = MultiplicityEstimationConfig::new(2);
//         (coverages, edges, num_unit, config, answer)
//     }
//     fn create_mock_data_4() -> (
//         Vec<u64>,
//         Vec<f64>,
//         usize,
//         MultiplicityEstimationConfig,
//         Vec<usize>,
//     ) {
//         let coverage = 20.;
//         let mut rng: Xoshiro256PlusPlus = SeedableRng::seed_from_u64(239);
//         let answer = vec![1, 1, 2, 1, 1, 1, 1, 1];
//         let num_unit = answer.len();
//         let coverages: Vec<_> = answer
//             .iter()
//             .map(|&s| Poisson::new(coverage * s as f64).unwrap().sample(&mut rng))
//             .collect();
//         let mut edges: Vec<_> = vec![];
//         //                0   1   2   3   4   5   6   7
//         edges.extend(vec![0., 1., 0., 0., 0., 0., 0., 0.]); // 0
//         edges.extend(vec![1., 0., 1., 0., 0., 0., 0., 0.]); // 1
//         edges.extend(vec![0., 1., 0., 1., 0., 1., 1., 0.]); // 2
//         edges.extend(vec![0., 0., 1., 0., 1., 0., 0., 0.]); // 3
//         edges.extend(vec![0., 0., 0., 1., 0., 1., 0., 0.]); // 4
//         edges.extend(vec![0., 0., 1., 0., 1., 0., 0., 0.]); // 5
//         edges.extend(vec![0., 0., 1., 0., 0., 0., 0., 1.]); // 6
//         edges.extend(vec![0., 0., 0., 0., 0., 0., 1., 0.]); // 7
//         edges.chunks_mut(num_unit).for_each(|es| {
//             let sum = es.iter().sum::<f64>();
//             es.iter_mut().for_each(|x| *x /= sum)
//         });
//         let config = MultiplicityEstimationConfig::new(2);
//         (coverages, edges, num_unit, config, answer)
//     }
//     #[test2]
//     fn test_mock_data_2() {
//         let (coverages, edges, dim, config, answer) = create_mock_data_1();
//         let _result = estimate_multiplicity_by_em(&coverages, &edges, dim, &config);
//         let (coverages, edges, dim, config, answer) = create_mock_data_2();
//         let _result = estimate_multiplicity_by_em(&coverages, &edges, dim, &config);
//         let (coverages, edges, dim, config, answer) = create_mock_data_3();
//         let _result = estimate_multiplicity_by_em(&coverages, &edges, dim, &config);
//         let (coverages, edges, dim, config, answer) = create_mock_data_4();
//         let _result = estimate_multiplicity_by_em(&coverages, &edges, dim, &config);
//     }
// }
