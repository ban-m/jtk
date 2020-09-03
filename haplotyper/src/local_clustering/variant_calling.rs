use super::create_model::get_models;
use super::ChunkedUnit;
use super::ClusteringConfig;
// use nalgebra::DMatrix;
use rand::Rng;
use rayon::prelude::*;
type VariantInformation = (Vec<Vec<Vec<f64>>>, Vec<bool>, Vec<Vec<f64>>, f64);
pub fn get_variants<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
    variant_num: usize,
) -> VariantInformation {
    let (betas, lks) = variant_calling(data, chain_len, rng, c);
    let (betas, position_in_use, margin) = select_variants(betas, chain_len, variant_num);
    let pos: Vec<_> = position_in_use
        .iter()
        .enumerate()
        .filter(|&(_, &b)| b)
        .map(|x| x.0)
        .collect();
    trace!("POS\t{:?}", pos);
    (betas, position_in_use, lks, margin)
}

// Calculate LK matrices for each read. Total LK would be also returned.
// Note that each matrix is flattened in row-order.
fn calc_matrices_poa<F, R>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> Vec<Vec<f64>>
//-> (Vec<Vec<f64>>, Vec<Vec<f64>>)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
    R: Rng,
{
    // let num_cluster = c.cluster_num;
    // LK matrix of each read, i.e., lk_matrix[i] would
    // be lk matrix for the i-th read.
    // Note that the matrix is flattened.
    // To access the likelihood of the j-th position of the k-th cluster,
    // lk_matrix[i][chain_len * k + j] would work.
    let poss = vec![true; chain_len];
    let models = get_models(data, chain_len, rng, c, &poss, None);
    let lk_matrices: Vec<Vec<f64>> = data
        .par_iter()
        .map(|read| lks_poa(&models, chain_len, c, read))
        .collect();
    lk_matrices
}

// Return likelihoods of the read.
// the cl * i + j -th element is the likelihood for the i-th cluster at the j-th position.
// fn lks_poa<F, R>(
fn lks_poa<F>(
    models: &[Vec<poa_hmm::POA>],
    chain_len: usize,
    c: &ClusteringConfig<F>,
    read: &ChunkedUnit,
) -> Vec<f64>
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
{
    let config = &c.poa_config;
    let mut res = vec![0.; c.cluster_num * chain_len];
    for (i, ms) in models.iter().enumerate() {
        for c in read.chunks.iter() {
            res[i * chain_len + c.pos] += ms[c.pos].forward(&c.seq, config);
        }
    }
    res
}

// Call variant of each pairs of cluster.
// Note that returned matrix is lower-triangled. In other words,
// xs[i][j] would be varid only if j < i.
fn variant_calling<F, R>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
) -> (Vec<Vec<Vec<f64>>>, Vec<Vec<f64>>)
where
    F: Fn(u8, u8) -> i32 + std::marker::Sync,
    R: Rng,
{
    let matrices = calc_matrices_poa(data, chain_len, rng, c);
    let betas: Vec<Vec<Vec<f64>>> = (0..c.cluster_num)
        .map(|i| {
            (0..i)
                .map(|j| call_variants(i, j, data, &matrices, chain_len))
                .collect()
        })
        .collect();
    (betas, matrices)
}

// fn normalize(mut xs: Vec<f64>) -> Vec<f64> {
//     let sign = xs
//         .iter()
//         .max_by(|a, b| {
//             a.abs()
//                 .partial_cmp(&b.abs())
//                 .unwrap_or(std::cmp::Ordering::Equal)
//         })
//         .map(|x| x.signum())
//         .unwrap_or(0.);
//     xs.iter_mut().for_each(|x| *x *= sign);
//     xs
// }

// Call varinants between cluster i and cluster j.
fn call_variants(
    i: usize,
    j: usize,
    data: &[ChunkedUnit],
    matrices: &[Vec<f64>],
    column: usize,
) -> Vec<f64> {
    // for (idx, d) in data.iter().enumerate() {
    //     trace!("ASN\t{}\t{}", idx, d.cluster);
    // }
    // Extract focal rows for each read.
    // let matrices: Vec<Vec<_>> = matrices
    //     .iter()
    //     .zip(data.iter())
    //     .filter(|(_, d)| d.cluster == i || d.cluster == j)
    //     .map(|(matrix, _)| {
    //         let class_i = matrix[i * column..(i + 1) * column].iter();
    //         let class_j = matrix[j * column..(j + 1) * column].iter();
    //         class_i.chain(class_j).copied().collect()
    //     })
    //     .collect();
    // let matrices = centrize(matrices, 2, column);
    // let margins = matrices
    //     .iter()
    //     .zip(data.iter())
    //     .filter(|(_, d)| d.cluster == i || d.cluster == j)
    //     .fold(vec![vec![]; column], |mut margins, (matrix, d)| {
    //         let cluster = if d.cluster == i { 1f64 } else { -1f64 };
    //         for pos in 0..column {
    //             if matrix[pos] < -0.001 || matrix[column + pos] < -0.001 {
    //                 let diff = matrix[i * column + pos] - matrix[j * column + pos];
    //                 margins[pos].push((cluster, diff));
    //             }
    //         }
    //         margins
    //     });
    let mut differences: Vec<_> = vec![vec![]; column];
    let mut indices = vec![vec![]; column];
    for (idx, (matrix, data)) in matrices.iter().zip(data.iter()).enumerate() {
        for chunk in data.chunks.iter() {
            let diff = matrix[i * column + chunk.pos] - matrix[j * column + chunk.pos];
            differences[chunk.pos].push(diff);
            indices[chunk.pos].push(idx);
        }
    }
    let result: Vec<_> = differences
        .iter()
        .zip(indices)
        .enumerate()
        .map(|(pos, (xs, idx))| {
            let (assignment, lk) = clustering(xs, pos as u64, 2);
            let (_, lk_0) = clustering(xs, pos as u64, 1);
            let gain = lk - lk_0 - 3.;
            let cluster = {
                let mut slot = vec![0; 2];
                for &asn in assignment.iter() {
                    slot[asn] += 1;
                }
                slot
            };
            // for ((asn, lk), i) in assignment.iter().zip(xs).zip(idx) {
            //     trace!("LK\t{}\t{}\t{}\t{}", pos, i, lk, asn);
            // }
            trace!(
                "VAL\t{}\t{:?}\t{:.0}\t{:.0}\t{:.0}",
                pos,
                cluster,
                lk,
                lk_0,
                gain,
            );
            gain.max(0.)
            // let result: Vec<_> = margins
            //     .into_iter()
            //     .map(|mut margin| {
            //         margin.sort_by(|x, y| x.1.abs().partial_cmp(&y.1.abs()).unwrap());
            //         //let cut = 3;
            //         // if margin.len() < 3 {
            //         let len = margin.len() as f64;
            //         let mean = margin.iter().map(|x| x.1).sum::<f64>() / len;
            //         margin
            //             .iter()
            //             .map(|(cluster, margin)| cluster * (margin - mean))
            //             .sum::<f64>()
            // } else {
            //     let len = margin.len() as f64 - cut as f64;
            //     let mean = margin.iter().rev().map(|x| x.1).skip(cut).sum::<f64>() / len;
            //     margin
            //         .iter()
            //         .rev()
            //         .skip(3)
            //         .map(|(cluster, margin)| cluster * (margin - mean))
            //         .sum::<f64>()
            // };
            // let mean = margin.iter().map(|x| x.1).sum::<f64>() / len;
            // let line: Vec<_> = margin
            //     .iter()
            //     .map(|(c, m)| format!("{:.3}", c * (m - mean)))
            //     .collect();
            // trace!("{}", line.join("\t"));
            // trace!("{}\t{}", idx, sum);
        })
        .collect();
    result
    // match maximize_margin_of(&matrices, 2, column) {
    //     Some(res) => normalize(res),
    //     None => vec![0.; column],
    // }
}

struct Gaussian {
    mean: f64,
    variance: f64,
    weight: f64,
}

impl std::fmt::Display for Gaussian {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}", self.mean, self.variance, self.weight)
    }
}

impl Gaussian {
    fn lk(&self, x: f64) -> f64 {
        if self.weight < 0.0001 || self.variance < 0.00001 {
            std::f64::NEG_INFINITY
        } else {
            -(2. * std::f64::consts::PI * self.variance).ln() / 2.
                - (x - self.mean).powi(2i32) / 2. / self.variance
                + self.weight.ln()
        }
    }
    fn update(xs: &[f64], weights: &[f64]) -> Self {
        let len = xs.len() as f64;
        let sum = weights.iter().sum::<f64>();
        let mean = xs
            .iter()
            .zip(weights.iter())
            .map(|(x, w)| x * w)
            .sum::<f64>()
            / sum;
        // TODO: Make this value as a parameter.
        let variance = 0.3;
        let weight = sum / len;
        Self {
            mean,
            variance,
            weight,
        }
    }
}

fn clustering(xs: &[f64], seed: u64, k: usize) -> (Vec<usize>, f64) {
    let variance = {
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        xs.iter().map(|x| (x - mean).powi(2i32)).sum::<f64>() / xs.len() as f64
    };
    if variance < 0.0001 {
        return (vec![0; xs.len()], 0.);
    }
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed);
    let mut weights: Vec<Vec<_>> = xs
        .iter()
        .map(|_| {
            let mut weights: Vec<_> = vec![0.; k];
            weights[rng.gen::<usize>() % k] = 1.;
            weights
        })
        .collect();
    let mut diff = 10000000000000000.;
    let mut lk = -10000000000000000.;
    while diff > 0.0000001 {
        let gaussian: Vec<_> = (0..k)
            .map(|i| {
                let weights: Vec<_> = weights.iter().map(|x| x[i]).collect();
                Gaussian::update(&xs, &weights)
            })
            .collect();
        weights = xs
            .iter()
            .map(|&x| {
                let lks: Vec<_> = gaussian.iter().map(|g| g.lk(x)).collect();
                let sum = logsumexp(&lks);
                lks.iter().map(|x| (x - sum).exp()).collect::<Vec<_>>()
            })
            .inspect(|x| assert!((1. - x.iter().sum::<f64>()).abs() < 0.01))
            .collect();
        let likelihood: f64 = xs
            .iter()
            .map(|&x| {
                let lks: Vec<_> = gaussian.iter().map(|g| g.lk(x)).collect();
                logsumexp(&lks)
            })
            .sum::<f64>();
        diff = likelihood - lk;
        lk = likelihood;
    }
    let gaussian: Vec<_> = (0..k)
        .map(|i| {
            let weights: Vec<_> = weights.iter().map(|x| x[i]).collect();
            Gaussian::update(&xs, &weights)
        })
        .collect();
    let result: Vec<_> = weights
        .iter()
        .map(|xs| {
            xs.iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(&(b.1)).unwrap())
                .map(|x| x.0)
                .unwrap()
        })
        .collect();
    let likelihood: f64 = xs
        .iter()
        .map(|&x| {
            let lks: Vec<_> = gaussian.iter().map(|g| g.lk(x)).collect();
            logsumexp(&lks)
        })
        .sum::<f64>();
    (result, likelihood)
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

fn select_variants(
    mut variants: Vec<Vec<Vec<f64>>>,
    chain_len: usize,
    variant_number: usize,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>, f64) {
    let mut position = vec![false; chain_len];
    // First, remove all the negative elements.
    let (thr, max) = {
        let mut var: Vec<_> = variants
            .iter()
            .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
            .copied()
            .collect();
        var.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let cluster_num = variants.len();
        let variant_number = variant_number * cluster_num * (cluster_num - 1) / 2;
        let pos = var.len() - variant_number.min(var.len());
        let max = *var.last().unwrap();
        (var[pos].max(0.2 * max + 0.05), max)
        // let margin = 0.2;
        // let sum = variants
        //     .iter()
        //     .map(|bss| {
        //         bss.iter()
        //             .map(|bs| bs.iter().map(|&x| x * x).sum::<f64>().sqrt())
        //             .sum::<f64>()
        //     })
        //     .sum::<f64>();
        // let num = variants.iter().map(|bss| bss.len()).sum::<usize>();
        // let mean = sum / num as f64;
        // let margin = mean * margin;
        // let max = variants
        //     .iter()
        //     .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
        //     .map(|x| x.abs())
        //     .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        //     .unwrap_or(0.);
        // trace!("MAX:{:.3},MARGIN:{:.3}", max, margin);
        // ((max - margin).max(0.05), max)
    };
    for bss in variants.iter_mut() {
        for bs in bss.iter_mut() {
            for (idx, b) in bs.iter_mut().enumerate() {
                if b.abs() < thr {
                    *b = 0.;
                } else {
                    *b = 1.;
                    position[idx] = true;
                }
            }
        }
    }
    (variants, position, max)
}
