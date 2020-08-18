use super::create_model::get_models;
use super::ChunkedUnit;
use super::ClusteringConfig;
use super::SMALL_WEIGHT;
use nalgebra::DMatrix;
use poa_hmm::POA;
use rand::Rng;
use rayon::prelude::*;
pub fn get_variants<F: Fn(u8, u8) -> i32 + std::marker::Sync, R: Rng>(
    data: &[ChunkedUnit],
    chain_len: usize,
    rng: &mut R,
    c: &ClusteringConfig<F>,
    //_max_value: f64,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>, f64) {
    let update_data = vec![false; data.len()];
    let usepos = vec![true; chain_len];
    let ws = get_cluster_fraction(data, &update_data, c.cluster_num);
    let (mut betas, mut lks) = {
        let ms = get_models(data, chain_len, rng, c, &usepos, None);
        let (vars, lk) = variant_calling(&data, &c.poa_config, &ws, &ms);
        (vars, vec![lk])
    };
    for _ in 0..c.repeat_num - 1 {
        let ms = get_models(data, chain_len, rng, c, &usepos, None);
        let (vars, lk) = variant_calling(&data, &c.poa_config, &ws, &ms);
        lks.push(lk);
        betas.iter_mut().zip(vars).for_each(|(bss, vss)| {
            bss.iter_mut().zip(vss).for_each(|(bs, vs)| {
                bs.iter_mut().zip(vs).for_each(|(b, v)| *b += v);
            })
        });
    }
    betas.iter_mut().for_each(|bss| {
        bss.iter_mut().for_each(|bs| {
            let sum = bs.iter().map(|x| x * x).sum::<f64>().sqrt();
            bs.iter_mut().for_each(|b| *b = *b / sum)
        })
    });
    let lk = logsumexp(&lks);
    // for bss in betas.iter() {
    //     for bs in bss.iter() {
    //         let line: Vec<_> = bs
    //             .iter()
    //             .enumerate()
    //             .map(|(i, b)| format!("{}:{:.3}", i, b))
    //             .collect();
    //         debug!("{:?}", line.join(","));
    //     }
    // }
    let (betas, position_in_use) = select_variants(betas, chain_len, c.variant_num);
    // let betas = normalize_weights(betas, max_value);
    (betas, position_in_use, lk)
}

// Calculate LK matrices for each read. Total LK would be also returned.
// Note that each matrix is flattened in row-order.
fn calc_matrices_poa(
    data: &[ChunkedUnit],
    c: &poa_hmm::Config,
    ws: &[f64],
    models: &[Vec<POA>],
) -> (Vec<Vec<f64>>, f64) {
    let num_cluster = models.len();
    let chain_len = models[0].len();
    // LK matrix of each read, i.e., lk_matrix[i] would
    // be lk matrix for the i-th read.
    // Note that the matrix is flattened.
    // To access the likelihood of the j-th position of the k-th cluster,
    // lk_matrix[i][chain_len * k + j] would work.
    let lk_matrices: Vec<Vec<f64>> = data
        .par_iter()
        .map(|read| lks_poa(models, read, c, chain_len))
        .collect();
    let lk = lk_matrices
        .iter()
        .map(|matrix| {
            assert_eq!(matrix.len() / chain_len, num_cluster);
            let lks: Vec<_> = matrix
                .chunks_exact(chain_len)
                .zip(ws.iter())
                .map(|(chunks, w)| w.ln() + chunks.iter().sum::<f64>())
                .collect();
            logsumexp(&lks)
        })
        .sum::<f64>();
    (lk_matrices, lk)
}

// Return likelihoods of the read.
// the cl * i + j -th element is the likelihood for the i-th cluster at the j-th position.
fn lks_poa(
    models: &[Vec<POA>],
    read: &ChunkedUnit,
    config: &poa_hmm::Config,
    cl: usize,
) -> Vec<f64> {
    let mut res = vec![0.; cl * models.len()];
    for (i, ms) in models.iter().enumerate() {
        for c in read.chunks.iter() {
            res[i * cl + c.pos] = ms[c.pos].forward(&c.seq, config);
        }
    }
    res
}

fn centrize(mut matrices: Vec<Vec<f64>>, row: usize, column: usize) -> Vec<Vec<f64>> {
    let cv = centrize_vector_of(&matrices, row, column);
    // Centrize the matrices. In other words,
    // we add the offset cv for each position.
    // Note that we only care at the position where the likelihood is non-zero value.
    matrices.iter_mut().for_each(|matrix| {
        for pos in 0..column {
            if (0..row).any(|cluster| matrix[column * cluster + pos].abs() > 0.001) {
                for cluster in 0..row {
                    matrix[column * cluster + pos] += cv[cluster][pos];
                }
            }
        }
    });
    matrices
}

fn maximize_margin_of(matrices: &[Vec<f64>], row: usize, column: usize) -> Option<Vec<f64>> {
    let matrix = matrices
        .iter()
        .map(|matrix| DMatrix::from_row_slice(row, column, &matrix))
        .fold(DMatrix::zeros(column, column), |x, l| {
            let trans = l.transpose();
            let ones = DMatrix::repeat(row, row, 1.);
            let unit = DMatrix::identity(row, row);
            let reg = unit - ones / row as f64;
            x + trans * reg * l
        });
    let eigens = matrix.symmetric_eigen();
    use std::cmp::Ordering::Equal;
    let max = eigens
        .eigenvalues
        .iter()
        .map(|&e| e as f64)
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Equal))?;
    if max.1.is_nan() || max.1.is_infinite() {
        None
    } else {
        Some(eigens.eigenvectors.column(max.0).iter().copied().collect())
    }
}

// Return the centrize vector for matrices.
// In other words, it calculates arrays of vector p_1,...,p_{row}, where
// each c_i is column-dimensional vector.
fn centrize_vector_of(matrices: &[Vec<f64>], row: usize, column: usize) -> Vec<Vec<f64>> {
    let mut centrize_vector = vec![vec![0.; column]; row];
    // How many units is considered at position i.
    let mut counts = vec![0; column];
    // Compute the avarage vector at each **column(position)**
    for matrix in matrices.iter() {
        for pos in 0..column {
            if (0..row).any(|r| matrix[column * r + pos].abs() > 0.001) {
                counts[pos] += 1;
                for cluster in 0..row {
                    centrize_vector[cluster][pos] += matrix[column * cluster + pos];
                }
            }
        }
    }
    centrize_vector.iter_mut().for_each(|row_vec| {
        row_vec
            .iter_mut()
            .zip(counts.iter())
            .for_each(|(sum, &count)| *sum /= count as f64)
    });
    // Compute projection to axis, i.e., <c,1>/K 1 - c.
    // Here, c is the column vector at each position.
    // This is a centrize vector.
    // First, compute the <c,1> for each column.
    // We could have this by just add by columnwise.
    let mut inner_product = centrize_vector
        .iter()
        .fold(vec![0.; column], |mut acc, row_vec| {
            acc.iter_mut()
                .zip(row_vec.iter())
                .for_each(|(x, y)| *x += y);
            acc
        });
    // Then, normalize it by dividing by length of row(number of cluster).
    inner_product.iter_mut().for_each(|x| *x /= row as f64);
    // Lastly, compute <c,1>/K 1 - c for each position.
    centrize_vector.iter_mut().for_each(|cv| {
        cv.iter_mut()
            .zip(inner_product.iter())
            .for_each(|(x, prod)| *x = prod - *x);
    });
    centrize_vector
}

// Call variant of each pairs of cluster.
// Note that returned matrix is lower-triangled. In other words,
// xs[i][j] would be varid only if j < i.
fn variant_calling(
    data: &[ChunkedUnit],
    c: &poa_hmm::Config,
    ws: &[f64],
    models: &[Vec<POA>],
) -> (Vec<Vec<Vec<f64>>>, f64) {
    let (matrices, lk) = calc_matrices_poa(data, c, ws, models);
    let cluster_num = models.len();
    let chain_len = models[0].len();
    let betas: Vec<Vec<Vec<f64>>> = (0..cluster_num)
        .map(|i| {
            (0..i)
                .map(|j| call_variants(i, j, &matrices, chain_len))
                .map(normalize)
                .collect()
        })
        .collect();
    (betas, lk)
}

fn normalize(mut xs: Vec<f64>) -> Vec<f64> {
    let sign = xs
        .iter()
        .max_by(|a, b| {
            a.abs()
                .partial_cmp(&b.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|x| x.signum())
        .unwrap_or(0.);
    xs.iter_mut().for_each(|x| *x *= sign);
    xs
}

// Call varinants between cluster i and cluster j.
fn call_variants(i: usize, j: usize, matrices: &[Vec<f64>], column: usize) -> Vec<f64> {
    // Extract focal rows for each read.
    let matrices: Vec<Vec<_>> = matrices
        .iter()
        .map(|matrix| {
            let use_pos: Vec<_> = matrix[i * column..(i + 1) * column]
                .iter()
                .zip(matrix[j * column..(j + 1) * column].iter())
                .map(|(&l1, &l2)| l1 > poa_hmm::DEFAULT_LK && l2 > poa_hmm::DEFAULT_LK)
                .collect();
            let class_i = matrix[i * column..(i + 1) * column]
                .iter()
                .zip(use_pos.iter())
                .map(|(&lk, &b)| if b { lk } else { 0. });
            let class_j = matrix[j * column..(j + 1) * column]
                .iter()
                .zip(use_pos.iter())
                .map(|(&lk, &b)| if b { lk } else { 0. });
            class_i.chain(class_j).collect()
        })
        .collect();
    let matrices = centrize(matrices, 2, column);
    match maximize_margin_of(&matrices, 2, column) {
        Some(res) => res,
        None => vec![0.; column],
    }
}
fn logsumexp(xs: &[f64]) -> f64 {
    let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}

fn select_variants(
    mut variants: Vec<Vec<Vec<f64>>>,
    chain_len: usize,
    variant_number: usize,
) -> (Vec<Vec<Vec<f64>>>, Vec<bool>) {
    let mut position = vec![false; chain_len];
    let thr = {
        let mut var: Vec<_> = variants
            .iter()
            .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
            .copied()
            .collect();
        var.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let pos = var.len() - variant_number;
        var[pos].max(0.05)
    };
    for bss in variants.iter_mut() {
        for bs in bss.iter_mut() {
            for (idx, b) in bs.iter_mut().enumerate() {
                if *b < thr {
                    *b = 0.;
                } else {
                    position[idx] = true;
                }
            }
        }
    }
    (variants, position)
}

fn get_cluster_fraction(
    data: &[ChunkedUnit],
    update_data: &[bool],
    cluster_num: usize,
) -> Vec<f64> {
    let total = update_data.iter().filter(|&b| !b).count() as f64;
    (0..cluster_num)
        .map(|cl| {
            data.iter()
                .zip(update_data)
                .filter(|&(d, &b)| !b && d.cluster == cl)
                .count()
        })
        .map(|count| (count as f64).max(SMALL_WEIGHT) / total)
        .collect()
}

// fn normalize_weights(mut betas: Vec<Vec<Vec<f64>>>, target: f64) -> Vec<Vec<Vec<f64>>> {
//     let max = betas
//         .iter()
//         .flat_map(|bss| bss.iter().flat_map(|bs| bs.iter()))
//         .fold(0., |x, &y| if x < y { y } else { x });
//     betas.iter_mut().for_each(|bss| {
//         bss.iter_mut().for_each(|betas| {
//             betas.iter_mut().for_each(|b| *b *= target / max);
//         })
//     });
//     betas
// }
