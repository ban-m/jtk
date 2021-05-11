//! A small K-means clustering algorithm.
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    band_width: usize,
    probe_num: usize,
    retry_num: usize,
    cluster_num: u8,
    subbatch_num: usize,
}

impl ClusteringConfig {
    pub fn new(
        band_width: usize,
        probe_num: usize,
        retry_num: usize,
        cluster_num: u8,
        subbatch_num: usize,
    ) -> Self {
        Self {
            band_width,
            probe_num,
            retry_num,
            cluster_num,
            subbatch_num,
        }
    }
}

/// Return the assignments and the consensus sequence.
pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<(Vec<u8>, Vec<u8>)> {
    let ClusteringConfig {
        band_width,
        cluster_num,
        ..
    } = *config;
    let cons_template = kiley::consensus(reads, rng.gen(), 10, band_width);
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::clr();
    let template = kiley::padseq::PadSeq::new(cons_template.as_slice());
    let reads: Vec<_> = reads
        .iter()
        .map(|r| kiley::padseq::PadSeq::new(r.borrow()))
        .collect();
    let (hmm, _) = hmm.fit_banded_inner(&template, &reads, 100);
    let profiles: Vec<Vec<_>> = reads
        .iter()
        .map(|read| {
            let prof = banded::ProfileBanded::new(&hmm, &template, &read, 100).unwrap();
            let lk = prof.lk();
            prof.to_modification_table()
                .iter()
                .take(9 * template.len())
                .map(|x| x - lk)
                .collect()
        })
        .collect();
    let var_num = 2 * cluster_num as usize;
    let probes: Vec<_> = (0..9 * template.len())
        .map(|pos| {
            let mut xs: Vec<_> = profiles.iter().map(|xs| xs[pos]).collect();
            xs.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let sum: f64 = xs.iter().filter(|&&x| 0f64 < x).sum();
            let median = xs[xs.len() / 2];
            (pos, sum, median)
        })
        .collect();
    let mut probes: Vec<_> = probes
        .chunks_exact(9)
        .filter_map(|xs| xs.iter().max_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap()))
        .copied()
        .collect();
    probes.sort_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap());
    probes.reverse();
    probes.truncate(var_num);
    let norm: f64 = probes.iter().map(|(_, s, _)| s * s).sum();
    let norm = norm.sqrt();
    let profiles: Vec<Vec<_>> = profiles
        .iter()
        .map(|xs| {
            probes
                .iter()
                .map(|&(pos, scale, _)| (xs[pos] * scale / norm))
                .collect()
        })
        .collect();
    let (assignments, _score) = (0..10)
        .map(|_| {
            let asn = kmeans_f64_plusplus(&profiles, cluster_num, rng);
            let score = score(&profiles, &asn, cluster_num as usize);
            (asn, score)
        })
        .fold(
            (vec![], std::f64::NEG_INFINITY),
            |(argmax, max), (asn, score)| {
                if max < score {
                    (asn, score)
                } else {
                    (argmax, max)
                }
            },
        );
    Some((assignments, cons_template))
}

// Return the gain of likelihood for a given dataset.
// To calculate the gain,we compute the following metrics:
// 1. Sum up the vectors for each cluster.
// 2. For each sum, the element of the positve values would be selected,
// or we acceept the edit operation at that point.
// 3. Sum up the positive values of summing-upped vector for each cluster.
fn score(data: &[Vec<f64>], asn: &[u8], k: usize) -> f64 {
    let dim = data[0].len();
    let mut sums = vec![vec![0f64; dim]; k];
    for (xs, &asn) in data.iter().zip(asn.iter()) {
        xs.iter()
            .zip(sums[asn as usize].iter_mut())
            .for_each(|(x, y)| *y += x);
    }
    sums.iter()
        .map(|xs| -> f64 { xs.iter().filter(|&&x| 0f64 < x).sum() })
        .sum()
}

fn kmeans_f64_with_init(data: &[Vec<f64>], assignments: &mut [u8], k: u8) {
    let mut is_updated = true;
    while is_updated {
        let centers: Vec<Vec<f64>> = (0..k)
            .filter_map(|cl| {
                let (mut count, mut slots) = (0, vec![0f64; data[0].len()]);
                let filtered = data
                    .iter()
                    .zip(assignments.iter())
                    .filter_map(|(d, &a)| (a == cl).then(|| d));
                for datum in filtered {
                    assert_eq!(slots.len(), datum.len());
                    slots.iter_mut().zip(datum).for_each(|(acc, x)| *acc += x);
                    count += 1;
                }
                let center: Vec<_> = slots.iter().map(|&x| x as f64 / count as f64).collect();
                (count != 0).then(|| center)
            })
            .collect();
        is_updated = false;
        for (x, asn) in data.iter().zip(assignments.iter_mut()) {
            let (new_asn, _) = centers
                .iter()
                .enumerate()
                .map(|(i, center)| (i as u8, euclid_norm_f64(center, x)))
                .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
                .unwrap();
            if new_asn != *asn {
                is_updated = true;
                *asn = new_asn;
            }
        }
    }
}

#[allow(dead_code)]
fn kmeans_f64_plusplus<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
    let mut centers: Vec<&[f64]> = vec![];
    let indices: Vec<_> = (0..data.len()).collect();
    // Choosing centers.
    use rand::seq::SliceRandom;
    while centers.len() < k as usize {
        // calculate distance to the most nearest centers.
        let mut dists: Vec<_> = data
            .iter()
            .map(|xs| {
                centers
                    .iter()
                    .map(|c| euclid_norm_f64(xs, c))
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap_or(1f64)
                    .powi(2)
            })
            .collect();
        let total: f64 = dists.iter().sum();
        dists.iter_mut().for_each(|x| *x /= total);
        let idx = *indices.choose_weighted(rng, |&idx| dists[idx]).unwrap();
        centers.push(&data[idx]);
    }
    let mut assignments: Vec<_> = data
        .iter()
        .map(|xs| {
            centers
                .iter()
                .enumerate()
                .map(|(i, center)| (i as u8, euclid_norm_f64(center, xs)))
                .min_by(|x, y| x.1.partial_cmp(&(y.1)).unwrap())
                .unwrap()
                .0
        })
        .collect();
    kmeans_f64_with_init(data, &mut assignments, k);
    assignments
}

#[allow(dead_code)]
fn kmeans_f64<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
    kmeans_f64_with_init(data, &mut assignments, k);
    assignments
}

// Return the distance between xs and ys.
fn euclid_norm_f64(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(ys.len(), xs.len());
    xs.iter()
        .zip(ys.iter())
        .map(|(x, &y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}
