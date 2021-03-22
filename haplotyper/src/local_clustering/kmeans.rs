//! A small K-means clustering algorithm.
//! It run several times, pick the best one.
use kiley::alignment::bialignment::edit_dist;
use kiley::alignment::bialignment::polish_until_converge_banded;
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct ClusteringConfig {
    band_width: usize,
    probe_num: usize,
    retry_num: usize,
    cluster_num: u8,
    subbatch_num: usize,
    // Maybe dist?
    // Or, is there some canonical meaning?
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

pub fn clustering<R: Rng>(reads: &[Vec<u8>], rng: &mut R, config: &ClusteringConfig) -> Vec<u8> {
    let cons_template = kiley::consensus(&reads, 232, 4, 30).unwrap();
    let ClusteringConfig {
        band_width,
        probe_num,
        retry_num,
        cluster_num,
        subbatch_num,
    } = *config;
    for _ in 0..retry_num {
        let probes: Vec<_> = (0..probe_num)
            .filter_map(|_| {
                use rand::seq::SliceRandom;
                let picked: Vec<_> = reads
                    .choose_multiple(rng, subbatch_num)
                    .map(|x| x.as_slice())
                    .collect();
                polish_until_converge_banded(&cons_template, &picked, band_width)
            })
            .collect();
        let profiles: Vec<Vec<_>> = reads
            .iter()
            .map(|r| {
                let center = edit_dist(r, &cons_template) as i32;
                probes
                    .iter()
                    .map(|p| edit_dist(p, r) as i32 - center)
                    .collect()
            })
            .collect();
        let selected_position: Vec<_> = (0..probes.len())
            .map(|i| {
                let first = profiles[0][i];
                profiles.iter().any(|p| p[i] != first)
            })
            .collect();
        if selected_position.iter().any(|&x| x) {
            let profiles: Vec<Vec<_>> = profiles
                .iter()
                .map(|xs| {
                    xs.iter()
                        .zip(selected_position.iter())
                        .filter_map(|(&x, &b)| b.then(|| x))
                        .collect()
                })
                .collect();
            return kmeans(&profiles, cluster_num, rng);
        }
    }
    vec![0; reads.len()]
}

// TODO:Maybe we need more sophisticated clustering algorithm here.
// DBSCAN?
// Kmeans with random seeds?
fn kmeans<R: Rng>(data: &[Vec<i32>], k: u8, rng: &mut R) -> Vec<u8> {
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
    loop {
        let centers: Vec<Vec<f64>> = (0..k)
            .filter_map(|cl| {
                let (mut count, mut slots) = (0, vec![0; data[0].len()]);
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
        let new_assignments: Vec<_> = data
            .iter()
            .filter_map(|x| {
                centers
                    .iter()
                    .enumerate()
                    .map(|(i, center)| (i as u8, euclid_norm(center, x)))
                    .min_by(|x, y| (x.1).partial_cmp(&(y.1)).unwrap())
                    .map(|x| x.0)
            })
            .collect();
        if new_assignments == assignments {
            break assignments;
        } else {
            assignments = new_assignments;
        }
    }
}
fn euclid_norm(xs: &[f64], ys: &[i32]) -> f64 {
    assert_eq!(ys.len(), xs.len());
    xs.iter()
        .zip(ys.iter())
        .map(|(x, &y)| (x - y as f64).powi(2))
        .sum()
}
