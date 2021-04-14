//! A small K-means clustering algorithm.
//! It run several times, pick the best one.
// use kiley::bialignment::edit_dist;
use kiley::bialignment::polish_until_converge_banded;
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

pub fn clustering<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Vec<u8> {
    // Polish clustering in `retry_num` times.
    let reads: Vec<_> = reads.iter().map(|x| x.borrow()).collect();
    (0..config.retry_num)
        .find_map(|_| try_clustering(&reads, rng, config))
        .unwrap_or(vec![0; reads.len()])
}

pub fn clustering_rep<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<Vec<u8>> {
    // Polish clustering in `retry_num` times.
    let reads: Vec<_> = reads.iter().map(|x| x.borrow()).collect();
    let band_width = config.band_width;
    let cons_template = kiley::consensus_bounded(&reads, rng.gen(), 3, band_width, 2 * band_width)?;
    let mut assignments: Vec<_> = (0..reads.len())
        .map(|i| (i % config.cluster_num as usize) as u8)
        .collect();
    // TODO: parametrize here.
    for _ in 0..5 {
        let profiles: Vec<Vec<_>> = (0..config.retry_num)
            .find_map(|_| get_profiles(&cons_template, &assignments, &reads, rng, config))?;
        assignments = kmeans(&profiles, config.cluster_num, rng);
    }
    Some(assignments)
}

pub fn clustering_hmm<R: Rng, T: std::borrow::Borrow<[u8]>>(
    reads: &[T],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<Vec<u8>> {
    // Polish clustering in `retry_num` times.
    let reads: Vec<_> = reads.iter().map(|x| x.borrow()).collect();
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::new_three_state(0.9, 0.05, 0.10, 0.98);
    let picked = &reads[..reads.len().min(config.subbatch_num)];
    let cons_template = kiley::ternary_consensus(&picked, rng.gen(), 3, config.band_width);
    let cons_template =
        hmm.correct_until_convergence_banded(&cons_template, picked, config.band_width)?;
    debug!("Polished\t{}\t{}", cons_template.len(), reads.len());
    let mut assignments: Vec<_> = (0..reads.len())
        .map(|i| (i % config.cluster_num as usize) as u8)
        .collect();
    for _ in 0..5 {
        let profiles = get_profiles_hmm(&cons_template, &assignments, &reads, rng, config);
        assignments = kmeans_f64(&profiles, config.cluster_num, rng);
    }
    Some(assignments)
}

fn get_profiles_hmm<R: Rng>(
    template: &[u8],
    asn: &[u8],
    reads: &[&[u8]],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Vec<Vec<f64>> {
    use kiley::gphmm::*;
    let hmm = kiley::gphmm::GPHMM::<Cond>::new_three_state(0.9, 0.05, 0.10, 0.98);
    let probes: Vec<_> = (0..config.probe_num)
        .filter_map(|_| {
            let picked = pick_reads(asn, reads, config, rng);
            let cons = hmm.correct_until_convergence_banded(&template, &picked, config.band_width);
            debug!("Polished");
            cons
        })
        .collect();
    reads
        .iter()
        .map(|r| {
            let center = hmm
                .likelihood_banded(&template, &r, config.band_width)
                .unwrap();
            probes
                .iter()
                .map(|p| hmm.likelihood_banded(&p, &r, config.band_width).unwrap() - center)
                .collect::<Vec<_>>()
        })
        .collect()
}

fn get_profiles<R: Rng>(
    template: &[u8],
    asn: &[u8],
    reads: &[&[u8]],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<Vec<Vec<i32>>> {
    let probes: Vec<_> = (0..config.probe_num)
        .filter_map(|_| {
            let picked = pick_reads(asn, reads, config, rng);
            polish_until_converge_banded(&template, &picked, config.band_width)
        })
        .collect();
    let mismatch = |xs: &[u8], ys: &[u8]| {
        edlib_sys::global(xs, ys)
            .iter()
            .filter(|&&x| x == 3)
            .count()
    };
    let profiles: Vec<Vec<_>> = reads
        .iter()
        .map(|r| {
            // let center = edlib_sys::global_dist(r, &template) as i32;
            let center = mismatch(r, &template) as i32;
            probes
                .iter()
                //.map(|p| edlib_sys::global_dist(p, r) as i32 - center)
                .map(|p| mismatch(p, r) as i32 - center)
                .collect()
        })
        .collect();
    let selected_position: Vec<bool> = {
        // Transpose profiles.
        let mut t_profiles: Vec<_> = vec![vec![]; probes.len()];
        for prf in profiles.iter() {
            for (i, x) in prf.iter().enumerate() {
                t_profiles[i].push(x);
            }
        }
        t_profiles
            .iter()
            .map(|prf| prf.iter().any(|&x| x != prf[0]))
            .collect()
    };
    selected_position.iter().any(|&x| x).then(|| {
        profiles
            .iter()
            .map(|xs| {
                xs.iter()
                    .zip(selected_position.iter())
                    .filter_map(|(&x, &b)| b.then(|| x))
                    .collect()
            })
            .collect()
    })
}

fn pick_reads<'a, R: Rng>(
    asn: &[u8],
    reads: &[&'a [u8]],
    config: &ClusteringConfig,
    rng: &mut R,
) -> Vec<&'a [u8]> {
    use rand::seq::SliceRandom;
    let max = asn.iter().max().copied().unwrap() as usize + 1;
    let mut bucket = vec![vec![]; max];
    for (&asn, read) in asn.iter().zip(reads.iter()) {
        bucket[asn as usize].push(read);
    }
    bucket.iter_mut().for_each(|x| x.shuffle(rng));
    let mut counts = vec![3; max];
    for (i, b) in bucket.iter().enumerate() {
        if b.is_empty() {
            counts[i] = 0;
        }
    }
    let choises: Vec<_> = (0..max).collect();
    (0..config.subbatch_num)
        .filter_map(|_| {
            let denom: usize = counts.iter().sum();
            if denom == 0 {
                return None;
            }
            let chosen = *choises
                .choose_weighted(rng, |&a| counts[a] as f64 / denom as f64)
                .ok()?;
            let picked = bucket[chosen].pop().unwrap();
            counts[chosen] += 1;
            if bucket[chosen].is_empty() {
                counts[chosen] = 0;
            }
            Some(*picked)
        })
        .collect()
}

fn try_clustering<R: Rng>(
    reads: &[&[u8]],
    rng: &mut R,
    config: &ClusteringConfig,
) -> Option<Vec<u8>> {
    let ClusteringConfig {
        band_width,
        probe_num,
        cluster_num,
        subbatch_num,
        ..
    } = *config;
    let cons_template = kiley::consensus_bounded(&reads, rng.gen(), 3, band_width, 2 * band_width)?;
    let probes: Vec<_> = (0..probe_num)
        .filter_map(|_| {
            use rand::seq::SliceRandom;
            let picked: Vec<_> = reads.choose_multiple(rng, subbatch_num).copied().collect();
            polish_until_converge_banded(&cons_template, &picked, band_width)
        })
        .collect();
    let profiles: Vec<Vec<_>> = reads
        .iter()
        .map(|r| {
            let center = edlib_sys::global_dist(r, &cons_template) as i32;
            probes
                .iter()
                .map(|p| edlib_sys::global_dist(p, r) as i32 - center)
                .collect()
        })
        .collect();
    let selected_position: Vec<_> = (0..probes.len())
        .map(|i| {
            let first = profiles[0][i];
            profiles.iter().any(|p| p[i] != first)
        })
        .collect();
    selected_position.iter().any(|&x| x).then(|| {
        let num_probe = selected_position.iter().filter(|&&x| x).count();
        debug!("PROBE\t{}", num_probe);
        let profiles: Vec<Vec<_>> = profiles
            .iter()
            .map(|xs| {
                xs.iter()
                    .zip(selected_position.iter())
                    .filter_map(|(&x, &b)| b.then(|| x))
                    .collect()
            })
            .collect();
        for (i, prf) in profiles.iter().enumerate() {
            let prf: Vec<_> = prf.iter().map(|x| format!("{}", x)).collect();
            debug!("{}\t[{}]", i, prf.join("\t"));
        }
        kmeans(&profiles, cluster_num, rng)
    })
}

// TODO:Maybe we need more sophisticated clustering algorithm here.
// DBSCAN?
// Kmeans with random seeds?
// Or, maybe we can tune the number of clustering?
fn kmeans<R: Rng>(data: &[Vec<i32>], k: u8, rng: &mut R) -> Vec<u8> {
    use rand::seq::SliceRandom;
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
                let mut dists: Vec<_> = centers
                    .iter()
                    .enumerate()
                    .map(|(i, center)| (i as u8, euclid_norm(center, x)))
                    .collect();
                let min = dists
                    .iter()
                    .map(|x| x.1)
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap();
                dists.retain(|x| x.1 - min < 0.000001);
                assert!(!dists.is_empty());
                dists.choose(rng).map(|x| x.0)
            })
            .collect();
        assert_eq!(new_assignments.len(), assignments.len());
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

fn kmeans_f64<R: Rng>(data: &[Vec<f64>], k: u8, rng: &mut R) -> Vec<u8> {
    use rand::seq::SliceRandom;
    let mut assignments: Vec<_> = (0..data.len()).map(|_| rng.gen_range(0..k)).collect();
    loop {
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
        let new_assignments: Vec<_> = data
            .iter()
            .filter_map(|x| {
                let mut dists: Vec<_> = centers
                    .iter()
                    .enumerate()
                    .map(|(i, center)| (i as u8, euclid_norm_f64(center, x)))
                    .collect();
                let min = dists
                    .iter()
                    .map(|x| x.1)
                    .min_by(|x, y| x.partial_cmp(y).unwrap())
                    .unwrap();
                dists.retain(|x| x.1 - min < 0.000001);
                assert!(!dists.is_empty());
                dists.choose(rng).map(|x| x.0)
            })
            .collect();
        assert_eq!(new_assignments.len(), assignments.len());
        if new_assignments == assignments {
            break assignments;
        } else {
            assignments = new_assignments;
        }
    }
}

fn euclid_norm_f64(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(ys.len(), xs.len());
    xs.iter()
        .zip(ys.iter())
        .map(|(x, &y)| (x - y).powi(2))
        .sum()
}
