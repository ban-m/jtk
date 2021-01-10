//! A small K-means clustering algorithm.
//! It run several times, pick the best one.
use rand::Rng;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
pub fn kmeans(data: &[Vec<f64>], k: usize, seed: Option<u64>) -> Vec<usize> {
    let data: Vec<_> = data.iter().map(|x| x.as_slice()).collect();
    kmeans_slice(&data, k, seed)
}

pub fn kmeans_slice(data: &[&[f64]], k: usize, seed: Option<u64>) -> Vec<usize> {
    let mut rng: Xoshiro256StarStar = match seed {
        Some(seed) => SeedableRng::seed_from_u64(seed),
        None => {
            let mut rng = rand::thread_rng();
            SeedableRng::seed_from_u64(rng.gen::<u64>())
        }
    };
    kmeans_with_rng(data, k, &mut rng)
}

fn kmeans_with_rng<R: Rng>(data: &[&[f64]], k: usize, rng: &mut R) -> Vec<usize> {
    let mut assignments: Vec<_> = data.iter().map(|_| rng.gen::<usize>() % k).collect();
    let len = data[0].len();
    for d in data.iter() {
        assert_eq!(d.len(), len);
    }
    let mut residual = std::f64::INFINITY;
    let mut diff = std::f64::INFINITY;
    while diff > 0.001 {
        let mut centers: Vec<_> = vec![vec![0.; len]; k];
        let mut counts: Vec<_> = vec![0; k];
        for (&asn, xs) in assignments.iter().zip(data.iter()) {
            for (i, x) in xs.iter().enumerate() {
                centers[asn][i] += x;
            }
            counts[asn] += 1;
        }
        centers.iter_mut().zip(counts.iter()).for_each(|(xs, &c)| {
            xs.iter_mut().for_each(|x| *x /= c as f64);
        });
        let mut next_residual = 0.;
        for (asn, xs) in assignments.iter_mut().zip(data.iter()) {
            let (new_asn, min_distance) = centers
                .iter()
                .map(|cs| {
                    cs.iter()
                        .zip(xs.iter())
                        .map(|(c, x)| (c - x).powi(2i32))
                        .sum::<f64>()
                })
                .enumerate()
                .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
                .unwrap();
            *asn = new_asn;
            next_residual += min_distance;
        }
        diff = residual - next_residual;
        assert!(
            diff.is_sign_positive(),
            "{},{},{}",
            residual,
            next_residual,
            diff
        );
        residual = next_residual;
    }
    assignments
}

#[cfg(test)]
mod test {
    use super::*;
    use rand_distr::{Distribution, Normal};
    #[test]
    fn works() {}
    #[test]
    fn k_means_test() {
        let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(239);
        let data_num = 100;
        let cluster = 4;
        let models: Vec<Vec<_>> = vec![
            vec![(0., 0.5), (0., 0.5), (0., 1.), (2., 1.)],
            vec![(2., 0.5), (0., 0.5), (2., 1.), (0., 1.)],
            vec![(0., 0.5), (2., 0.5), (0., 1.), (0., 1.)],
            vec![(1., 0.5), (0., 0.5), (0., 1.), (-3., 1.)],
        ]
        .into_iter()
        .map(|params| {
            params
                .into_iter()
                .map(|(mean, sd)| Normal::new(mean, sd).unwrap())
                .collect()
        })
        .collect();
        let data: Vec<Vec<f64>> = (0..cluster)
            .flat_map(|k| {
                (0..data_num)
                    .map(|_| models[k].iter().map(|m| m.sample(&mut rng)).collect())
                    .collect::<Vec<_>>()
            })
            .collect();
        let asn = kmeans(&data, cluster, None);
        let total_missed = asn
            .chunks(data_num)
            .map(|ax| {
                let repr = ax[0];
                ax.iter().filter(|&&x| x != repr).count()
            })
            .sum::<usize>();
        let error = total_missed as f64 / asn.len() as f64;
        if error > 0.05 {
            eprintln!("{}\t{}\t{}", total_missed, asn.len(), error);
            // for (i, a) in asn.iter().enumerate() {
            //     eprintln!("{}\t{}", i, a,);
            // }
            panic!();
        }
    }
}
