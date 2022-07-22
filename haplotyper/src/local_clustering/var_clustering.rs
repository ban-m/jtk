pub fn clustering_variant_to(vars: &[Vec<i8>], k: usize, seed: u64) -> (Vec<u64>, f64) {
    let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(seed + 390240);
    let var_size = vars[0].len();
    if var_size == 1 {
        return (brute_split(vars), 1.);
    } else if var_size == 0 {
        return (vec![0; vars.len()], 1.);
    }
    let mut model = Model::new(vars, &mut rng, var_size, k);
    let mut lk = std::f64::NEG_INFINITY;
    loop {
        model.update(vars);
        let new_lk = model.lk(vars);
        if new_lk - lk < 0.0001 {
            break;
        }
        lk = new_lk;
    }
    let assignments: Vec<_> = vars.iter().map(|d| model.assign(d) as u64).collect();
    let num_params = var_size * k + (k - 1);
    (assignments, lk - num_params as f64)
}

fn brute_split(vars: &[Vec<i8>]) -> Vec<u64> {
    vars.iter()
        .enumerate()
        .map(|(i, vs)| {
            assert_eq!(vs.len(), 1);
            match vs[0] {
                1 => 1u64,
                0 => 0u64,
                _ => (i % 2) as u64,
            }
        })
        .collect()
}

use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
#[derive(Debug)]
struct Model {
    cluster_num: usize,
    var_size: usize,
    weights: Vec<Vec<f64>>,
    // Cluster -> Position -> Fraction of variant
    params: Vec<Vec<f64>>,
    fractions: Vec<f64>,
}

const SMALL: f64 = 0.0000000001;
impl Model {
    pub fn new<R: rand::Rng>(
        data: &[Vec<i8>],
        rng: &mut R,
        var_size: usize,
        cluster_num: usize,
    ) -> Self {
        let weights: Vec<_> = data
            .iter()
            .map(|_| {
                let mut ws = vec![0.; cluster_num];
                ws[rng.gen::<usize>() % cluster_num] = 1.;
                ws
            })
            .collect();
        let params: Vec<Vec<_>> = Self::calc_params(data, &weights, cluster_num, var_size);
        let fractions: Vec<_> = Self::calc_fractions(&weights, cluster_num);
        Self {
            cluster_num,
            weights,
            var_size,
            fractions,
            params,
        }
    }
    fn calc_fractions(weights: &[Vec<f64>], cluster_num: usize) -> Vec<f64> {
        let mut fractions = vec![SMALL; cluster_num];
        for ws in weights.iter() {
            assert!((1. - ws.iter().sum::<f64>()).abs() < 0.001);
            for (i, &w) in ws.iter().enumerate() {
                fractions[i] += w;
            }
        }
        fractions
            .iter_mut()
            .for_each(|x| *x /= weights.len() as f64);
        fractions
    }
    fn calc_params(
        data: &[Vec<i8>],
        weights: &[Vec<f64>],
        cluster_num: usize,
        var_size: usize,
    ) -> Vec<Vec<f64>> {
        (0..cluster_num)
            .map(|cl| {
                let mut total = vec![SMALL; var_size];
                let mut count = vec![SMALL; var_size];
                for (v, w) in data.iter().zip(weights.iter()) {
                    let w = w[cl];
                    for (pos, base) in v.iter().enumerate() {
                        match base {
                            1 => {
                                count[pos] += w;
                                total[pos] += w;
                            }
                            0 => total[pos] += w,
                            _ => {}
                        }
                    }
                }
                count.into_iter().zip(total).map(|(c, t)| c / t).collect()
            })
            .collect()
    }
    // Return log lk.
    fn calc_lks(&self, datum: &[i8]) -> Vec<f64> {
        self.fractions
            .iter()
            .zip(self.params.iter())
            .map(|(f, param)| {
                f.ln()
                    + datum
                        .iter()
                        .zip(param)
                        .filter(|&(&x, _)| x != -1)
                        .map(|(x, p)| match x {
                            1 => p.ln(),
                            0 => (1. - p).ln(),
                            _ => unreachable!(),
                        })
                        .sum::<f64>()
            })
            .collect()
    }
    fn assign(&self, datum: &[i8]) -> usize {
        self.calc_lks(datum)
            .iter()
            .enumerate()
            .max_by(|x, y| (x.1).partial_cmp(y.1).unwrap())
            .unwrap()
            .0
    }
    fn update(&mut self, data: &[Vec<i8>]) {
        self.weights = data
            .iter()
            .map(|datum| {
                let lks = self.calc_lks(datum);
                let total = crate::misc::logsumexp(&lks);
                let weights: Vec<_> = lks.iter().map(|lk| (lk - total).exp()).collect();
                assert!((1. - weights.iter().sum::<f64>()).abs() < 0.001);
                weights
            })
            .collect();
        self.params = Self::calc_params(data, &self.weights, self.cluster_num, self.var_size);
        self.fractions = Self::calc_fractions(&self.weights, self.cluster_num);
    }
    fn lk(&self, data: &[Vec<i8>]) -> f64 {
        data.iter()
            .map(|d| crate::misc::logsumexp(&self.calc_lks(d)))
            .sum()
    }
}
