fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for line in File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|e| e.ok())
    {
        let line: Vec<_> = line.split('\t').collect();
        let read_id: usize = line[1].parse().unwrap();
        let pos: usize = line[2].parse().unwrap();
        let lk1: f64 = line[3].parse().unwrap();
        let lk2: f64 = line[4].parse().unwrap();
        let diff = lk1 - lk2;
        buckets.entry(pos).or_default().push((read_id, diff))
    }
    let mut buckets: Vec<_> = buckets.into_iter().collect();
    buckets.sort_by_key(|x| x.0);
    for (pos, data) in buckets {
        let (assignment, lk) = clustering(&data, pos as u64, 2);
        let (_, lk_0) = clustering(&data, pos as u64, 1);
        let cluster = {
            let mut slot = vec![0; 2];
            for &asn in assignment.iter() {
                slot[asn] += 1;
            }
            slot
        };
        let gain = lk - lk_0 - 3.;
        eprintln!(
            "{}\t{}\t{}\t{}\t{:?}",
            pos,
            lk - 6.,
            lk_0 - 3.,
            gain,
            cluster
        );
        for ((_, x), asn) in data.iter().zip(assignment) {
            println!("{}\t{}\t{}", pos, asn, x);
        }
    }
    //     use rand_distr::Distribution;
    // use rand_distr::Normal;
    // let mut rng: Xoshiro256StarStar = SeedableRng::seed_from_u64(239812903);
    // let normal1 = Normal::new(10., 5.).unwrap();
    // let normal2 = Normal::new(-10., 2.).unwrap();
    // let data: Vec<_> = (0..1000)
    //     .map(|i| {
    //         let x = if i < 500 {
    //             normal1.sample(&mut rng)
    //         } else {
    //             normal2.sample(&mut rng)
    //         };
    //         (0, x)
    //     })
    //     .collect();
    // let (result, lk) = (0..10u64)
    //     .map(|seed| clustering(&data, seed))
    //     .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
    //     .unwrap();
    // eprintln!("MAXLK:{}", lk);
    // for ((_, x), asn) in data.iter().zip(result) {
    //     println!("{}\t{}", asn, x);
    // }
    Ok(())
}

struct Gaussian {
    mean: f64,
    variance: f64,
    weight: f64,
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
        // let variance = xs
        //     .iter()
        //     .zip(weights.iter())
        //     .map(|(x, w)| (x - mean).powi(2i32) * w)
        //     .sum::<f64>()
        //     / sum;
        let variance = 0.05;
        let weight = sum / len;
        // eprintln!("{}\t{}\t{}", mean, variance, weight);
        Self {
            mean,
            variance,
            weight,
        }
    }
}

fn clustering(data: &[(usize, f64)], seed: u64, k: usize) -> (Vec<usize>, f64) {
    let xs: Vec<_> = data.iter().map(|x| x.1).collect();
    let variance = {
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        xs.iter().map(|x| (x - mean).powi(2i32)).sum::<f64>() / xs.len() as f64
    };
    if variance < 0.0001 {
        return (vec![0; xs.len()], 0.);
    }
    use rand::Rng;
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
    // for (cl, g) in gaussian.iter().enumerate() {
    //     eprintln!("Model:{}\t{}\t{}", cl, g.mean, g.variance);
    // }
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
    let max = xs
        .iter()
        .max_by(|x, y| x.partial_cmp(&y).unwrap_or_else(|| panic!("{},{}", x, y)))
        .unwrap_or_else(|| panic!("{:?}", xs));
    let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
    assert!(sum >= 0., "{:?}->{}", xs, sum);
    max + sum
}
