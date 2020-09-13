use definitions::*;
fn main() -> std::io::Result<()> {
    env_logger::init();
    use std::io::BufReader;
    let args: Vec<_> = std::env::args().collect();
    let ds = BufReader::new(std::fs::File::open(&args[1])?);
    let ds: DataSet = serde_json::de::from_reader(ds).unwrap();
    use haplotyper::Encode;
    let ds = ds.encode(24);
    use haplotyper::MultiplicityEstimation;
    let config = haplotyper::MultiplicityEstimationConfig::new(4, 201, None);
    ds.estimate_multiplicity(&config);
    // Remove all clustering information.
    // for read in ds.encoded_reads.iter_mut() {
    //     for node in read.nodes.iter_mut() {
    //         node.cluster = 0;
    //     }
    // }
    // ds.assignments = ds
    //     .encoded_reads
    //     .iter()
    //     .map(|r| Assignment::new(r.id, 0))
    //     .collect();
    // use haplotyper::Assemble;
    // let config = haplotyper::AssembleConfig::default();
    // let graphs = ds.assemble_as_graph(&config);
    // println!("ID\tCoverage\tMean\tMax\tLen");
    // for graph in graphs {
    //     for node in graph.nodes.iter() {
    //         let len = node.segments.len();
    //         let unit: HashSet<_> = node.segments.iter().map(|t| t.unit).collect();
    //         let coverage = ds
    //             .encoded_reads
    //             .iter()
    //             .map(|r| r.nodes.iter().filter(|n| unit.contains(&n.unit)).count())
    //             .sum::<usize>();
    //         let mut covs: Vec<_> = unit
    //             .iter()
    //             .map(|&unit| {
    //                 ds.encoded_reads
    //                     .iter()
    //                     .filter(|r| r.nodes.iter().any(|n| n.unit == unit))
    //                     .count()
    //             })
    //             .collect();
    //         covs.sort();
    //         let med = covs.last().unwrap();
    //         let mean = coverage / len;
    //         println!("{}\t{}\t{}\t{}\t{}", node.id, coverage, mean, med, len);
    //     }
    //     let covs: Vec<_> = graph
    //         .nodes
    //         .iter()
    //         .map(|node| {
    //             let len = node.segments.len();
    //             let unit: HashSet<_> = node.segments.iter().map(|t| t.unit).collect();
    //             let coverage = ds
    //                 .encoded_reads
    //                 .iter()
    //                 .map(|r| r.nodes.iter().filter(|n| unit.contains(&n.unit)).count())
    //                 .sum::<usize>();
    //             (coverage / len) as u64
    //         })
    //         .collect();
    //     for k in 1..4 {
    //         let (asns, lk) = clustering(&covs, k);
    //         let aic = -2. * lk + 2. * k as f64;
    //         println!("{:?}\t{:?}", lk, aic);
    //         for (asn, cov) in asns.iter().zip(covs.iter()) {
    //             println!("{}\t{}", asn, cov);
    //         }
    //     }
    // }
    Ok(())
}

// struct Model {
//     cluster: usize,
//     fractions: Vec<f64>,
//     lambdas: Vec<f64>,
// }

// impl Model {
//     fn new(data: &[u64], weight: &[Vec<f64>], k: usize) -> Self {
//         let sum: Vec<_> = (0..k)
//             .map(|cl| weight.iter().map(|ws| ws[cl]).sum::<f64>())
//             .collect();
//         let fractions: Vec<_> = sum.iter().map(|w| w / data.len() as f64).collect();
//         let lambdas: Vec<_> = sum
//             .iter()
//             .enumerate()
//             .map(|(cl, sum)| {
//                 weight
//                     .iter()
//                     .zip(data)
//                     .map(|(ws, &x)| x as f64 * ws[cl])
//                     .sum::<f64>()
//                     / sum
//             })
//             .collect();
//         let cluster = k;
//         Self {
//             cluster,
//             fractions,
//             lambdas,
//         }
//     }
//     fn lk(&self, data: &[u64]) -> f64 {
//         data.iter().map(|&d| self.lk_data(d)).sum::<f64>()
//     }
//     fn lk_data(&self, data: u64) -> f64 {
//         let lks: Vec<_> = (0..self.cluster)
//             .map(|cl| {
//                 self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//                     - self.lambdas[cl]
//                     - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
//             })
//             .collect();
//         logsumexp(&lks)
//     }
//     fn new_weight(&self, data: u64) -> Vec<f64> {
//         let lks: Vec<_> = (0..self.cluster)
//             .map(|cl| {
//                 self.fractions[cl].ln() + data as f64 * self.lambdas[cl].ln()
//                     - self.lambdas[cl]
//                     - (0..data).map(|x| ((x + 1) as f64).ln()).sum::<f64>()
//             })
//             .collect();
//         let lk = logsumexp(&lks);
//         assert!((1. - lks.iter().map(|x| (x - lk).exp()).sum::<f64>()).abs() < 0.0001);
//         lks.iter().map(|x| (x - lk).exp()).collect()
//     }
// }

// fn logsumexp(xs: &[f64]) -> f64 {
//     if xs.is_empty() {
//         return 0.;
//     }
//     let max = xs.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
//     let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
//     assert!(sum >= 0., "{:?}->{}", xs, sum);
//     max + sum
// }

// fn clustering(data: &[u64], k: usize) -> (Vec<usize>, f64) {
//     let mut weight: Vec<_> = (0..data.len())
//         .map(|idx| {
//             let mut ws = vec![0.; k];
//             ws[idx % k] = 1.;
//             ws
//         })
//         .collect();
//     let mut lk = std::f64::NEG_INFINITY;
//     loop {
//         let model = Model::new(data, &weight, k);
//         let new_lk = model.lk(data);
//         let diff = new_lk - lk;
//         if diff < 0.00001 {
//             break;
//         }
//         lk = new_lk;
//         weight = data.iter().map(|&d| model.new_weight(d)).collect();
//     }
//     let asns: Vec<_> = weight
//         .iter()
//         .map(|ws| {
//             ws.iter()
//                 .enumerate()
//                 .max_by(|x, y| (x.1).partial_cmp(&y.1).unwrap())
//                 .unwrap()
//                 .0
//         })
//         .collect();
//     let model = Model::new(data, &weight, k);
//     for (f, lambda) in model.fractions.iter().zip(model.lambdas.iter()) {
//         println!("{}\t{}", f, lambda);
//     }
//     (asns, lk)
// }
