use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    // for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
    // let sum = logsumexp(&node.posterior);
    // node.posterior
    //     .iter_mut()
    //     .for_each(|x| *x = (*x - sum).exp());
    // assert!((1f64 - node.posterior.iter().sum::<f64>()).abs() < 0.000001);
    // }
    use haplotyper::dirichlet_mixture;
    let config = dirichlet_mixture::ClusteringConfig::default();
    use dirichlet_mixture::DirichletMixtureCorrection;
    ds.correct_clustering(&config);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}

// fn logsumexp(xs: &[f64]) -> f64 {
//     if xs.is_empty() {
//         0.
//     } else if xs.len() == 1 {
//         xs[0]
//     } else {
//         let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
//         let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
//         assert!(sum >= 0., "{:?}->{}", xs, sum);
//         max + sum
//     }
// }
