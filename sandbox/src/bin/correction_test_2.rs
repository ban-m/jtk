use definitions::*;
use std::io::BufReader;
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let mut ds: DataSet =
        serde_json::de::from_reader(BufReader::new(std::fs::File::open(&args[1]).unwrap()))
            .unwrap();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        let sum = logsumexp(&node.posterior);
        node.posterior.iter_mut().for_each(|x| *x = *x - sum);
    }
    // use haplotyper::dirichlet_correction::*;
    // let config = Config::default();
    // ds.correct_clustering(&config);
    use haplotyper::read_clustering::*;
    let config = ReadClusteringConfig::default();
    ds.read_clustering(&config);
    println!("{}", serde_json::ser::to_string(&ds).unwrap());
    Ok(())
}

fn logsumexp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        0.
    } else if xs.len() == 1 {
        xs[0]
    } else {
        let max = xs.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        let sum = xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln();
        assert!(sum >= 0., "{:?}->{}", xs, sum);
        max + sum
    }
}
