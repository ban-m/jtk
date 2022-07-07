use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let (answer, preds): (Vec<_>, Vec<_>) = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|line| line.ok())
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let mut line = line.split('\t');
            let answer: usize = line.next().and_then(|x| x.parse().ok()).unwrap();
            let phase: usize = line.next().and_then(|x| x.parse().ok()).unwrap();
            (answer, phase)
        })
        .unzip();
    let rand_idx = haplotyper::local_clustering::rand_index(&answer, &preds);
    let accuracy = accuracy(&answer, &preds);
    let adj_rand_idx = haplotyper::local_clustering::adjusted_rand_index(&answer, &preds);
    println!("{rand_idx}\t{accuracy}\t{adj_rand_idx}");
    Ok(())
}

fn accuracy(label: &[usize], pred: &[usize]) -> f64 {
    let matches = label.iter().zip(pred.iter()).filter(|(x, y)| x == y);
    let num_matches = matches.count();
    num_matches.max(label.len() - num_matches) as f64 / label.len() as f64
}
