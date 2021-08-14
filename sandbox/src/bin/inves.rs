use std::io::{BufRead, BufReader};
fn main() -> std::io::Result<()> {
    env_logger::init();
    let args: Vec<_> = std::env::args().collect();
    let seqs = std::fs::File::open(&args[1])
        .map(BufReader::new)?
        .lines()
        .filter_map(|line| line.ok());
    let template: Vec<_> = std::fs::File::open(&args[2])
        .map(BufReader::new)?
        .lines()
        .find_map(|line| line.ok())
        .unwrap()
        .into_bytes();
    let features: Vec<_> = std::fs::File::open(&args[3])
        .map(BufReader::new)?
        .lines()
        .filter_map(|line| line.ok())
        .map(|line| {
            line.split(',')
                .map(|x| -> f64 { x.parse().unwrap() })
                .collect::<Vec<_>>()
        })
        .collect();
    let queries: Vec<_> = seqs
        .map(|line| {
            let mut line = line.split('\t');
            let name = line.next().unwrap().to_string();
            let is_dbb = line.next().unwrap() == "true";
            let asn: u8 = line.next().unwrap().parse().unwrap();
            let seq = line.next().unwrap().to_string();
            (is_dbb, asn, name, seq)
        })
        .collect();
    let (answer, queries): (Vec<_>, Vec<_>) = queries
        .into_iter()
        //.filter_map(|(is_dbb, asn, seq)| (asn != 0).then(|| (is_dbb, seq)))
        .map(|(is_dbb, _, name, seq)| ((is_dbb, name), seq))
        .unzip();
    use haplotyper::local_clustering;
    let mut config = local_clustering::kmeans::ClusteringConfig::new(100, 4, 28f64);
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128Plus;
    let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(293);
    let queries: Vec<_> = queries.iter().map(|x| x.as_bytes()).collect();
    let asn = local_clustering::kmeans::clustering_w_template(
        &template,
        &features,
        &queries,
        &mut rng,
        &mut config,
    )
    .unwrap();
    for (asn, (id, name)) in asn.iter().zip(answer.iter()) {
        println!("{}\t{}\t{}", asn, id, name);
    }
    Ok(())
}
