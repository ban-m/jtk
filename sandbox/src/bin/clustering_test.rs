use clap::Parser;
use haplotyper::local_clustering::kmeans::ClusteringConfig;
use rand::SeedableRng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(name = "GenMockData")]
#[clap(author = "Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>")]
#[clap(version = "1.0")]
#[clap(author,version,about,long_about=None)]
struct Args {
    /// Set the seed of a pseudorandom number generator.
    #[clap(short, long, default_value_t = 7)]
    seed: u64,
    /// Set the number of the cluster.
    #[clap(short, long, default_value_t = 2)]
    cluster_num: usize,
    /// The path to the reads.
    #[clap(short, long, parse(from_os_str))]
    reads: PathBuf,
    /// The path to the template.
    #[clap(short, long, parse(from_os_str))]
    template: PathBuf,
    #[clap(short, long, parse(from_occurrences))]
    verbose: usize,
}
fn main() -> std::io::Result<()> {
    let command_arg = Args::parse();
    let level = match command_arg.verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let dataset = bio_utils::fasta::parse_into_vec(&command_arg.reads)?;
    let mut reference = bio_utils::fasta::parse_into_vec(&command_arg.template)?;
    let reference = reference.pop().unwrap();
    let clusters = command_arg.cluster_num;
    let coverage = (dataset.len() / clusters) as f64;
    let clusters = clusters as u8;
    use haplotyper::local_clustering::kmeans;
    let seqs: Vec<_> = dataset
        .iter()
        .map(|x| match x.desc().unwrap().contains("true") {
            true => bio_utils::revcmp(x.seq()),
            false => x.seq().to_vec(),
        })
        .collect();
    let template = reference.seq();
    let mut rng: Xoroshiro128PlusPlus = SeedableRng::seed_from_u64(command_arg.seed);
    let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let radius = template.len() / 50;
    let config = ClusteringConfig::new(radius, clusters, coverage, 1.8);
    let strand: Vec<_> = (0..seqs.len()).map(|x| x % 2 == 0).collect();
    let mut ops: Vec<_> = seqs
        .iter()
        .map(|seq| kiley::bialignment::global_banded(template, seq, 1, -1, -1, -1, radius).1)
        .collect();
    let (preds, _, _, _) =
        kmeans::clustering_dev(&template, &seqs, &mut ops, &strand, &mut rng, &hmm, &config)
            .unwrap();
    for (ans, record) in preds.iter().zip(dataset.iter()) {
        println!("{}\t{ans}", record.id());
    }
    Ok(())
}
