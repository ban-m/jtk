use clap::Parser;
use haplotyper::local_clustering::kmeans::ClusteringConfig;
use kiley::gen_seq::*;
use log::*;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
#[derive(Parser, Debug)]
#[clap(name = "CheckClustering")]
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
    /// The number of the read on each cluster.
    #[clap(short, long, default_value_t = 20)]
    coverage: usize,
    /// The error rate of the reads
    #[clap(short, long, default_value_t = 0.15)]
    error_rate: f64,
    /// The number of variants in each cluster.
    #[clap(short, long, default_value_t = 1)]
    variant_num: usize,
    /// The length of the template.
    #[clap(short, long, default_value_t = 2000)]
    template_len: usize,
    /// The band length of the alignment.
    #[clap(short, long, default_value_t = 50)]
    radius: usize,
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
    let cluster_num = command_arg.cluster_num;
    let coverage = command_arg.coverage;
    let seed = command_arg.seed;
    let band = command_arg.radius;
    let error_rate = command_arg.error_rate;
    let profile = Profile {
        sub: command_arg.error_rate / 3f64,
        del: command_arg.error_rate / 3f64,
        ins: command_arg.error_rate / 3f64,
    };
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(command_arg.seed);
    let template: Vec<_> = generate_seq(&mut rng, command_arg.template_len);
    let mut templates = vec![];
    for _ in 0..cluster_num - 1 {
        let mut seq = template.clone();
        for _ in 0..command_arg.variant_num {
            let var_pos = rng.gen_range(0..seq.len());
            let var_type = loop {
                let var_type = rng.gen::<usize>() % 9;
                if 4 <= var_type || b"ACGT"[var_type] != seq[var_pos] {
                    break var_type;
                }
            };
            match var_type {
                0..=3 => {
                    debug!("Sub\t{var_pos}");
                    seq[var_pos] = b"ACGT"[var_type];
                }
                4..=7 => {
                    debug!("Ins\t{var_pos}");
                    seq.insert(var_pos, b"ACGT"[var_type - 4]);
                }
                _ => {
                    debug!("Del\t{var_pos}");
                    seq.remove(var_pos);
                }
            };
        }
        templates.push(seq);
    }
    templates.push(template);
    let answer: Vec<_> = (0..cluster_num).flat_map(|k| vec![k; coverage]).collect();
    let reads: Vec<_> = templates
        .iter()
        .flat_map(|template| -> Vec<Vec<u8>> {
            (0..coverage)
                .map(|_| introduce_randomness(template, &mut rng, &profile))
                .collect()
        })
        .collect();
    let mut draft = kiley::ternary_consensus_by_chunk(&reads, band);
    let hmm = kiley::hmm::guided::PairHiddenMarkovModel::default();
    let mut ops: Vec<_> = reads.iter().map(|x| hmm.align(&draft, x, band).1).collect();
    draft = hmm.polish_until_converge_with(&draft, &reads, &mut ops, band);
    let gains = haplotyper::likelihood_gains::estimate_gain(&hmm, 4283094, 100, 20, 5);
    let config = ClusteringConfig::new(band, cluster_num as u8, coverage as f64, &gains);
    let strands = vec![true; reads.len()];
    let start = std::time::Instant::now();
    use haplotyper::local_clustering::kmeans;
    let clustering =
        kmeans::clustering_dev(&draft, &reads, &mut ops, &strands, &mut rng, &hmm, &config);
    let (preds, _, _, _) = clustering.unwrap();
    let end = std::time::Instant::now();
    let time = (end - start).as_millis();
    let rand_idx = haplotyper::local_clustering::rand_index(&preds, &answer);
    let adj_rand = haplotyper::local_clustering::adjusted_rand_index(&preds, &answer);
    println!("RESULT\t{seed}\t{time}\t{rand_idx}\t{adj_rand}\t{coverage}\t{error_rate}");
    Ok(())
}
