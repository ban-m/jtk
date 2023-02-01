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
    /// The number of variants in each cluster.
    #[clap(long, default_value_t = 1)]
    variant_num: usize,
    /// The band length of the alignment.
    #[clap(short, long, default_value_t = 20)]
    radius: usize,
    #[clap(short, long, parse(from_occurrences))]
    verbose: usize,
}

const ERROR_RATE: f64 = 0.10;
const TEMPLATE_LEN: usize = 1_000;
fn main() -> std::io::Result<()> {
    let command_arg = Args::parse();
    let level = match command_arg.verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .unwrap();
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let cluster_num = command_arg.cluster_num;
    let coverage = command_arg.coverage;
    let seed = command_arg.seed;
    let band = command_arg.radius;
    let profile = Profile {
        sub: ERROR_RATE / 3f64,
        del: ERROR_RATE / 3f64,
        ins: ERROR_RATE / 3f64,
    };
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(command_arg.seed);
    let template: Vec<_> = generate_seq(&mut rng, TEMPLATE_LEN);
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
    templates.push(template.clone());
    let reads: Vec<_> = templates
        .iter()
        .flat_map(|template| -> Vec<Vec<u8>> {
            (0..coverage)
                .map(|_| introduce_randomness(template, &mut rng, &profile))
                .collect()
        })
        .collect();
    let hmm = kiley::hmm::PairHiddenMarkovModel::default();
    let template = hmm.polish_until_converge(&template, &reads, command_arg.radius);
    let gains = haplotyper::likelihood_gains::estimate_gain(&hmm, 4283094, 100, 20, 5);
    let coverage = coverage as f64;
    let config = ClusteringConfig::new(band, cluster_num, coverage, coverage, &gains);
    let strands = vec![true; reads.len()];
    use haplotyper::local_clustering::kmeans;
    let ops: Vec<_> = reads
        .iter()
        .map(|x| hmm.align_guided(&template, x, band).1)
        .collect();
    let feature_vectors =
        kmeans::search_variants(&template, &reads, &ops, &strands, &mut rng, &hmm, &config);
    let mcmc_start = std::time::Instant::now();
    let mcmc_score = kmeans::cluster_filtered_variants(&feature_vectors, &config, &mut rng);
    let mcmc_score = mcmc_score.2;
    let mcmc_time = (std::time::Instant::now() - mcmc_start).as_millis();
    let exact_start = std::time::Instant::now();
    let exact_score = kmeans::cluster_filtered_variants_exact(&feature_vectors, &config);
    let exact_score = exact_score.2;
    let exact_time = (std::time::Instant::now() - exact_start).as_millis();
    let var_num = command_arg.variant_num;
    print!("RESULT\t{seed}\t{coverage}\t{var_num}\t{cluster_num}");
    println!("\t{mcmc_score}\t{exact_score}\t{mcmc_time}\t{exact_time}");
    Ok(())
}
