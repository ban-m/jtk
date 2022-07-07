use clap::Parser;
use kiley::gen_seq;
use rand::Rng;
use rand_xoshiro::Xoroshiro128PlusPlus;
use std::{io::BufWriter, path::PathBuf};

#[derive(Parser, Debug)]
#[clap(name = "GenMockData")]
#[clap(author = "Bansho Masutani<ban-m@g.ecc.u-tokyo.ac.jp>")]
#[clap(version = "1.0")]
#[clap(author,version,about,long_about=None)]
struct Args {
    /// Set the seed of a pseudorandom number generator.
    #[clap(short, long, default_value_t = 7)]
    seed: u64,
    /// Set the number of reads to be simulated.
    #[clap(short, long, default_value_t = 40)]
    read_num: usize,
    /// Set the number of the cluster.
    #[clap(short, long, default_value_t = 2)]
    cluster_num: usize,
    /// Set the number of variants in each cluster.
    #[clap(short, long, default_value_t = 1)]
    variant_num: usize,
    #[clap(short, long, default_value_t = 2000)]
    /// Set the length of the template
    length: usize,
    #[clap(short, long, default_value_t = 0.15)]
    /// Set the error rate.
    error_rate: f64,
    #[clap(short, long, parse(from_os_str), default_value = "./")]
    output_dir: PathBuf,
}
fn main() -> std::io::Result<()> {
    let command_arg = Args::parse();
    env_logger::init();
    std::fs::create_dir_all(&command_arg.output_dir)?;
    let profile = gen_seq::Profile {
        sub: command_arg.error_rate / 3f64,
        del: command_arg.error_rate / 3f64,
        ins: command_arg.error_rate / 3f64,
    };
    let mut rng: Xoroshiro128PlusPlus = rand::SeedableRng::seed_from_u64(command_arg.seed);
    let template: Vec<_> = gen_seq::generate_seq(&mut rng, command_arg.length);
    let mut templates = vec![template.clone()];
    for _ in 0..command_arg.cluster_num - 1 {
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
                    seq[var_pos] = b"ACGT"[var_type];
                }
                4..=7 => {
                    seq.insert(var_pos, b"ACGT"[var_type - 4]);
                }
                _ => {
                    seq.remove(var_pos);
                }
            };
        }
        templates.push(seq);
    }
    use sandbox::generate_test_data;
    let probs: Vec<_> = vec![(command_arg.cluster_num as f64).recip(); command_arg.cluster_num];
    let (dataset, answer) =
        generate_test_data(&templates, command_arg.read_num, &mut rng, &probs, &profile);
    // Write templates.
    use std::io::Write;
    for (i, template) in templates.iter().enumerate() {
        let mut path = command_arg.output_dir.clone();
        path.push(&format!("{}.fa", i));
        let mut wtr = std::fs::File::create(path).map(BufWriter::new)?;
        let seq = std::str::from_utf8(&template).unwrap();
        writeln!(&mut wtr, ">{i}\n{seq}")?;
    }
    // Write reads.
    let mut path = command_arg.output_dir.clone();
    path.push(&format!("reads.fa"));
    let mut wtr = std::fs::File::create(path).map(BufWriter::new)?;
    for (i, (seq, answer)) in dataset.iter().zip(answer.iter()).enumerate() {
        let seq = match i % 2 == 0 {
            true => seq.clone(),
            false => bio_utils::revcmp(seq),
        };
        let seq = std::str::from_utf8(&seq).unwrap();
        let is_rev = i % 2 == 1;
        writeln!(&mut wtr, ">{i}_{answer} is_rev={is_rev}\n{seq}")?;
    }
    Ok(())
}
