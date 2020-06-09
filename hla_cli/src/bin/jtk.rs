use clap::{App, Arg, SubCommand};
use definitions::*;
use haplotyper::*;
#[macro_use]
extern crate log;
fn subcommand_entry() -> App<'static, 'static> {
    SubCommand::with_name("entry")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Entry point. It encodes a fasta file into HLA-class file.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
}

fn subcommand_extract() -> App<'static, 'static> {
    SubCommand::with_name("extract")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Exit point. It extract fasta/q file from a HLA-class file.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("target")
                .short("t")
                .long("target")
                .takes_value(true)
                .value_name("TARGET")
                .required(true)
                .possible_values(&["raw_reads", "hic_reads", "units"]),
        )
        .arg(
            Arg::with_name("format")
                .short("f")
                .required(true)
                .takes_value(true)
                .value_name("FORMAT")
                .possible_values(&["fasta", "fastq"]),
        )
}

fn subcommand_stats() -> App<'static, 'static> {
    SubCommand::with_name("stats")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Write stats to the specified file. It passes through the stdin to the stdout")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("file")
                .long("file")
                .value_name("FILE")
                .short("f")
                .required(true)
                .takes_value(true),
        )
}

fn entry(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Entry");
    use std::io::BufReader;
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    let seqs = bio_utils::fasta::parse_into_vec_from(reader)?;
    debug!("Encoding {} reads", seqs.len());
    let dataset: DataSet = DataSet::new(seqs);
    use std::io::BufWriter;
    let stdout = std::io::stdout();
    let mut wtr = BufWriter::new(stdout.lock());
    match serde_json::ser::to_writer_pretty(&mut wtr, &dataset) {
        Err(why) => {
            eprintln!("{:?}", why);
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(()) => Ok(()),
    }
}
fn extract(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    use std::io::BufReader;
    debug!("Extract");
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    let dataset: DataSet = match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            std::process::exit(1);
        }
        Ok(res) => res,
    };
    debug!("Target is {}", matches.value_of("target").unwrap());
    debug!("Format is {}", matches.value_of("format").unwrap());
    let target = match matches.value_of("target").unwrap() {
        "raw_reads" => ExtractTarget::RawReads,
        "hic_reads" => ExtractTarget::HiCReads,
        "units" => ExtractTarget::Units,
        &_ => unreachable!(),
    };
    let stdout = std::io::stdout();
    let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
    if matches.value_of("format").unwrap() == "fasta" {
        for seq in dataset.extract_fasta(target) {
            wtr.write_record(&seq)?;
        }
    }
    Ok(())
}

fn stats(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    use std::io::BufReader;
    debug!("Stats");
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    let dataset: DataSet = match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            std::process::exit(1);
        }
        Ok(res) => res,
    };
    let wtr = std::io::BufWriter::new(std::fs::File::create(matches.value_of("file").unwrap())?);
    dataset.stats(wtr)?;
    let stdout = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(stdout.lock());
    match serde_json::ser::to_writer_pretty(&mut wtr, &dataset) {
        Err(why) => {
            eprintln!("{:?}", why);
            std::process::exit(1);
        }
        Ok(()) => Ok(()),
    }
}
fn main() -> std::io::Result<()> {
    let matches = App::new("jtk")
        .version("0.1")
        .author("Bansho Masutani")
        .about("HLA toolchain")
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .get_matches();
    match matches.subcommand() {
        ("entry", Some(sub_m)) => entry(sub_m),
        ("extract", Some(sub_m)) => extract(sub_m),
        ("stats", Some(sub_m)) => stats(sub_m),
        _ => Ok(()),
    }
}
