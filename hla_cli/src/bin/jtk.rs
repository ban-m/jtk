use bio_utils::lasttab::LastTAB;
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
                .long("format")
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

fn subcommand_select_unit() -> App<'static, 'static> {
    SubCommand::with_name("select_unit")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Selecting units")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("chunk_len")
                .short("l")
                .long("chunk_len")
                .takes_value(true)
                .default_value(&"1000")
                .help("Length of a chunk"),
        )
        .arg(
            Arg::with_name("chunk_num")
                .short("n")
                .long("chunk_num")
                .takes_value(true)
                .default_value(&"500")
                .help("Number of chunks"),
        )
        .arg(
            Arg::with_name("skip_len")
                .short("s")
                .long("skip_len")
                .takes_value(true)
                .default_value(&"4000")
                .help("Margin between units"),
        )
        .arg(
            Arg::with_name("margin")
                .short("m")
                .long("margin")
                .takes_value(true)
                .default_value(&"2000")
                .help("Margin at the both end of a read."),
        )
}

fn subcommand_encode() -> App<'static, 'static> {
    SubCommand::with_name("encode")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Encode reads by alignments.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("alignment")
                .short("a")
                .takes_value(true)
                .required(true)
                .long("alignment")
                .value_name("ALIGNMENT<LastTAB>")
                .help("alignment between units and reads (units are references)."),
        )
}

fn subcommand_view() -> App<'static, 'static> {
    SubCommand::with_name("view")
        .version("0.1")
        .author("BanshoMasutani")
        .about("View reads")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("name")
                .short("n")
                .long("name")
                .takes_value(true)
                .required(true)
                .value_name("NAME")
                .help("Name of the read to be showed."),
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

fn select_unit(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    use std::io::BufReader;
    debug!("Select Units");
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
    let chunk_num: usize = matches
        .value_of("chunk_num")
        .and_then(|e| e.parse().ok())
        .expect("Chunk Num");
    let chunk_len: usize = matches
        .value_of("chunk_len")
        .and_then(|e| e.parse().ok())
        .expect("Chunk len");
    let margin: usize = matches
        .value_of("margin")
        .and_then(|e| e.parse().ok())
        .expect("Margin");
    let skip_len: usize = matches
        .value_of("skip_len")
        .and_then(|e| e.parse().ok())
        .expect("Skip Len");
    let config = UnitConfig {
        chunk_len,
        chunk_num,
        margin,
        skip_len,
    };
    let dataset = dataset.select_chunks(&config);
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

fn parse_tab_file<P: AsRef<std::path::Path>>(tab_file: P) -> std::io::Result<Vec<LastTAB>> {
    let lines = std::fs::read_to_string(tab_file)?;
    Ok(lines
        .lines()
        .filter(|e| !e.starts_with('#'))
        .filter_map(LastTAB::from_line)
        .collect())
}

fn encode(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Encode");
    use std::io::BufReader;
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
    let alignment = parse_tab_file(matches.value_of("alignment").unwrap())?;
    let dataset = dataset.encode(&alignment);
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

fn view(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("View");
    use std::io::BufReader;
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
    let name = matches.value_of("name").unwrap();
    dataset.view(name);
    Ok(())
}

fn main() -> std::io::Result<()> {
    let matches = App::new("jtk")
        .version("0.1")
        .author("Bansho Masutani")
        .about("HLA toolchain")
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .subcommand(subcommand_select_unit())
        .subcommand(subcommand_encode())
        .subcommand(subcommand_view())
        .get_matches();
    match matches.subcommand() {
        ("entry", Some(sub_m)) => entry(sub_m),
        ("extract", Some(sub_m)) => extract(sub_m),
        ("stats", Some(sub_m)) => stats(sub_m),
        ("select_unit", Some(sub_m)) => select_unit(sub_m),
        ("encode", Some(sub_m)) => encode(sub_m),
        ("view", Some(sub_m)) => view(sub_m),
        _ => Ok(()),
    }
}
