// use bio_utils::lasttab::LastTAB;
use clap::{App, Arg, SubCommand};
use definitions::*;
use haplotyper::*;
use std::io::Write;
use std::io::{BufReader, BufWriter};
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
    let targets = ["raw_reads", "hic_reads", "units", "assignments"];
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
                .possible_values(&targets),
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
                .default_value(&"2000")
                .help("Length of a chunk"),
        )
        .arg(
            Arg::with_name("skip_len")
                .short("s")
                .long("skip_len")
                .takes_value(true)
                .default_value(&"2000")
                .help("Margin between units"),
        )
        .arg(
            Arg::with_name("margin")
                .short("m")
                .long("margin")
                .takes_value(true)
                .default_value(&"500")
                .help("Margin at the both end of a read."),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value(&"1")
                .help("number of threads"),
        )
        .arg(
            Arg::with_name("read_type")
                .long("read_type")
                .takes_value(true)
                .default_value(&"CCS")
                .possible_values(&["CCS", "CLR", "ONT"])
                .help("Read type. CCS, CLR, or ONT."),
        )
}

fn subcommand_encode() -> App<'static, 'static> {
    SubCommand::with_name("encode")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Encode reads by alignments (Internally invoke `LAST` tools).")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .short("t")
                .help("Number of threads")
                .takes_value(true)
                .default_value(&"1"),
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
        .arg(
            Arg::with_name("type")
                .short("t")
                .long("type")
                .takes_value(true)
                .required(true)
                .value_name("TYPE")
                .possible_values(&["read", "unit"]),
        )
}

fn subcommand_local_clustering() -> App<'static, 'static> {
    SubCommand::with_name("local_clustering")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Clustering reads. (Local)")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("limit")
                .short("l")
                .long("limit")
                .required(false)
                .value_name("LIMIT")
                .help("Maximum Execution time(sec)")
                .default_value(&"3000")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cluster_num")
                .short("c")
                .long("cluster_num")
                .required(false)
                .value_name("CLUSTER_NUM")
                .help("Minimum cluster number.")
                .default_value(&"2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("subchunk_len")
                .short("s")
                .long("subchunk_len")
                .required(false)
                .value_name("SubChunkLength")
                .help("The length of sub-chunks")
                .default_value(&"100")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("read_type")
                .long("read_type")
                .takes_value(true)
                .default_value(&"CCS")
                .help("Read type. CCS, CLR, or ONT."),
        )
}

fn subcommand_global_clustering() -> App<'static, 'static> {
    SubCommand::with_name("global_clustering")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Clustering reads (Global).")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("k")
                .short("k")
                .long("kmer_size")
                .required(false)
                .value_name("KMER_SIZE")
                .help("The size of the kmer")
                .default_value(&"3")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_cluster_size")
                .short("m")
                .long("min_cluster_size")
                .required(false)
                .value_name("MIN_CLUSTER_SIZE")
                .help("The minimum size of a cluster")
                .default_value(&"50")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mat_score")
                .short("p")
                .long("match_score")
                .required(false)
                .value_name("MATCH_SCORE")
                .help("The match score")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mismat_score")
                .short("q")
                .long("mismatch_score")
                .required(false)
                .value_name("MISMATCH_SCORE")
                .help("The mismatch score")
                .default_value(&"-1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gap_score")
                .short("g")
                .long("gap_score")
                .required(false)
                .value_name("GAP_SCORE")
                .help("The gap penalty(< 0)")
                .default_value(&"-2")
                .takes_value(true),
        )
}

fn subcommand_polish_clustering() -> App<'static, 'static> {
    SubCommand::with_name("polish_clustering")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Polish local clustering.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mat_score")
                .short("p")
                .long("match_score")
                .required(false)
                .value_name("MATCH_SCORE")
                .help("The match score")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mismat_score")
                .short("q")
                .long("mismatch_score")
                .required(false)
                .value_name("MISMATCH_SCORE")
                .help("The mismatch score")
                .default_value(&"-1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gap_score")
                .short("g")
                .long("gap_score")
                .required(false)
                .value_name("GAP_SCORE")
                .help("The gap penalty(< 0)")
                .default_value(&"-2")
                .takes_value(true),
        )
}

fn subcommand_assembly() -> App<'static, 'static> {
    SubCommand::with_name("assemble")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Assemble reads.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
}

fn subcommand_pipeline() -> App<'static, 'static> {
    SubCommand::with_name("pipeline")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Exec JTK pipeline.")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::with_name("threads")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("chunk_len")
                .long("chunk_len")
                .takes_value(true)
                .default_value(&"2000")
                .help("Length of a chunk"),
        )
        .arg(
            Arg::with_name("skip_len")
                .long("skip_len")
                .takes_value(true)
                .default_value(&"2000")
                .help("Margin between units"),
        )
        .arg(
            Arg::with_name("margin")
                .long("margin")
                .takes_value(true)
                .default_value(&"500")
                .help("Margin at the both end of a read."),
        )
        .arg(
            Arg::with_name("read_type")
                .long("read_type")
                .takes_value(true)
                .default_value(&"CCS")
                .possible_values(&["CCS", "CLR", "ONT"])
                .help("Read type. CCS, CLR, or ONT."),
        )
        .arg(
            Arg::with_name("limit")
                .long("limit")
                .required(false)
                .value_name("LIMIT")
                .help("Maximum Execution time(sec)")
                .default_value(&"3000")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("cluster_num")
                .long("cluster_num")
                .required(false)
                .value_name("CLUSTER_NUM")
                .help("Minimum cluster number.")
                .default_value(&"2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("subchunk_len")
                .long("subchunk_len")
                .required(false)
                .value_name("SubChunkLength")
                .help("The length of sub-chunks")
                .default_value(&"100")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("k")
                .long("kmer_size")
                .required(false)
                .value_name("KMER_SIZE")
                .help("The size of the kmer")
                .default_value(&"3")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_cluster_size")
                .long("min_cluster_size")
                .required(false)
                .value_name("MIN_CLUSTER_SIZE")
                .help("The minimum size of a cluster")
                .default_value(&"50")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mat_score")
                .long("match_score")
                .required(false)
                .value_name("MATCH_SCORE")
                .help("The match score")
                .default_value(&"1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mismat_score")
                .long("mismatch_score")
                .required(false)
                .value_name("MISMATCH_SCORE")
                .help("The mismatch score")
                .default_value(&"-1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gap_score")
                .long("gap_score")
                .required(false)
                .value_name("GAP_SCORE")
                .help("The gap penalty(< 0)")
                .default_value(&"-2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gfa")
                .long("gfa_output_path")
                .required(true)
                .value_name("PATH")
                .help("path to the output gfa file")
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
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    let seqs = bio_utils::fasta::parse_into_vec_from(reader)?;
    debug!("Encoding {} reads", seqs.len());
    let dataset: DataSet = DataSet::new(seqs);
    flush_file(&dataset)
}

fn extract(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Extract");
    let dataset = get_input_file()?;
    debug!("Target is {}", matches.value_of("target").unwrap());
    match matches.value_of("target").unwrap() {
        "raw_reads" => {
            let stdout = std::io::stdout();
            let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
            for seq in dataset.extract_fasta(ExtractTarget::RawReads) {
                wtr.write_record(&seq)?;
            }
        }
        "hic_reads" => {
            let stdout = std::io::stdout();
            let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
            for seq in dataset.extract_fasta(ExtractTarget::HiCReads) {
                wtr.write_record(&seq)?;
            }
        }
        "units" => {
            let stdout = std::io::stdout();
            let mut wtr = bio_utils::fasta::Writer::new(stdout.lock());
            for seq in dataset.extract_fasta(ExtractTarget::Units) {
                wtr.write_record(&seq)?;
            }
        }
        "assignments" => {
            let asn_name_desc = dataset.extract_assignments();
            let stdout = std::io::stdout();
            let mut wtr = std::io::BufWriter::new(stdout.lock());
            for (asn, name, desc) in asn_name_desc {
                writeln!(&mut wtr, "{}\t{}\t{}", asn, name, desc)?;
            }
        }
        &_ => unreachable!(),
    };
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
    debug!("Start Stats step");
    let dataset = get_input_file()?;
    let wtr = std::io::BufWriter::new(std::fs::File::create(matches.value_of("file").unwrap())?);
    dataset.stats(wtr)?;
    flush_file(&dataset)
}

fn select_unit(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Selecting Units");
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
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse().ok())
        .expect("threads");
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let config = match matches.value_of("read_type").unwrap() {
        "CCS" => UnitConfig::new_ccs(chunk_len, skip_len, margin),
        "CLR" => UnitConfig::new_clr(chunk_len, skip_len, margin),
        "ONT" => UnitConfig::new_ont(chunk_len, skip_len, margin),
        _ => unreachable!(),
    };
    let dataset = get_input_file()?.select_chunks_freq(&config);
    flush_file(&dataset)
}

fn encode(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Encoding step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let dataset = get_input_file()?.encode(threads);
    flush_file(&dataset)
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
    let dataset = get_input_file()?;
    let name = matches.value_of("name").unwrap();
    match matches.value_of("type") {
        Some(x) if x == "read" => dataset.view(name),
        Some(x) if x == "unit" => dataset.view_unit(name),
        _ => Some(()),
    };
    Ok(())
}

fn local_clustering(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Local Clustering step");
    let cluster_num: usize = matches
        .value_of("cluster_num")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let limit: u64 = matches
        .value_of("limit")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let length: usize = matches
        .value_of("subchunk_len")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let dataset = get_input_file()?;
    let config = match matches.value_of("read_type").unwrap() {
        "CCS" => ClusteringConfig::ccs(&dataset, cluster_num, length, limit),
        "CLR" => ClusteringConfig::clr(&dataset, cluster_num, length, limit),
        "ONT" => ClusteringConfig::ont(&dataset, cluster_num, length, limit),
        _ => unreachable!(),
    };
    let dataset = dataset.local_clustering(&config);
    flush_file(&dataset)
}

fn global_clustering(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Global Clustering step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let kmer: usize = matches
        .value_of("k")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let min_cluster_size = matches
        .value_of("min_cluster_size")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mat_score: i32 = matches
        .value_of("mat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mismat_score: i32 = matches
        .value_of("mismat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let gap_score: i32 = matches
        .value_of("gap_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let config = haplotyper::GlobalClusteringConfig::new(
        kmer,
        min_cluster_size,
        mat_score,
        mismat_score,
        gap_score,
    );
    let dataset = get_input_file()?.global_clustering(&config);
    flush_file(&dataset)
}

fn polish_clustering(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Polish Clustering step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mat_score: i32 = matches
        .value_of("mat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mismat_score: i32 = matches
        .value_of("mismat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let gap_score: i32 = matches
        .value_of("gap_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let config = haplotyper::PolishClusteringConfig::new(mat_score, mismat_score, gap_score);
    let dataset = get_input_file()?.polish_clustering(&config);
    flush_file(&dataset)
}

fn assembly(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let level = match matches.occurrences_of("verbose") {
        0 => "warn",
        1 => "info",
        2 => "debug",
        3 | _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    debug!("Start Assembly step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let config = AssembleConfig::default();
    let gfa = get_input_file()?.assemble_as_gfa(&config);
    let stdout = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(stdout.lock());
    writeln!(&mut wtr, "{}", gfa)
}

fn pipeline(matches: &clap::ArgMatches) -> std::io::Result<()> {
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    let seqs = bio_utils::fasta::parse_into_vec_from(reader)?;
    debug!("Encoding {} reads", seqs.len());
    let dataset: DataSet = DataSet::new(seqs);
    debug!("Start Selecting Units");
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
    let unit_config = match matches.value_of("read_type").unwrap() {
        "CCS" => UnitConfig::new_ccs(chunk_len, skip_len, margin),
        "CLR" => UnitConfig::new_clr(chunk_len, skip_len, margin),
        "ONT" => UnitConfig::new_ont(chunk_len, skip_len, margin),
        _ => unreachable!(),
    };
    let dataset = dataset.select_chunks_freq(&unit_config).encode(threads);
    debug!("Start Local Clustering step");
    let cluster_num: usize = matches
        .value_of("cluster_num")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let limit: u64 = matches
        .value_of("limit")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let length: usize = matches
        .value_of("subchunk_len")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let local_config = match matches.value_of("read_type").unwrap() {
        "CCS" => ClusteringConfig::ccs(&dataset, cluster_num, length, limit),
        "CLR" => ClusteringConfig::clr(&dataset, cluster_num, length, limit),
        "ONT" => ClusteringConfig::ont(&dataset, cluster_num, length, limit),
        _ => unreachable!(),
    };
    let kmer: usize = matches
        .value_of("k")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let min_cluster_size = matches
        .value_of("min_cluster_size")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mat_score: i32 = matches
        .value_of("mat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let mismat_score: i32 = matches
        .value_of("mismat_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let gap_score: i32 = matches
        .value_of("gap_score")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let global_config = haplotyper::GlobalClusteringConfig::new(
        kmer,
        min_cluster_size,
        mat_score,
        mismat_score,
        gap_score,
    );
    let config = haplotyper::PolishClusteringConfig::new(mat_score, mismat_score, gap_score);
    let dataset = dataset
        .local_clustering(&local_config)
        .global_clustering(&global_config)
        .polish_clustering(&config);
    let config = AssembleConfig::default();
    let gfa = dataset.assemble_as_gfa(&config);
    let gfa_file = matches.value_of("gfa").unwrap();
    let mut gfa_path = std::path::PathBuf::from(&gfa_file);
    gfa_path.pop();
    match std::fs::create_dir_all(gfa_path) {
        Err(why) => error!("{:?}. Please check gfa path.", why),
        Ok(_) => match std::fs::File::create(gfa_file).map(BufWriter::new) {
            Ok(mut wtr) => writeln!(&mut wtr, "{}", gfa)?,
            _ => {}
        },
    }
    flush_file(&dataset)
}

fn get_input_file() -> std::io::Result<DataSet> {
    let stdin = std::io::stdin();
    let reader = BufReader::new(stdin.lock());
    match serde_json::de::from_reader(reader) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid Input from STDIN.");
            Err(std::io::Error::from(std::io::ErrorKind::Other))
        }
        Ok(res) => Ok(res),
    }
}

fn flush_file(dataset: &DataSet) -> std::io::Result<()> {
    let stdout = std::io::stdout();
    let mut wtr = std::io::BufWriter::new(stdout.lock());
    match serde_json::ser::to_writer_pretty(&mut wtr, dataset) {
        Err(why) => {
            eprintln!("{:?}", why);
            eprintln!("Invalid output to the STDOUT.");
            std::process::exit(1);
        }
        _ => Ok(()),
    }
}

fn main() -> std::io::Result<()> {
    let matches = App::new("jtk")
        .version("0.1")
        .author("Bansho Masutani")
        .about("HLA toolchain")
        .setting(clap::AppSettings::ArgRequiredElseHelp)
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .subcommand(subcommand_select_unit())
        .subcommand(subcommand_encode())
        .subcommand(subcommand_view())
        .subcommand(subcommand_local_clustering())
        .subcommand(subcommand_global_clustering())
        .subcommand(subcommand_polish_clustering())
        .subcommand(subcommand_assembly())
        .subcommand(subcommand_pipeline())
        .get_matches();
    match matches.subcommand() {
        ("entry", Some(sub_m)) => entry(sub_m),
        ("extract", Some(sub_m)) => extract(sub_m),
        ("stats", Some(sub_m)) => stats(sub_m),
        ("select_unit", Some(sub_m)) => select_unit(sub_m),
        ("encode", Some(sub_m)) => encode(sub_m),
        ("view", Some(sub_m)) => view(sub_m),
        ("local_clustering", Some(sub_m)) => local_clustering(sub_m),
        ("global_clustering", Some(sub_m)) => global_clustering(sub_m),
        ("polish_clustering", Some(sub_m)) => polish_clustering(sub_m),
        ("assemble", Some(sub_m)) => assembly(sub_m),
        ("pipeline", Some(sub_m)) => pipeline(sub_m),
        _ => Ok(()),
    }
}
