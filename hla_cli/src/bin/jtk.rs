use clap::{App, Arg, SubCommand};
use definitions::*;
use haplotyper::*;
use std::io::BufReader;
use std::io::{BufWriter, Write};
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
        .arg(
            Arg::with_name("input")
                .long("input")
                .short("r")
                .value_name("READS")
                .takes_value(true)
                .required(true)
                .help("Input FASTA file."),
        )
        .arg(
            Arg::with_name("read_type")
                .long("read_type")
                .takes_value(true)
                .default_value("CLR")
                .possible_values(&["CCS", "CLR", "ONT"])
                .help("Read type. CCS, CLR, or ONT."),
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
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .value_name("PATH")
                .required(true),
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
        .about("Pick subsequence from raw reads.")
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
                .default_value("2000")
                .help("Length of a chunk"),
        )
        .arg(
            Arg::with_name("skip_len")
                .short("s")
                .long("skip_len")
                .takes_value(true)
                .default_value("2000")
                .help("Margin between units"),
        )
        .arg(
            Arg::with_name("take_num")
                .short("n")
                .long("take_num")
                .takes_value(true)
                .default_value("3000")
                .help("Number of units;4*Genome size/chunk_len would be nice."),
        )
        .arg(
            Arg::with_name("margin")
                .short("m")
                .long("margin")
                .takes_value(true)
                .default_value("500")
                .help("Margin at the both end of a read."),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("number of threads"),
        )
        .arg(
            Arg::with_name("exclude")
                .long("exclude")
                .takes_value(true)
                .default_value("0.4")
                .help("filter out unit having more than [exclude] repetitiveness."),
        )
        .arg(
            Arg::with_name("upper")
                .short("u")
                .long("upper")
                .help("Discard units with occurence more than or equal to [upper].")
                .takes_value(true)
                .default_value("200"),
        )
        .arg(
            Arg::with_name("lower")
                .short("l")
                .long("lower")
                .help("Discard units with occurence less than or equal to [upper].")
                .takes_value(true)
                .default_value("4"),
        )
}

fn subcommand_polish_unit() -> App<'static, 'static> {
    SubCommand::with_name("polish_unit")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Polishing units by consuming encoded reads")
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
                .takes_value(true)
                .default_value("1")
                .help("number of threads"),
        )
        .arg(
            Arg::with_name("filter_size")
                .long("filter_size")
                .takes_value(true)
                .default_value("10")
                .help("Unit with coverage less than this value would be discarded."),
        )
        .arg(
            Arg::with_name("consensus_size")
                .long("consensus_size")
                .takes_value(true)
                .default_value("20")
                .help("The number of string to take consensus"),
        )
}

fn subcommand_mask_repeats() -> App<'static, 'static> {
    SubCommand::with_name("mask_repeats")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Mask Repeat(i.e., frequent k-mer)")
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
                .default_value("1"),
        )
        .arg(
            Arg::with_name("k")
                .short("k")
                .help("K-mer size(<32)")
                .takes_value(true)
                .default_value("15"),
        )
        .arg(
            Arg::with_name("freq")
                .short("f")
                .long("freq")
                .help("Mask top [freq] k-mer")
                .takes_value(true)
                .default_value("0.0008"),
        )
        .arg(
            Arg::with_name("min")
                .short("m")
                .long("min")
                .help("Prevent k-mer occuring less than [min] times from masking.")
                .takes_value(true)
                .default_value("10"),
        )
}

fn subcommand_encode() -> App<'static, 'static> {
    SubCommand::with_name("encode")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Encode reads by alignments (Internally invoke `minimap2` tools).")
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
                .default_value("1"),
        )
}

fn subcommand_polish_encoding() -> App<'static, 'static> {
    SubCommand::with_name("polish_encoding")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Remove nodes from reads.")
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
                .default_value("1"),
        )
}

fn subcommand_pick_components() -> App<'static, 'static> {
    SubCommand::with_name("pick_components")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Take top n largest components, discarding the rest and empty reads.")
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
                .default_value("1"),
        )
        .arg(
            Arg::with_name("component_num")
                .short("c")
                .long("component_num")
                .value_name("COMP")
                .help("Take top [COMP] largest connected-components.")
                .takes_value(true)
                .default_value("1"),
        )
}

fn subcommand_estimate_multiplicity() -> App<'static, 'static> {
    SubCommand::with_name("estimate_multiplicity")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Determine multiplicities of units.")
        .arg(Arg::with_name("verbose").short("v").multiple(true))
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("max_cluster_size")
                .short("m")
                .long("max_cluster_size")
                .required(false)
                .value_name("MAXCLUSTER")
                .help("Maximum number of cluster")
                .default_value("2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("seed")
                .short("s")
                .long("seed")
                .required(false)
                .value_name("SEED")
                .help("Seed for pseudorandon number generators.")
                .default_value("24")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("draft_assembly")
                .short("o")
                .long("draft_assembly")
                .required(false)
                .value_name("PATH")
                .help("If given, output draft GFA to PATH.")
                .takes_value(true),
        )
}

fn subcommand_partition_local() -> App<'static, 'static> {
    SubCommand::with_name("partition_local")
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("subchunk_len")
                .short("s")
                .long("subchunk_len")
                .required(false)
                .value_name("SubChunkLength")
                .help("The length of sub-chunks")
                .default_value("100")
                .takes_value(true),
        )
}

fn subcommand_correct_deletion() -> App<'static, 'static> {
    SubCommand::with_name("correct_deletion")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Correct unit deletion")
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("dissim")
                .short("d")
                .long("dis_similarity")
                .required(false)
                .value_name("DISSIM")
                .help("Acceptable Dissimilarity threshold.")
                .default_value("0.35")
                .takes_value(true),
        )
}

fn subcommand_partition_global() -> App<'static, 'static> {
    SubCommand::with_name("partition_global")
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
            Arg::with_name("graph")
                .long("graph")
                .help("Invoke graph-WhatsHap instead of de Bruijn."),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("k")
                .short("k")
                .long("kmer_size")
                .required(false)
                .value_name("KMER_SIZE")
                .help("The size of the kmer")
                .default_value("4")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_cluster_size")
                .short("m")
                .long("min_cluster_size")
                .required(false)
                .value_name("MIN_CLUSTER_SIZE")
                .help("The minimum size of a cluster")
                .default_value("50")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mat_score")
                .short("p")
                .long("match_score")
                .required(false)
                .value_name("MATCH_SCORE")
                .help("The match score")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mismat_score")
                .short("q")
                .long("mismatch_score")
                .required(false)
                .value_name("MISMATCH_SCORE")
                .help("The mismatch score")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gap_score")
                .short("g")
                .long("gap_score")
                .required(false)
                .value_name("GAP_SCORE")
                .help("The gap penalty")
                .default_value("2")
                .takes_value(true),
        )
}

fn subcommand_resolve_tangle() -> App<'static, 'static> {
    SubCommand::with_name("resolve_tangle")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Resolve tangle by re-estimate copy numbers.")
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("repeat_num")
                .short("r")
                .long("repeat_num")
                .required(false)
                .value_name("REPEAT_NUM")
                .help("Do EM algorithm for REPEAT_NUM times.")
                .default_value("7")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("coverage_threshold")
                .short("x")
                .long("threshold")
                .required(false)
                .value_name("THRESHOLD")
                .help("Unit with less that this coverage would be ignored.")
                .default_value("5")
                .takes_value(true),
        )
}

fn subcommand_correct_clustering() -> App<'static, 'static> {
    SubCommand::with_name("correct_clustering")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Correct local clustering by EM algorithm.")
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("repeat_num")
                .short("r")
                .long("repeat_num")
                .required(false)
                .value_name("REPEAT_NUM")
                .help("Do EM algorithm for REPEAT_NUM times.")
                .default_value("20")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("coverage_threshold")
                .short("x")
                .long("threshold")
                .required(false)
                .value_name("THRESHOLD")
                .help("Unit with less that this coverage would be ignored.")
                .default_value("5")
                .takes_value(true),
        )
}

fn subcommand_encode_densely() -> App<'static, 'static> {
    SubCommand::with_name("encode_densely")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Encoding homologoud diplotig in densely.")
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("length")
                .short("l")
                .long("length")
                .required(false)
                .value_name("LENGTH")
                .help("Contig shorter than this value would be compressed.")
                .default_value("15")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_span_reads")
                .short("m")
                .long("min_span_reads")
                .required(false)
                .value_name("MIN")
                .help("Minimum required number of reads to solve repeats")
                .default_value("4")
                .takes_value(true),
        )
}

fn subcommand_assemble() -> App<'static, 'static> {
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
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("window_size")
                .short("w")
                .long("window_size")
                .required(false)
                .value_name("WINDOW_SIZE")
                .help("Size of the window to take consensus sequences.")
                .default_value("2000")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("min_span_reads")
                .short("m")
                .long("min_span_reads")
                .required(false)
                .value_name("MIN")
                .help("Minimum required number of reads to solve repeats")
                .default_value("4")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("no_polish")
                .short("n")
                .long("no_polish")
                .help("If this flag is given, polishing stage would be skipped."),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .required(true)
                .value_name("PATH")
                .help("Output file name")
                .takes_value(true),
        )
}

fn entry(matches: &clap::ArgMatches) -> std::io::Result<DataSet> {
    debug!("START\tEntry");
    let file = matches.value_of("input").unwrap();
    let reader = std::fs::File::open(file).map(BufReader::new)?;
    let seqs = bio_utils::fasta::parse_into_vec_from(reader)?;
    debug!("Encoding {} reads", seqs.len());
    let read_type = matches.value_of("read_type").unwrap();
    Ok(DataSet::new(file, &seqs, read_type))
}

fn extract(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tExtract");
    debug!("Target is {}", matches.value_of("target").unwrap());
    let file = std::fs::File::create(matches.value_of("output").unwrap())?;
    match matches.value_of("target").unwrap() {
        "raw_reads" => {
            let mut wtr = bio_utils::fasta::Writer::new(file);
            for seq in dataset.extract_fasta(ExtractTarget::RawReads) {
                wtr.write_record(&seq)?;
            }
        }
        "hic_reads" => {
            let mut wtr = bio_utils::fasta::Writer::new(file);
            for seq in dataset.extract_fasta(ExtractTarget::HiCReads) {
                wtr.write_record(&seq)?;
            }
        }
        "units" => {
            let mut wtr = bio_utils::fasta::Writer::new(file);
            for seq in dataset.extract_fasta(ExtractTarget::Units) {
                wtr.write_record(&seq)?;
            }
        }
        "assignments" => {
            let asn_name_desc = dataset.extract_assignments();
            let mut wtr = BufWriter::new(file);
            for (asn, name, desc) in asn_name_desc {
                writeln!(&mut wtr, "{}\t{}\t{}", asn, name, desc)?;
            }
        }
        &_ => unreachable!(),
    };
    Ok(dataset)
}

fn stats(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tStats step");
    let wtr = std::io::BufWriter::new(std::fs::File::create(matches.value_of("file").unwrap())?);
    dataset.stats(wtr)?;
    Ok(dataset)
}

fn select_unit(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tSelecting Units");
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
    let take_num: usize = matches
        .value_of("take_num")
        .and_then(|e| e.parse().ok())
        .expect("Take num");
    let thrds: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse().ok())
        .expect("threads");
    let filter: f64 = matches
        .value_of("exclude")
        .and_then(|e| e.parse().ok())
        .expect("exclude");
    let upper: usize = matches
        .value_of("upper")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    let lower: usize = matches
        .value_of("lower")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(thrds)
        .build_global()
    {
        debug!("{:?}:If you run `pipeline` module, this is Harmless.", why);
    }
    use ReadType::*;
    let (cl, tn) = (chunk_len, take_num);
    let config = match dataset.read_type {
        CCS => UnitConfig::new_ccs(cl, tn, skip_len, margin, thrds, filter, upper, lower),
        CLR => UnitConfig::new_clr(cl, tn, skip_len, margin, thrds, filter, upper, lower),
        _ => UnitConfig::new_ont(cl, tn, skip_len, margin, thrds, filter, upper, lower),
    };
    Ok(dataset.select_chunks(&config))
}

fn repeat_masking(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tmasking repeat.");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    let k: usize = matches.value_of("k").and_then(|l| l.parse().ok()).unwrap();
    if k > 32 {
        panic!("K should be less than 32.");
    }
    let freq: f64 = matches
        .value_of("freq")
        .and_then(|l| l.parse().ok())
        .unwrap();
    let min: u32 = matches
        .value_of("min")
        .and_then(|l| l.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let config = haplotyper::RepeatMaskConfig::new(k, freq, min);
    Ok(dataset.mask_repeat(&config))
}

fn encode(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tEncoding step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    // TODO:Branching by read type.
    Ok(dataset.encode(threads, haplotyper::encode::CLR_CLR_SIM))
}

fn polish_encode(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tPolish encoding.");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    use haplotyper::remove_erroneous_nodes::RemoveErroneousNodes;
    Ok(dataset.remove_erroneous_nodes())
}

fn pick_components(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tpicking components.");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let component_num: usize = matches
        .value_of("component_num")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let config = ComponentPickingConfig::new(component_num);
    Ok(dataset.pick_top_n_component(&config))
}

fn polish_unit(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tpolishing units");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    let filter_size: usize = matches
        .value_of("filter_size")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    let consensus_size: usize = matches
        .value_of("consensus_size")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let config = PolishUnitConfig::new(dataset.read_type, filter_size, consensus_size);
    Ok(dataset.polish_unit(&config))
}
fn multiplicity_estimation(
    matches: &clap::ArgMatches,
    dataset: DataSet,
) -> std::io::Result<DataSet> {
    debug!("START\tmultiplicity estimation");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let max_cluster_size: usize = matches
        .value_of("max_cluster_size")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let seed: u64 = matches
        .value_of("seed")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let path = matches.value_of("draft_assembly");
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let config = MultiplicityEstimationConfig::new(threads, max_cluster_size, seed, path);
    Ok(dataset.estimate_multiplicity_graph(&config))
}

fn local_clustering(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tLocal Clustering step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let length: usize = matches
        .value_of("subchunk_len")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let config = ClusteringConfig::with_default(&dataset, 2, length);
    Ok(dataset.local_clustering(&config))
}

fn correct_deletion(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tCorrectDeletion");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let sim_thr: f64 = matches
        .value_of("dissim")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    use haplotyper::encode::deletion_fill;
    Ok(deletion_fill::correct_unit_deletion(dataset, sim_thr))
}

fn global_clustering(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tGlobal Clustering step");
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
    // Minus.
    let mismat_score: i32 = -matches
        .value_of("mismat_score")
        .and_then(|num| num.parse::<i32>().ok())
        .unwrap();
    let gap_score: i32 = -matches
        .value_of("gap_score")
        .and_then(|num| num.parse::<i32>().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let config = haplotyper::GlobalClusteringConfig::new(
        kmer,
        min_cluster_size,
        mat_score,
        mismat_score,
        gap_score,
    );
    if matches.is_present("graph") {
        Ok(dataset.global_clustering_graph(&config))
    } else {
        Ok(dataset.global_clustering(&config))
    }
}

fn clustering_correction(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tClustering Correction step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let repeat_num: usize = matches
        .value_of("repeat_num")
        .and_then(|num| num.parse::<usize>().ok())
        .unwrap();
    let threshold: usize = matches
        .value_of("coverage_threshold")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    Ok(dataset.correct_clustering_em(repeat_num, threshold, true))
}

fn resolve_tangle(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tClustering Correction step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let repeat_num: usize = matches
        .value_of("repeat_num")
        .and_then(|num| num.parse::<usize>().ok())
        .unwrap();
    let threshold: usize = matches
        .value_of("coverage_threshold")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    use haplotyper::re_clustering::*;
    let config = ReClusteringConfig::new(threads, repeat_num, threshold);
    Ok(dataset.re_clustering(&config))
}

fn encode_densely(matches: &clap::ArgMatches, dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tEncode densely");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let length: usize = matches
        .value_of("length")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let min_span_reads: usize = matches
        .value_of("min_span_reads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    use haplotyper::dense_encoding::*;
    let config = DenseEncodingConfig::new(length, min_span_reads);
    Ok(dataset.dense_encoding(&config))
}

fn assembly(matches: &clap::ArgMatches, mut dataset: DataSet) -> std::io::Result<DataSet> {
    debug!("START\tAssembly step");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let window_size: usize = matches
        .value_of("window_size")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let min_span_reads: usize = matches
        .value_of("min_span_reads")
        .and_then(|num| num.parse().ok())
        .unwrap();

    if let Err(why) = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
    {
        debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
    }
    let skip_polish = matches.is_present("no_polish");
    let file = matches.value_of("output").unwrap();
    let mut file = std::fs::File::create(file).map(BufWriter::new)?;
    let config = AssembleConfig::new(threads, window_size, !skip_polish, true, min_span_reads);
    if dataset.assignments.is_empty() {
        dataset.assignments = dataset
            .encoded_reads
            .iter()
            .map(|r| Assignment::new(r.id, 0))
            .collect();
    }
    // TODO:Parametrize here.
    let min_count = 6;
    let max_tig_length = 20;
    dataset.zip_up_suspicious_haplotig(&config, min_count, max_tig_length);
    // let squish_limit = 1;
    // dataset.squish_small_contig(&config, squish_limit);
    debug!("START\tFinal assembly");
    let gfa = dataset.assemble(&config);
    writeln!(&mut file, "{}", gfa)?;
    Ok(dataset)
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
    match serde_json::ser::to_writer(&mut wtr, dataset) {
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
        .subcommand(subcommand_polish_unit())
        .subcommand(subcommand_encode())
        .subcommand(subcommand_polish_encoding())
        .subcommand(subcommand_estimate_multiplicity())
        .subcommand(subcommand_resolve_tangle())
        .subcommand(subcommand_partition_local())
        .subcommand(subcommand_correct_deletion())
        .subcommand(subcommand_partition_global())
        .subcommand(subcommand_correct_clustering())
        .subcommand(subcommand_encode_densely())
        .subcommand(subcommand_assemble())
        .subcommand(subcommand_pick_components())
        .subcommand(subcommand_mask_repeats())
        .get_matches();
    if let Some(sub_m) = matches.subcommand().1 {
        let level = match sub_m.occurrences_of("verbose") {
            0 => "warn",
            1 => "info",
            2 => "debug",
            _ => "trace",
        };
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    }
    if let ("entry", Some(sub_m)) = matches.subcommand() {
        return entry(sub_m).and_then(|x| flush_file(&x));
    }
    let ds = get_input_file()?;
    let result = match matches.subcommand() {
        ("select_unit", Some(sub_m)) => select_unit(sub_m, ds),
        ("mask_repeats", Some(sub_m)) => repeat_masking(sub_m, ds),
        ("encode", Some(sub_m)) => encode(sub_m, ds),
        ("polish_unit", Some(sub_m)) => polish_unit(sub_m, ds),
        ("pick_components", Some(sub_m)) => pick_components(sub_m, ds),
        ("polish_encoding", Some(sub_m)) => polish_encode(sub_m, ds),
        ("partition_local", Some(sub_m)) => local_clustering(sub_m, ds),
        ("correct_deletion", Some(sub_m)) => correct_deletion(sub_m, ds),
        ("resolve_tangle", Some(sub_m)) => resolve_tangle(sub_m, ds),
        ("estimate_multiplicity", Some(sub_m)) => multiplicity_estimation(sub_m, ds),
        ("partition_global", Some(sub_m)) => global_clustering(sub_m, ds),
        ("correct_clustering", Some(sub_m)) => clustering_correction(sub_m, ds),
        ("encode_densely", Some(sub_m)) => encode_densely(sub_m, ds),
        ("assemble", Some(sub_m)) => assembly(sub_m, ds),
        ("extract", Some(sub_m)) => extract(sub_m, ds),
        ("stats", Some(sub_m)) => stats(sub_m, ds),
        _ => unreachable!(),
    };
    result.and_then(|x| flush_file(&x))
}
