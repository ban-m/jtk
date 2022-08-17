use clap::{Arg, Command};
use definitions::*;
use std::io::BufReader;
use std::io::{BufWriter, Write};
#[macro_use]
extern crate log;
fn subcommand_entry() -> Command<'static> {
    Command::new("entry")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Entry point. It encodes a fasta file into HLA-class file.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("input")
                .long("input")
                .short('r')
                .value_name("READS")
                .takes_value(true)
                .required(true)
                .help("Input FASTA/Q file."),
        )
        .arg(
            Arg::new("read_type")
                .long("read_type")
                .takes_value(true)
                .default_value("CLR")
                .possible_values(&["CCS", "CLR", "ONT"])
                .help("Read type. CCS, CLR, or ONT."),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("number of threads"),
        )
}

fn subcommand_extract() -> Command<'static> {
    Command::new("extract")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Extract (unfold) all the information in the packed file into one tsv file")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .takes_value(true)
                .value_name("PATH")
                .required(true),
        )
}

fn subcommand_stats() -> Command<'static> {
    Command::new("stats")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Write stats to the specified file.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("file")
                .long("file")
                .value_name("FILE")
                .short('f')
                .required(true)
                .takes_value(true),
        )
}

fn subcommand_select_unit() -> Command<'static> {
    Command::new("select_unit")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Pick subsequence from raw reads.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("chunk_len")
                .short('l')
                .long("chunk_len")
                .takes_value(true)
                .default_value("2000")
                .help("Length of a chunk"),
        )
        .arg(
            Arg::new("take_num")
                .short('n')
                .long("take_num")
                .takes_value(true)
                .default_value("500")
                .help("Number of units; Genome size/chunk_len would be nice."),
        )
        .arg(
            Arg::new("margin")
                .short('m')
                .long("margin")
                .takes_value(true)
                .default_value("500")
                .help("Margin at the both end of a read."),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("number of threads"),
        )
        .arg(
            Arg::new("exclude")
                .long("exclude")
                .takes_value(true)
                .default_value("0.6")
                .help("filter out unit having more than [exclude] repetitiveness."),
        )
        .arg(
            Arg::new("upper")
                .short('u')
                .long("upper")
                .help("Discard units with occurence more than or equal to [upper].")
                .takes_value(true)
                .default_value("150"),
        )
        .arg(
            Arg::new("lower")
                .short('l')
                .long("lower")
                .help("Discard units with occurence less than or equal to [upper].")
                .takes_value(true)
                .default_value("4"),
        )
}

fn subcommand_mask_repeats() -> Command<'static> {
    Command::new("mask_repeats")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Mask Repeat(i.e., frequent k-mer)")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .help("Number of threads")
                .takes_value(true)
                .default_value("1"),
        )
        .arg(
            Arg::new("k")
                .short('k')
                .help("K-mer size(<32)")
                .takes_value(true)
                .default_value("20"),
        )
        .arg(
            Arg::new("freq")
                .short('f')
                .long("freq")
                .help("Mask top [freq] k-mer")
                .takes_value(true)
                .default_value("0.0005"),
        )
        .arg(
            Arg::new("min")
                .short('m')
                .long("min")
                .help("Prevent k-mer occuring less than [min] times from masking.")
                .takes_value(true)
                .default_value("10"),
        )
}

fn subcommand_encode() -> Command<'static> {
    Command::new("encode")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Encode reads by alignments (Internally invoke `minimap2` tools).")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .help("Number of threads")
                .takes_value(true)
                .default_value("1"),
        )
}

fn subcommand_polish_encoding() -> Command<'static> {
    Command::new("polish_encoding")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Remove nodes from reads.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .help("Number of threads")
                .takes_value(true)
                .default_value("1"),
        )
}

fn subcommand_pick_components() -> Command<'static> {
    Command::new("pick_components")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Take top n largest components, discarding the rest and empty reads.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .short('t')
                .help("Number of threads")
                .takes_value(true)
                .default_value("1"),
        )
        .arg(
            Arg::new("component_num")
                .short('c')
                .long("component_num")
                .value_name("COMP")
                .help("Take top [COMP] largest connected-components.")
                .takes_value(true)
                .default_value("1"),
        )
}

fn subcommand_estimate_multiplicity() -> Command<'static> {
    Command::new("estimate_multiplicity")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Determine multiplicities of units.")
        .arg(Arg::new("verbose").short('v').multiple_occurrences(true))
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("seed")
                .short('s')
                .long("seed")
                .required(false)
                .value_name("SEED")
                .help("Seed for pseudorandon number generators.")
                .default_value("24")
                .takes_value(true),
        )
        .arg(
            Arg::new("draft_assembly")
                .short('o')
                .long("draft_assembly")
                .required(false)
                .value_name("PATH")
                .help("If given, output draft GFA to PATH.")
                .takes_value(true),
        )
        .arg(
            Arg::new("purge")
                .short('p')
                .long("purge_copy_num")
                .required(false)
                .value_name("COPY NUM")
                .help("If given, remove all the chunks with copy number >= COPY_NUM")
                .takes_value(true),
        )
        .arg(
            Arg::new("coverage")
                .short('c')
                .long("coverage")
                .required(false)
                .value_name("HAP COV")
                .help("If given, use this value as a haploid coverage estimate")
                .takes_value(true),
        )
}

fn subcommand_partition_local() -> Command<'static> {
    Command::new("partition_local")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Clustering reads. (Local)")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
}

fn subcommand_purge_diverged() -> Command<'static> {
    Command::new("purge_diverged")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Purge diverged clusters")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
}

fn subcommand_correct_deletion() -> Command<'static> {
    Command::new("correct_deletion")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Correct unit deletion")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("re_cluster")
                .short('r')
                .long("re_cluster")
                .required(false)
                .help("Re-calculate the posterior probability of newly encoded units."),
        )
}

fn subcommand_correct_clustering() -> Command<'static> {
    Command::new("correct_clustering")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Correct local clustering by EM algorithm.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("repeat_num")
                .short('r')
                .long("repeat_num")
                .required(false)
                .value_name("REPEAT_NUM")
                .help("Do EM algorithm for REPEAT_NUM times.")
                .default_value("5")
                .takes_value(true),
        )
        .arg(
            Arg::new("coverage_threshold")
                .short('x')
                .long("threshold")
                .required(false)
                .value_name("THRESHOLD")
                .help("Unit with less that this coverage would be ignored.")
                .default_value("5")
                .takes_value(true),
        )
}

fn subcommand_encode_densely() -> Command<'static> {
    Command::new("encode_densely")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Encoding homologoud diplotig in densely.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("length")
                .short('l')
                .long("length")
                .required(false)
                .value_name("LENGTH")
                .help("Contig shorter than this value would be compressed.")
                .default_value("15")
                .takes_value(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .required(false)
                .value_name("PATH")
                .help("Dump the intermediate assembly to the PATH")
                .takes_value(true),
        )
}

fn subcommand_assemble() -> Command<'static> {
    Command::new("assemble")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Assemble reads.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("window_size")
                .short('w')
                .long("window_size")
                .required(false)
                .value_name("WINDOW_SIZE")
                .help("Size of the window to take consensus sequences.")
                .default_value("2000")
                .takes_value(true),
        )
        .arg(
            Arg::new("no_polish")
                .short('n')
                .long("no_polish")
                .help("If this flag is given, polishing stage would be skipped."),
        )
        .arg(
            Arg::new("min_llr")
                .long("min_llr")
                .takes_value(true)
                .value_name("LK Ratio")
                .required(false)
                .help("Minimum likelihood ratio"),
        )
        .arg(
            Arg::new("min_span")
                .long("min_span")
                .takes_value(true)
                .value_name("# of reads")
                .required(false)
                .help("Minimum required reads to span repeats"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .required(true)
                .value_name("PATH")
                .help("Output file name")
                .takes_value(true),
        )
}

fn subcommand_polish() -> Command<'static> {
    Command::new("polish")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Polish contigs.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .required(false)
                .value_name("THREADS")
                .help("Number of Threads")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::new("window_size")
                .short('w')
                .long("window_size")
                .required(false)
                .value_name("WINDOW_SIZE")
                .help("Size of the window to take consensus sequences.")
                .default_value("2000")
                .takes_value(true),
        )
        .arg(
            Arg::new("reads")
                .short('r')
                .long("reads")
                .takes_value(true)
                .required(true)
                .value_name("Reads<FASTA|FASTQ>")
                .help("raw reads. The extension is used to decide format."),
        )
        .arg(
            Arg::new("contigs")
                .short('c')
                .long("contigs")
                .required(true)
                .takes_value(true)
                .value_name("Contig<GFA|FA>")
                .help("Contigs to be polished. The extension is used to decide format"),
        )
        .arg(
            Arg::new("alignments")
                .short('a')
                .long("alignments")
                .takes_value(true)
                .required(true)
                .value_name("Alignments<PAF|SAM>")
                .help("Alignments reads -> contigs. If -, read from the stdin."),
        )
        .arg(
            Arg::new("format")
                .short('f')
                .long("format")
                .takes_value(true)
                .required(true)
                .possible_values(&["sam", "paf"])
                .help("Format of the alignments")
                .default_value("sam"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .takes_value(true)
                .default_value("42"),
        )
}

fn subcommand_pipeline() -> Command<'static> {
    Command::new("pipeline")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Run pipeline based on the given TOML file.")
        .arg(
            Arg::new("profile")
                .short('p')
                .takes_value(true)
                .required(true)
                .help("TOML configuration file. See example.toml for an example."),
        )
}

fn entry(matches: &clap::ArgMatches) -> std::io::Result<DataSet> {
    use haplotyper::entry::Entry;
    debug!("START\tEntry");
    set_threads(matches);
    let file = matches.value_of("input").unwrap();
    let reader = std::fs::File::open(file).map(BufReader::new)?;
    debug!("Opening {}", file);
    let seqs: Vec<(String, Vec<u8>)> = match file.chars().last() {
        Some('a') => bio_utils::fasta::parse_into_vec_from(reader)?
            .into_iter()
            .map(|records| {
                let (id, _, seq) = records.into();
                let seq = seq.into_bytes();
                (id, seq)
            })
            .collect(),
        Some('q') => bio_utils::fastq::parse_into_vec_from(reader)?
            .into_iter()
            .map(|record| {
                if record.seq().iter().any(|x| !b"ACGT".contains(x)) {
                    debug!("{}", record);
                }
                let (id, seq, _) = record.into();
                (id, seq)
            })
            .collect(),
        _ => panic!("file type:{} not supported", file),
    };
    let read_type = matches.value_of("read_type").unwrap();
    Ok(DataSet::entry(file, seqs, read_type))
}

fn extract(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    use haplotyper::extract::Extract;
    debug!("START\tExtract");
    debug!("Target is {}", matches.value_of("target").unwrap());
    let file = std::fs::File::create(matches.value_of("output").unwrap())?;
    let mut wtr = std::io::BufWriter::new(file);
    dataset.extract(&mut wtr)?;
    Ok(())
}

fn stats(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    use haplotyper::stats::Stats;
    debug!("START\tStats step");
    let wtr = std::io::BufWriter::new(std::fs::File::create(matches.value_of("file").unwrap())?);
    dataset.sanity_check();
    dataset.stats(wtr)?;
    Ok(())
}

fn select_unit(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tSelecting Units");
    let chunk_len: usize = matches
        .value_of("chunk_len")
        .and_then(|e| e.parse().ok())
        .expect("Chunk len");
    let margin: usize = matches
        .value_of("margin")
        .and_then(|e| e.parse().ok())
        .expect("Margin");
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
    set_threads(matches);
    use haplotyper::determine_units::{DetermineUnit, DetermineUnitConfig};
    let (cl, tn) = (chunk_len, take_num);
    let config = DetermineUnitConfig::new(cl, tn, margin, thrds, filter, upper, lower);
    dataset.select_chunks(&config);
}

fn repeat_masking(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tmasking repeat.");
    set_threads(matches);
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
    use haplotyper::repeat_masking::*;
    let config = RepeatMaskConfig::new(k, freq, min);
    dataset.mask_repeat(&config)
}

fn encode(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tEncoding step");
    set_threads(matches);
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    use haplotyper::encode::Encode;
    dataset.encode(threads, dataset.read_type.sim_thr())
}

fn polish_encode(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tPolish encoding.");
    set_threads(matches);
    use haplotyper::remove_erroneous_nodes::RemoveErroneousNodes;
    dataset.remove_erroneous_nodes();
}

fn pick_components(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tpicking components.");
    set_threads(matches);
    let component_num: usize = matches
        .value_of("component_num")
        .and_then(|num| num.parse().ok())
        .unwrap();
    use haplotyper::pick_component::*;
    let config = ComponentPickingConfig::new(component_num);
    dataset.pick_top_n_component(&config);
}
fn multiplicity_estimation(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tmultiplicity estimation");
    let seed: u64 = matches
        .value_of("seed")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let cov: Option<f64> = matches.value_of("coverage").and_then(|e| e.parse().ok());
    dataset.coverage = match cov {
        Some(cov) => definitions::Coverage::Protected(cov),
        None => definitions::Coverage::NotAvailable,
    };
    let path = matches.value_of("draft_assembly");
    set_threads(matches);
    use haplotyper::multiplicity_estimation::*;
    let config = MultiplicityEstimationConfig::new(seed, path);
    dataset.estimate_multiplicity(&config);
    let purge: Option<usize> = matches.value_of("purge").and_then(|x| x.parse().ok());
    if let Some(upper_copy_num) = purge {
        dataset.purge_multiplicity(upper_copy_num);
    }
}

fn local_clustering(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tLocal Clustering step");
    set_threads(matches);
    use haplotyper::local_clustering::*;
    dataset.local_clustering();
}

fn purge_diverged(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tPurge diverged clusters");
    set_threads(matches);
    use haplotyper::purge_diverged::*;
    let config = PurgeDivConfig::new();
    dataset.purge(&config);
}

fn correct_deletion(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tCorrectDeletion");
    set_threads(matches);
    let to_recal = matches.is_present("re_cluster");
    use haplotyper::determine_units::calc_sim_thr;
    use haplotyper::determine_units::TAKE_THR;
    let sim_thr = calc_sim_thr(dataset, TAKE_THR);
    use haplotyper::encode::deletion_fill::*;
    let config = CorrectDeletionConfig::new(to_recal, sim_thr);
    debug!("SIMTHR\t{sim_thr}");
    dataset.correct_deletion(&config);
}

fn clustering_correction(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tClustering Correction step");
    let _repeat_num: usize = matches
        .value_of("repeat_num")
        .and_then(|num| num.parse::<usize>().ok())
        .unwrap();
    let _threshold: usize = matches
        .value_of("coverage_threshold")
        .and_then(|num| num.parse().ok())
        .unwrap();
    set_threads(matches);
    use haplotyper::phmm_likelihood_correction::*;
    let config = CorrectionConfig::default();
    dataset.correct_clustering(&config);
}

fn encode_densely(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tEncode densely");
    set_threads(matches);
    let length: usize = matches
        .value_of("length")
        .and_then(|num| num.parse().ok())
        .unwrap();
    use haplotyper::dense_encoding::*;
    let file = matches.value_of("output");
    let config = DenseEncodingConfig::new(length, file);
    dataset.dense_encoding_dev(&config);
}

fn assembly(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    debug!("START\tAssembly step");
    set_threads(matches);
    let window_size: usize = matches
        .value_of("window_size")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let min_llr: Option<f64> = matches.value_of("min_llr").map(|num| num.parse().unwrap());
    let min_span: Option<usize> = matches.value_of("min_span").map(|num| num.parse().unwrap());
    let skip_polish = matches.is_present("no_polish");
    let file = matches.value_of("output").unwrap();
    let mut file = std::fs::File::create(file).map(BufWriter::new)?;
    use haplotyper::assemble::*;
    let msr = min_span.unwrap_or_else(|| dataset.read_type.min_span_reads());
    let min_lk = min_llr.unwrap_or_else(|| dataset.read_type.min_llr_value());
    let config = AssembleConfig::new(window_size, !skip_polish, true, msr, min_lk, true);
    debug!("START\tFinal assembly");
    if !skip_polish {
        use haplotyper::model_tune::update_model;
        update_model(dataset);
    }
    let gfa = dataset.assemble(&config);
    writeln!(file, "{}", gfa)?;
    Ok(())
}

fn polish(matches: &clap::ArgMatches) -> std::io::Result<()> {
    set_threads(matches);
    let window_size: usize = matches
        .value_of("window_size")
        .and_then(|x| x.parse().ok())
        .unwrap();
    let reads = matches.value_of("reads").unwrap();
    let contig = matches.value_of("contigs").unwrap();
    let alignments = matches.value_of("alignments").unwrap();
    let format = matches.value_of("format").unwrap();
    let seed: u64 = matches
        .value_of("seed")
        .and_then(|x| x.parse().ok())
        .unwrap();
    use haplotyper::polish_segments::polish_segmnents;
    polish_segmnents(reads, contig, alignments, format, window_size, seed)
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
    let matches = Command::new("jtk")
        .version("0.1")
        .author("Bansho Masutani <ban-m@g.ecc.u-tokyo.ac.jp>")
        .about("JTK regional diploid assembler")
        .arg_required_else_help(true)
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .subcommand(subcommand_select_unit())
        .subcommand(subcommand_encode())
        .subcommand(subcommand_polish_encoding())
        .subcommand(subcommand_estimate_multiplicity())
        .subcommand(subcommand_partition_local())
        .subcommand(subcommand_purge_diverged())
        .subcommand(subcommand_correct_deletion())
        .subcommand(subcommand_correct_clustering())
        .subcommand(subcommand_encode_densely())
        .subcommand(subcommand_assemble())
        .subcommand(subcommand_pick_components())
        .subcommand(subcommand_mask_repeats())
        .subcommand(subcommand_polish())
        .subcommand(subcommand_pipeline())
        .get_matches();
    if let Some(("pipeline", sub_m)) = matches.subcommand() {
        let path = sub_m.value_of("profile").unwrap();
        use std::io::Read;
        let mut rdr = std::fs::File::open(path).map(std::io::BufReader::new)?;
        let mut file = String::new();
        rdr.read_to_string(&mut file)?;
        let config: hla_cli::PipelineConfig = toml::from_str(&file).unwrap();
        return hla_cli::run_pipeline(&config);
    }
    if let Some((_, sub_m)) = matches.subcommand() {
        let level = match sub_m.occurrences_of("verbose") {
            0 => "warn",
            1 => "info",
            2 => "debug",
            _ => "trace",
        };
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    }
    if let Some(("entry", sub_m)) = matches.subcommand() {
        return entry(sub_m).and_then(|x| flush_file(&x));
    }
    if let Some(("polish", sub_m)) = matches.subcommand() {
        return polish(sub_m);
    }
    let mut ds = get_input_file()?;
    let ds = &mut ds;
    match matches.subcommand() {
        Some(("select_unit", sub_m)) => select_unit(sub_m, ds),
        Some(("mask_repeats", sub_m)) => repeat_masking(sub_m, ds),
        Some(("encode", sub_m)) => encode(sub_m, ds),
        Some(("pick_components", sub_m)) => pick_components(sub_m, ds),
        Some(("polish_encoding", sub_m)) => polish_encode(sub_m, ds),
        Some(("partition_local", sub_m)) => local_clustering(sub_m, ds),
        Some(("purge_diverged", sub_m)) => purge_diverged(sub_m, ds),
        Some(("correct_deletion", sub_m)) => correct_deletion(sub_m, ds),
        Some(("estimate_multiplicity", sub_m)) => multiplicity_estimation(sub_m, ds),
        Some(("correct_clustering", sub_m)) => clustering_correction(sub_m, ds),
        Some(("encode_densely", sub_m)) => encode_densely(sub_m, ds),
        Some(("assemble", sub_m)) => assembly(sub_m, ds).unwrap(),
        Some(("extract", sub_m)) => extract(sub_m, ds).unwrap(),
        Some(("stats", sub_m)) => stats(sub_m, ds).unwrap(),
        _ => unreachable!(),
    };
    flush_file(ds)
}

fn set_threads(matches: &clap::ArgMatches) {
    if let Some(threads) = matches.value_of("threads").and_then(|num| num.parse().ok()) {
        debug!("Set Threads\t{}", threads);
        if let Err(why) = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
        {
            debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
        }
    }
}
