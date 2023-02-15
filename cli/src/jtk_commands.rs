use clap::{Arg, Command};
fn subcommand_entry() -> Command<'static> {
    Command::new("entry")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Entry point. It encodes a fasta file into JSON file.")
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
                .possible_values(["CCS", "CLR", "ONT"])
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
        .about("Extract all the information in the packed file into one tsv")
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

fn subcommand_select_chunk() -> Command<'static> {
    Command::new("select_chunks")
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
                .help("Number of chunks; Genome size/chunk_len would be nice."),
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
                .default_value("0.8")
                .help("filter out chunks having more than [exclude] repetitiveness."),
        )
        .arg(
            Arg::new("purge_copy_num")
                .short('u')
                .long("purge_copy_num")
                .help("Discard chunks with copy_number more than or equal to [upper].")
                .takes_value(true)
                .default_value("10"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .help("Seed value for random number generators")
                .takes_value(true)
                .default_value("42"),
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
                .default_value("12"),
        )
        .arg(
            Arg::new("freq")
                .short('f')
                .long("freq")
                .help("Mask top [freq] non-singleton k-mer")
                .takes_value(true)
                .default_value("0.001"),
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
        .arg(
            Arg::new("sim_thr")
                .long("sim_thr")
                .help("similarity threshold")
                .takes_value(true),
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
        .about("Determine multiplicities of chunks.")
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
        .about("Correct deletions of chunks inside the reads.")
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
                .help("Re-calculate the posterior probability of newly encoded chunks."),
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

fn subcommand_squish() -> Command<'static> {
    Command::new("squish")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Squish erroneous clusters")
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
            Arg::new("ari")
                .short('a')
                .long("ari")
                .required(false)
                .value_name("ARI")
                .help("clusters with ARI smaller than [ARI] would be compressed")
                .default_value("0.4")
                .takes_value(true),
        )
        .arg(
            Arg::new("match_score")
                .long("match_score")
                .required(false)
                .help("Match score for clusters with adj.rand index above [ARI]")
                .default_value("4.0")
                .takes_value(true),
        )
        .arg(
            Arg::new("mismatch_score")
                .long("mismatch_score")
                .required(false)
                .help("Mismatch score for clusters with adj.rand index below [ARI]")
                .default_value("-1.0")
                .takes_value(true),
        )
        .arg(
            Arg::new("count")
                .short('c')
                .long("count")
                .required(false)
                .value_name("COUNT")
                .help("Min req coverage to compute ARI.")
                .default_value("7")
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
                .default_value("1")
                .required(false)
                .help("Minimum likelihood ratio"),
        )
        .arg(
            Arg::new("min_span")
                .long("min_span")
                .takes_value(true)
                .value_name("# of reads")
                .default_value("2")
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
                .possible_values(["sam", "paf"])
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
pub fn jtk_parser() -> clap::Command<'static> {
    clap::Command::new("jtk")
        .version("0.1")
        .author("Bansho Masutani <ban-m@g.ecc.u-tokyo.ac.jp>")
        .about("JTK regional diploid assembler")
        .arg_required_else_help(true)
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .subcommand(subcommand_select_chunk())
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
        .subcommand(subcommand_squish())
        .subcommand(subcommand_polish())
        .subcommand(subcommand_pipeline())
}
