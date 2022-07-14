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
        // .arg(
        //     Arg::new("slag")
        //         .long("slag")
        //         .takes_value(true)
        //         .value_name("PATH")
        //         .help("Dump low-quality reads into PATH"),
        // )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .takes_value(true)
                .default_value("1")
                .help("number of threads"),
        )
}

const TARGETS: [&str; 4] = ["raw_reads", "hic_reads", "units", "assignments"];
fn subcommand_extract() -> Command<'static> {
    Command::new("extract")
        .version("0.1")
        .author("Bansho Masutani")
        .about("Exit point. It extract fasta/q file from a HLA-class file.")
        .arg(
            Arg::new("verbose")
                .short('v')
                .multiple_occurrences(true)
                .help("Debug mode"),
        )
        .arg(
            Arg::new("target")
                .short('t')
                .long("target")
                .takes_value(true)
                .value_name("TARGET")
                .required(true)
                .possible_values(&TARGETS),
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
            Arg::new("skip_len")
                .short('s')
                .long("skip_len")
                .takes_value(true)
                .default_value("2000")
                .help("Margin between units"),
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

// fn subcommand_polish_unit() -> Command<'static> {
//     Command::new("polish_unit")
//         .version("0.1")
//         .author("BanshoMasutani")
//         .about("Polishing units by consuming encoded reads")
//         .arg(
//             Arg::new("verbose")
//                 .short('v')
//                 .multiple_occurrences(true)
//                 .help("Debug mode"),
//         )
//         .arg(
//             Arg::new("threads")
//                 .short('t')
//                 .long("threads")
//                 .takes_value(true)
//                 .default_value("1")
//                 .help("number of threads"),
//         )
//         .arg(
//             Arg::new("filter_size")
//                 .long("filter_size")
//                 .takes_value(true)
//                 .default_value("10")
//                 .help("Unit with coverage less than this value would be discarded."),
//         )
//         .arg(
//             Arg::new("consensus_size")
//                 .long("consensus_size")
//                 .takes_value(true)
//                 .default_value("20")
//                 .help("The number of string to take consensus"),
//         )
// }

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
                .default_value("17"),
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

// fn subcommand_correct_multiplicity() -> Command<'static> {
//     Command::new("correct_multiplicity")
//         .version("0.1")
//         .author("Bansho Masutani")
//         .about("Fix multiplicities of units, re-clustering if needed.")
//         .arg(Arg::new("verbose").short('v').multiple_occurrences(true))
//         .arg(
//             Arg::new("threads")
//                 .short('t')
//                 .long("threads")
//                 .required(false)
//                 .value_name("THREADS")
//                 .help("Number of Threads")
//                 .default_value("1")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("seed")
//                 .short('s')
//                 .long("seed")
//                 .required(false)
//                 .value_name("SEED")
//                 .help("Seed for pseudorandon number generators.")
//                 .default_value("24")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("draft_assembly")
//                 .short('o')
//                 .long("draft_assembly")
//                 .required(false)
//                 .value_name("PATH")
//                 .help("If given, output draft GFA to PATH.")
//                 .takes_value(true),
//         )
// }

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

fn subcommand_squish() -> Command<'static> {
    Command::new("squish")
        .version("0.1")
        .author("BanshoMasutani")
        .about("Squish unreliable clusterings (Local).")
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
            Arg::new("supress_level")
                .short('s')
                .long("supress")
                .required(false)
                .value_name("THR")
                .help("Supression level from 0(=do not supress) to 1=(squish all the clusters).")
                .default_value("0.4")
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

// fn subcommand_partition_global() -> Command<'static> {
//     Command::new("partition_global")
//         .version("0.1")
//         .author("BanshoMasutani")
//         .about("Clustering reads (Global).")
//         .arg(
//             Arg::new("verbose")
//                 .short('v')
//                 .multiple_occurrences(true)
//                 .help("Debug mode"),
//         )
//         .arg(
//             Arg::new("graph")
//                 .long("graph")
//                 .help("Invoke graph-WhatsHap instead of de Bruijn."),
//         )
//         .arg(
//             Arg::new("threads")
//                 .short('t')
//                 .long("threads")
//                 .required(false)
//                 .value_name("THREADS")
//                 .help("Number of Threads")
//                 .default_value("1")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("k")
//                 .short('k')
//                 .long("kmer_size")
//                 .required(false)
//                 .value_name("KMER_SIZE")
//                 .help("The size of the kmer")
//                 .default_value("4")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("min_cluster_size")
//                 .short('m')
//                 .long("min_cluster_size")
//                 .required(false)
//                 .value_name("MIN_CLUSTER_SIZE")
//                 .help("The minimum size of a cluster")
//                 .default_value("50")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("mat_score")
//                 .short('p')
//                 .long("match_score")
//                 .required(false)
//                 .value_name("MATCH_SCORE")
//                 .help("The match score")
//                 .default_value("1")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("mismat_score")
//                 .short('q')
//                 .long("mismatch_score")
//                 .required(false)
//                 .value_name("MISMATCH_SCORE")
//                 .help("The mismatch score")
//                 .default_value("1")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("gap_score")
//                 .short('g')
//                 .long("gap_score")
//                 .required(false)
//                 .value_name("GAP_SCORE")
//                 .help("The gap penalty")
//                 .default_value("2")
//                 .takes_value(true),
//         )
// }

// fn subcommand_resolve_tangle() -> Command<'static> {
//     Command::new("resolve_tangle")
//         .version("0.1")
//         .author("BanshoMasutani")
//         .about("Resolve tangle by re-estimate copy numbers.")
//         .arg(
//             Arg::new("verbose")
//                 .short('v')
//                 .multiple_occurrences(true)
//                 .help("Debug mode"),
//         )
//         .arg(
//             Arg::new("threads")
//                 .short('t')
//                 .long("threads")
//                 .required(false)
//                 .value_name("THREADS")
//                 .help("Number of Threads")
//                 .default_value("1")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("repeat_num")
//                 .short('r')
//                 .long("repeat_num")
//                 .required(false)
//                 .value_name("REPEAT_NUM")
//                 .help("Do EM algorithm for REPEAT_NUM times.")
//                 .default_value("7")
//                 .takes_value(true),
//         )
//         .arg(
//             Arg::new("coverage_threshold")
//                 .short('x')
//                 .long("threshold")
//                 .required(false)
//                 .value_name("THRESHOLD")
//                 .help("Unit with less that this coverage would be ignored.")
//                 .default_value("5")
//                 .takes_value(true),
//         )
// }

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
            Arg::new("output")
                .short('o')
                .long("output")
                .required(true)
                .value_name("PATH")
                .help("Output file name")
                .takes_value(true),
        )
}

fn entry(matches: &clap::ArgMatches) -> std::io::Result<DataSet> {
    use haplotyper::entry::Entry;
    debug!("START\tEntry");
    set_threads(matches);
    let file = matches.value_of("input").unwrap();
    // let slag = matches.value_of("slag");
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
    use haplotyper::extract::{Extract, ExtractTarget};
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
                writeln!(wtr, "{}\t{}\t{}", asn, name, desc)?;
            }
        }
        &_ => unreachable!(),
    };
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
    set_threads(matches);
    use haplotyper::determine_units::{DetermineUnit, UnitConfig};
    use ReadType::*;
    let (cl, tn) = (chunk_len, take_num);
    let config = match dataset.read_type {
        CCS => UnitConfig::new_ccs(cl, tn, skip_len, margin, thrds, filter, upper, lower),
        CLR => UnitConfig::new_clr(cl, tn, skip_len, margin, thrds, filter, upper, lower),
        _ => UnitConfig::new_ont(cl, tn, skip_len, margin, thrds, filter, upper, lower),
    };
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
    // TODO:Branching by read type.
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

fn polish_unit(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tpolishing units");
    set_threads(matches);
    let filter_size: usize = matches
        .value_of("filter_size")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    let consensus_size: usize = matches
        .value_of("consensus_size")
        .and_then(|e| e.parse::<usize>().ok())
        .unwrap();
    use haplotyper::polish_units::*;
    let config = PolishUnitConfig::new(dataset.read_type, filter_size, consensus_size);
    dataset.polish_unit(&config);
}
fn multiplicity_estimation(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tmultiplicity estimation");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let seed: u64 = matches
        .value_of("seed")
        .and_then(|e| e.parse().ok())
        .unwrap();
    let cov: Option<f64> = matches.value_of("coverage").and_then(|e| e.parse().ok());
    let path = matches.value_of("draft_assembly");
    set_threads(matches);
    use haplotyper::multiplicity_estimation::*;
    let config = MultiplicityEstimationConfig::new(threads, seed, cov, path);
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

fn squish(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    set_threads(matches);
    let supress: f64 = matches
        .value_of("supress_level")
        .and_then(|e| e.parse().ok())
        .unwrap();
    use haplotyper::squish_erroneous_clusters::*;
    let config = SquishConfig::new(supress, 5);
    dataset.squish_erroneous_clusters(&config);
}

fn purge_diverged(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tPurge diverged clusters");
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    set_threads(matches);
    use haplotyper::purge_diverged::*;
    let config = PurgeDivConfig::new(threads);
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

// fn correct_multiplicity(matches: &clap::ArgMatches, dataset: &mut DataSet) {
//     debug!("START\tmultiplicity estimation");
//     let threads: usize = matches
//         .value_of("threads")
//         .and_then(|e| e.parse().ok())
//         .unwrap();
//     let seed: u64 = matches
//         .value_of("seed")
//         .and_then(|e| e.parse().ok())
//         .unwrap();
//     let path = matches.value_of("draft_assembly");
//     if let Err(why) = rayon::ThreadPoolBuilder::new()
//         .num_threads(threads)
//         .build_global()
//     {
//         debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
//     }
//     // Re-estimate the copy number, retry if needed.
//     use std::collections::{HashMap, HashSet};
//     let selection: HashSet<_> = {
//         let re_estimated_cluster_num: HashMap<_, _> = {
//             // TODO: This is very inefficient.
//             let mut ds: DataSet = dataset.clone();
//             use haplotyper::dirichlet_mixture::{ClusteringConfig, DirichletMixtureCorrection};
//             use haplotyper::multiplicity_estimation::*;
//             let config = ClusteringConfig::new(5, 10, 5);
//             ds.correct_clustering(&config);
//             let config = MultiplicityEstimationConfig::new(threads, seed, ds.coverage, path);
//             ds.estimate_multiplicity(&config);
//             ds.selected_chunks
//                 .iter()
//                 .map(|c| (c.id, c.copy_num))
//                 .collect()
//         };
//         dataset
//             .selected_chunks
//             .iter_mut()
//             .filter_map(|chunk| match re_estimated_cluster_num.get(&chunk.id) {
//                 Some(&new) if new != chunk.copy_num => {
//                     debug!("FIXMULTP\t{}\t{}\t{}", chunk.id, chunk.copy_num, new);
//                     chunk.copy_num = new;
//                     Some(chunk.id)
//                 }
//                 _ => None,
//             })
//             .collect()
//     };
//     // Re clustering.
//     debug!("FIXMULTP\tTargetNum\t{}", selection.len());
//     use haplotyper::local_clustering::*;
//     local_clustering_selected(dataset, &selection);
// }

// fn global_clustering(matches: &clap::ArgMatches, dataset: &mut DataSet) {
//     debug!("START\tGlobal Clustering step");
//     let threads: usize = matches
//         .value_of("threads")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     let kmer: usize = matches
//         .value_of("k")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     let min_cluster_size = matches
//         .value_of("min_cluster_size")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     let mat_score: i32 = matches
//         .value_of("mat_score")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     // Minus.
//     let mismat_score: i32 = -matches
//         .value_of("mismat_score")
//         .and_then(|num| num.parse::<i32>().ok())
//         .unwrap();
//     let gap_score: i32 = -matches
//         .value_of("gap_score")
//         .and_then(|num| num.parse::<i32>().ok())
//         .unwrap();
//     if let Err(why) = rayon::ThreadPoolBuilder::new()
//         .num_threads(threads)
//         .build_global()
//     {
//         debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
//     }
//     use haplotyper::global_clustering::*;
//     let config =
//         GlobalClusteringConfig::new(kmer, min_cluster_size, mat_score, mismat_score, gap_score);
//     if matches.is_present("graph") {
//         dataset.global_clustering_graph(&config);
//     } else {
//         dataset.global_clustering(&config);
//     }
// }

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

// fn resolve_tangle(matches: &clap::ArgMatches, _dataset: &mut DataSet) {
//     debug!("START\tClustering Correction step");
//     let threads: usize = matches
//         .value_of("threads")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     let _repeat_num: usize = matches
//         .value_of("repeat_num")
//         .and_then(|num| num.parse::<usize>().ok())
//         .unwrap();
//     let _threshold: usize = matches
//         .value_of("coverage_threshold")
//         .and_then(|num| num.parse().ok())
//         .unwrap();
//     if let Err(why) = rayon::ThreadPoolBuilder::new()
//         .num_threads(threads)
//         .build_global()
//     {
//         debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
//     }
//     todo!()
// }

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
    let threads: usize = matches
        .value_of("threads")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let window_size: usize = matches
        .value_of("window_size")
        .and_then(|num| num.parse().ok())
        .unwrap();
    let skip_polish = matches.is_present("no_polish");
    let file = matches.value_of("output").unwrap();
    let mut file = std::fs::File::create(file).map(BufWriter::new)?;
    use haplotyper::assemble::*;
    let msr = dataset.read_type.min_span_reads();
    let min_lk = dataset.read_type.min_llr_value();
    let config = AssembleConfig::new(threads, window_size, !skip_polish, true, msr, min_lk);
    debug!("START\tFinal assembly");
    let gfa = dataset.assemble(&config);
    writeln!(file, "{}", gfa)?;
    Ok(())
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
        .author("Bansho Masutani")
        .about("HLA toolchain")
        .arg_required_else_help(true)
        .subcommand(subcommand_entry())
        .subcommand(subcommand_extract())
        .subcommand(subcommand_stats())
        .subcommand(subcommand_select_unit())
        .subcommand(subcommand_encode())
        .subcommand(subcommand_polish_encoding())
        .subcommand(subcommand_estimate_multiplicity())
        .subcommand(subcommand_partition_local())
        .subcommand(subcommand_squish())
        .subcommand(subcommand_purge_diverged())
        .subcommand(subcommand_correct_deletion())
        .subcommand(subcommand_correct_clustering())
        .subcommand(subcommand_encode_densely())
        .subcommand(subcommand_assemble())
        .subcommand(subcommand_pick_components())
        .subcommand(subcommand_mask_repeats())
        .get_matches();
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
    let mut ds = get_input_file()?;
    let ds = &mut ds;
    match matches.subcommand() {
        Some(("select_unit", sub_m)) => select_unit(sub_m, ds),
        Some(("mask_repeats", sub_m)) => repeat_masking(sub_m, ds),
        Some(("encode", sub_m)) => encode(sub_m, ds),
        Some(("polish_unit", sub_m)) => polish_unit(sub_m, ds),
        Some(("pick_components", sub_m)) => pick_components(sub_m, ds),
        Some(("polish_encoding", sub_m)) => polish_encode(sub_m, ds),
        Some(("partition_local", sub_m)) => local_clustering(sub_m, ds),
        Some(("squish", sub_m)) => squish(sub_m, ds),
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
