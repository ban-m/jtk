use definitions::*;
use haplotyper::model_tune::ModelFit;
use std::io::BufReader;
use std::io::{BufWriter, Write};
#[macro_use]
extern crate log;

fn main() -> std::io::Result<()> {
    let matches = jtk_cli::jtk_commands::jtk_parser().get_matches();
    if let Some(("pipeline", sub_m)) = matches.subcommand() {
        let path: &String = sub_m.get_one("profile").unwrap();
        use std::io::Read;
        let mut rdr = std::fs::File::open(path).map(std::io::BufReader::new)?;
        let mut file = String::new();
        rdr.read_to_string(&mut file)?;
        let config: jtk_cli::pipeline::PipelineConfig = toml::from_str(&file).unwrap();
        return jtk_cli::pipeline::run_pipeline(&config);
    }
    if let Some((_, sub_m)) = matches.subcommand() {
        let level = match sub_m.get_count("verbose") {
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
        Some(("select_chunks", sub_m)) => select_chunks(sub_m, ds),
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
        Some(("squish", sub_m)) => squish(sub_m, ds),
        _ => unreachable!(),
    };
    flush_file(ds)
}

fn entry(matches: &clap::ArgMatches) -> std::io::Result<DataSet> {
    use haplotyper::entry::Entry;
    debug!("START\tEntry");
    set_threads(matches);
    let file: &String = matches.get_one("input").unwrap();
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
    let read_type: &String = matches.get_one("read_type").unwrap();
    Ok(DataSet::entry(file, seqs, read_type))
}

fn extract(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    use haplotyper::extract::Extract;
    debug!("START\tExtract");
    debug!("Target is {}", matches.get_one::<String>("target").unwrap());
    let file = std::fs::File::create(matches.get_one::<String>("output").unwrap())?;
    let mut wtr = std::io::BufWriter::new(file);
    dataset.extract(&mut wtr)?;
    Ok(())
}

fn stats(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    use haplotyper::stats::Stats;
    debug!("START\tStats step");
    let wtr = matches.get_one::<String>("file").unwrap().as_str();
    let wtr = std::io::BufWriter::new(std::fs::File::create(wtr)?);
    dataset.sanity_check();
    dataset.stats(wtr)?;
    Ok(())
}

fn select_chunks(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tSelecting Units");
    let chunk_len: usize = matches
        .get_one("chunk_len")
        .and_then(|e: &String| e.parse().ok())
        .expect("Chunk len");
    let margin: usize = matches
        .get_one("margin")
        .and_then(|e: &String| e.parse().ok())
        .expect("Margin");
    let take_num: usize = matches
        .get_one("take_num")
        .and_then(|e: &String| e.parse().ok())
        .expect("Take num");
    let thrds: usize = matches
        .get_one("threads")
        .and_then(|e: &String| e.parse().ok())
        .expect("threads");
    let filter: f64 = matches
        .get_one("exclude")
        .and_then(|e: &String| e.parse().ok())
        .expect("exclude");
    let purge_copy_num: usize = matches
        .get_one("purge_copy_num")
        .and_then(|e: &String| e.parse::<usize>().ok())
        .unwrap();
    let seed: u64 = matches
        .get_one("seed")
        .and_then(|e: &String| e.parse::<u64>().ok())
        .unwrap();
    set_threads(matches);
    use haplotyper::determine_chunks::{DetermineUnit, DetermineUnitConfig};
    let (cl, tn) = (chunk_len, take_num);
    let config = DetermineUnitConfig::new(cl, tn, margin, thrds, filter, purge_copy_num, seed);
    dataset.select_chunks(&config);
}

fn repeat_masking(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tmasking repeat.");
    set_threads(matches);
    let k: usize = matches
        .get_one("k")
        .and_then(|l: &String| l.parse().ok())
        .unwrap();
    if k > 32 {
        panic!("K should be less than 32.");
    }
    let freq: f64 = matches
        .get_one("freq")
        .and_then(|l: &String| l.parse().ok())
        .unwrap();
    let min: u32 = matches
        .get_one("min")
        .and_then(|l: &String| l.parse().ok())
        .unwrap();
    use haplotyper::repeat_masking::*;
    let config = RepeatMaskConfig::new(k, freq, min);
    dataset.mask_repeat(&config)
}

fn encode(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tEncoding step");
    set_threads(matches);
    let threads: usize = matches
        .get_one("threads")
        .and_then(|e: &String| e.parse::<usize>().ok())
        .unwrap();
    use haplotyper::encode::Encode;
    let rt = dataset.read_type;
    let sim_thr = match matches
        .get_one("sim_thr")
        .and_then(|e: &String| e.parse().ok())
    {
        Some(res) => res,
        None => rt.sim_thr(),
    };
    dataset.encode(threads, sim_thr, rt.sd_of_error())
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
        .get_one("component_num")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    use haplotyper::pick_component::*;
    let config = ComponentPickingConfig::new(component_num);
    dataset.pick_top_n_component(&config);
}
fn multiplicity_estimation(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tmultiplicity estimation");
    let seed: u64 = matches
        .get_one("seed")
        .and_then(|e: &String| e.parse().ok())
        .unwrap();
    let cov: Option<f64> = matches
        .get_one("coverage")
        .and_then(|e: &String| e.parse().ok());
    dataset.coverage = match cov {
        Some(cov) => definitions::Coverage::Protected(cov),
        None => definitions::Coverage::NotAvailable,
    };
    let path: Option<&String> = matches.get_one("draft_assembly");
    let path = path.map(|x| x.as_str());
    set_threads(matches);
    use haplotyper::multiplicity_estimation::*;
    let purge: Option<usize> = matches
        .get_one("purge")
        .and_then(|x: &String| x.parse().ok());
    let config = MultiplicityEstimationConfig::new(seed, path);
    dataset.estimate_multiplicity(&config);
    if let Some(upper) = purge {
        dataset.purge_multiplicity(upper);
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
    let to_recal = matches.get_flag("re_cluster");
    use haplotyper::determine_chunks::calc_sim_thr;
    use haplotyper::determine_chunks::TAKE_THR;
    let sim_thr = calc_sim_thr(dataset, TAKE_THR);
    use haplotyper::encode::deletion_fill::*;
    let rt = dataset.read_type;
    let config = CorrectDeletionConfig::new(to_recal, Some(sim_thr), Some(rt.sd_of_error()));
    debug!("SIMTHR\t{sim_thr}");
    dataset.correct_deletion(&config);
}

fn clustering_correction(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tClustering Correction step");
    let _repeat_num: usize = matches
        .get_one("repeat_num")
        .and_then(|num: &String| num.parse::<usize>().ok())
        .unwrap();
    let _threshold: usize = matches
        .get_one("coverage_threshold")
        .and_then(|num: &String| num.parse().ok())
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
        .get_one("length")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    use haplotyper::dense_encoding::*;
    let file: Option<&String> = matches.get_one("output");
    let file = file.map(|x| x.as_str());
    let config = DenseEncodingConfig::new(length, file);
    dataset.dense_encoding(&config);
}
fn squish(matches: &clap::ArgMatches, dataset: &mut DataSet) {
    debug!("START\tSquish Clustering");
    set_threads(matches);
    let ari_thr: f64 = matches
        .get_one("ari")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let count: usize = matches
        .get_one("count")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let match_score: f64 = matches
        .get_one("match_score")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let mismatch_score: f64 = matches
        .get_one("mismatch_score")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    use haplotyper::{SquishConfig, SquishErroneousClusters};
    let config = SquishConfig::new(ari_thr, count, match_score, mismatch_score);
    dataset.squish_erroneous_clusters(&config);
}

fn assembly(matches: &clap::ArgMatches, dataset: &mut DataSet) -> std::io::Result<()> {
    debug!("START\tAssembly step");
    set_threads(matches);
    let window_size: usize = matches
        .get_one("window_size")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let min_llr: f64 = matches
        .get_one("min_llr")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let min_span: usize = matches
        .get_one("min_span")
        .and_then(|num: &String| num.parse().ok())
        .unwrap();
    let skip_polish = matches.get_flag("no_polish");
    let file: Option<&String> = matches.get_one("output");
    let file = file.map(|x| x.as_str()).unwrap();
    use haplotyper::assemble::*;
    let config = AssembleConfig::new(
        window_size,
        !skip_polish,
        true,
        min_span,
        min_llr,
        true,
        Some(file),
    );
    debug!("START\tFinal assembly");
    if !skip_polish {
        dataset.fit_models_on_both_strands();
    }
    let gfa = dataset.assemble(&config);
    let mut file = std::fs::File::create(format!("{file}.gfa")).map(BufWriter::new)?;
    writeln!(file, "{}", gfa)?;
    Ok(())
}

fn polish(matches: &clap::ArgMatches) -> std::io::Result<()> {
    set_threads(matches);
    let window_size: usize = matches
        .get_one("window_size")
        .and_then(|x: &String| x.parse().ok())
        .unwrap();
    let reads: &String = matches.get_one("reads").unwrap();
    let contig: &String = matches.get_one("contigs").unwrap();
    let alignments: &String = matches.get_one("alignments").unwrap();
    let format: &String = matches.get_one("format").unwrap();
    let seed: u64 = matches
        .get_one("seed")
        .and_then(|x: &String| x.parse().ok())
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

fn set_threads(matches: &clap::ArgMatches) {
    if let Some(threads) = matches
        .get_one("threads")
        .and_then(|num: &String| num.parse().ok())
    {
        debug!("Set Threads\t{}", threads);
        if let Err(why) = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
        {
            debug!("{:?} If you run `pipeline` module, this is Harmless.", why);
        }
    }
}
