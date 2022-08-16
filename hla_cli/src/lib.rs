use definitions::DataSet;
use serde::{Deserialize, Serialize};
extern crate log;
use log::*;
use std::path::PathBuf;
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PipelineConfig {
    input_file: String,
    read_type: String,
    out_dir: std::path::PathBuf,
    prefix: String,
    verbose: usize,
    threads: usize,
    seed: u64,
    chunk_len: usize,
    take_num: usize,
    margin: usize,
    exclude: f64,
    upper: usize,
    lower: usize,
    kmersize: usize,
    top_freq: f64,
    min_count: u32,
    component_num: usize,
    purge_copy_num: usize,
    haploid_coverage: Option<f64>,
    compress_contig: usize,
    polish_window_size: usize,
    to_polish: bool,
    min_span: usize,
    min_llr: f64,
    resume: bool,
}
use haplotyper::{local_clustering::LocalClustering, *};
use std::io::{BufReader, BufWriter, Write};
pub fn run_pipeline(config: &PipelineConfig) -> std::io::Result<()> {
    let PipelineConfig {
        input_file,
        read_type,
        out_dir,
        prefix,
        threads,
        seed,
        verbose,
        chunk_len,
        take_num,
        margin,
        exclude,
        upper,
        lower,
        kmersize,
        top_freq,
        min_count,
        component_num,
        purge_copy_num,
        haploid_coverage,
        compress_contig,
        polish_window_size,
        to_polish,
        min_span,
        min_llr,
        resume,
    } = config.clone();
    let level = match verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let file_stem = out_dir.join(prefix);
    let encoded = file_stem.with_extension("encoded.json");
    let clustered = file_stem.with_extension("clustered.json");
    let dense_encoded = file_stem.with_extension("de.json");
    let corrected = file_stem.with_extension("json");
    assert!(kmersize < 32);
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    std::fs::create_dir_all(&out_dir)?;
    assert!(out_dir.is_dir());
    // Configurations.
    let repeat_mask_config = RepeatMaskConfig::new(kmersize, top_freq, min_count);
    let select_unit_config =
        DetermineUnitConfig::new(chunk_len, take_num, margin, threads, exclude, upper, lower);
    let pick_component_config = ComponentPickingConfig::new(component_num);
    let draft = file_stem.with_extension("draft.gfa");
    let multp_config = MultiplicityEstimationConfig::new(seed, draft.as_os_str().to_str());
    let purge_config = PurgeDivConfig::new();
    let de = file_stem.with_extension("draft2.gfa");
    let dense_encode_config = DenseEncodingConfig::new(compress_contig, de.as_os_str().to_str());
    let correction_config = CorrectionConfig::default();
    let assemble_config =
        AssembleConfig::new(polish_window_size, to_polish, true, min_span, min_llr, true);
    let mut asm_file =
        std::fs::File::create(file_stem.with_extension("gfa")).map(BufWriter::new)?;
    if resume {
        warn!("Currently resuming is not supported....");
    }
    // Pipeline.
    let mut ds = parse_input(&input_file, &read_type)?;
    if let Some(hap) = haploid_coverage {
        ds.coverage = definitions::Coverage::Protected(hap);
    }
    ds.mask_repeat(&repeat_mask_config);
    ds.select_chunks(&select_unit_config);
    ds.pick_top_n_component(&pick_component_config);
    let correct_deletion_config = CorrectDeletionConfig::new(false, get_sim_thr(&ds));
    ds.correct_deletion(&correct_deletion_config);
    ds.remove_erroneous_nodes();
    ds.estimate_multiplicity(&multp_config);
    ds.purge_multiplicity(purge_copy_num);
    log(&ds, encoded)?;
    ds.local_clustering();
    log(&ds, clustered)?;
    ds.purge(&purge_config);
    ds.purge(&purge_config);
    let correct_deletion_config = CorrectDeletionConfig::new(true, get_sim_thr(&ds));
    ds.correct_deletion(&correct_deletion_config);
    ds.dense_encoding_dev(&dense_encode_config);
    log(&ds, dense_encoded)?;
    ds.correct_deletion(&correct_deletion_config);
    ds.correct_clustering(&correction_config);
    let gfa = ds.assemble(&assemble_config);
    log(&ds, corrected)?;
    writeln!(asm_file, "{gfa}")
}

fn get_sim_thr(ds: &DataSet) -> f64 {
    haplotyper::determine_units::calc_sim_thr(&ds, haplotyper::determine_units::TAKE_THR)
}

fn log(ds: &DataSet, path: PathBuf) -> std::io::Result<()> {
    let mut wtr = std::fs::File::create(path).map(BufWriter::new)?;
    serde_json::ser::to_writer(&mut wtr, ds).unwrap();
    Ok(())
}

fn parse_input(input_file: &str, read_type: &str) -> std::io::Result<DataSet> {
    debug!("Opening {}", input_file);
    let reader = std::fs::File::open(input_file).map(BufReader::new)?;
    let seqs: Vec<(String, Vec<u8>)> = match input_file.chars().last() {
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
        _ => panic!("file type:{} not supported", input_file),
    };
    Ok(DataSet::entry(input_file, seqs, read_type))
}
