//! Pipelines -- the whole pipeline of the JTK.
//!
//! This module defines the pipeline of the JTK to assemble a genomic region in a diploid resolution.
use definitions::DataSet;
use definitions::ReadType;
use serde::{Deserialize, Serialize};
extern crate log;
use log::*;
use std::path::Path;
use std::path::PathBuf;

/// The configuration of the pipeline.
/// This struct is a comprehensive list of the parameters that can be
/// set by a user. All other parameters in the pipeline or algorithm would be determined automatically or
/// hard-coded to the values that work well for most of the case.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PipelineConfig {
    /// The path to the input file.
    input_file: PathBuf,
    /// The type of the read. One of the
    read_type: ReadType,
    /// The path to the output directory.
    out_dir: PathBuf,
    prefix: String,
    verbose: usize,
    threads: usize,
    seed: u64,
    region_size: String,
    chunk_len: usize,
    margin: usize,
    exclude: f64,
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
    supress_ari: f64,
    match_ari: f64,
    mismatch_ari: f64,
    required_count: usize,
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
        region_size,
        margin,
        exclude,
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
        supress_ari,
        match_ari,
        mismatch_ari,
        required_count,
    } = config.clone();
    let level = match verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(level)).init();
    let file_stem = out_dir.join(prefix);
    let entry = file_stem.with_extension("entry.json");
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
    assert!(std::path::Path::new(&out_dir).is_dir());
    // Configurations.
    let genome_size = match parse_si(&region_size) {
        Ok(res) => res,
        Err(why) => {
            error!("{} can not be parsed!", region_size);
            panic!("{}", why);
        }
    };
    let take_num = 3 * genome_size / chunk_len / 2;
    let repeat_mask_config = RepeatMaskConfig::new(kmersize, top_freq, min_count);
    let select_chunk_config = DetermineUnitConfig::new(
        chunk_len,
        take_num,
        margin,
        threads,
        exclude,
        purge_copy_num,
        seed,
    );
    let pick_component_config = ComponentPickingConfig::new(component_num);
    let draft = file_stem.with_extension("draft.gfa");
    let multp_config = MultiplicityEstimationConfig::new(seed, Some(&draft));
    let purge_config = PurgeDivConfig::new();
    let de = file_stem.with_extension("draft2.gfa");
    let dense_encode_config = DenseEncodingConfig::new(compress_contig, Some(&de));
    let correction_config = CorrectionConfig::default();
    use haplotyper::determine_chunks::STDDEV_OR_ERROR;
    // TODO: Make this to as-is.
    let dump = file_stem.as_os_str().to_str();
    let assemble_config = AssembleConfig::new(
        polish_window_size,
        to_polish,
        true,
        min_span,
        min_llr,
        true,
        dump,
    );
    let correct_deletion_config = CorrectDeletionConfig::new(false, None, Some(STDDEV_OR_ERROR));
    let correct_deletion_config_recluster =
        CorrectDeletionConfig::new(true, None, Some(STDDEV_OR_ERROR));
    let squish_config = SquishConfig::new(supress_ari, required_count, match_ari, mismatch_ari);
    // Pipeline.
    let mut ds = match resume && matches!(std::path::Path::new(&entry).try_exists(), Ok(true)) {
        false => {
            let mut ds = parse_input(&input_file, read_type)?;
            if let Some(hap) = haploid_coverage {
                ds.coverage = definitions::Coverage::Protected(hap);
            }
            log(&ds, &entry)?;
            ds
        }
        true => parse_input(&input_file, read_type)?,
    };
    if resume && matches!(std::path::Path::new(&encoded).try_exists(), Ok(true)) {
        ds = parse_json(&encoded)?
    } else {
        ds.mask_repeat(&repeat_mask_config);
        ds.select_chunks(&select_chunk_config);
        ds.pick_top_n_component(&pick_component_config);
        ds.correct_deletion(&correct_deletion_config);
        ds.remove_erroneous_nodes();
        ds.estimate_multiplicity(&multp_config);
        ds.purge_multiplicity(purge_copy_num);
        log(&ds, &encoded)?;
    }
    if resume && matches!(std::path::Path::new(&clustered).try_exists(), Ok(true)) {
        ds = parse_json(&clustered)?
    } else {
        ds.local_clustering();
        log(&ds, &clustered)?;
    }
    if resume && matches!(std::path::Path::new(&dense_encoded).try_exists(), Ok(true)) {
        ds = parse_json(&dense_encoded)?;
    } else {
        ds.purge(&purge_config);
        ds.purge(&purge_config);
        ds.correct_deletion(&correct_deletion_config_recluster);
        ds.dense_encoding(&dense_encode_config);
        ds.correct_deletion(&correct_deletion_config_recluster);
        log(&ds, &dense_encoded)?;
    }
    if resume && matches!(std::path::Path::new(&corrected).try_exists(), Ok(true)) {
        ds = parse_json(&corrected)?;
    } else {
        ds.squish_erroneous_clusters(&squish_config);
        ds.correct_clustering(&correction_config);
        log(&ds, &corrected)?;
    }
    // Flush the result.
    let gfa = ds.assemble(&assemble_config);
    let mut asm_file =
        std::fs::File::create(file_stem.clone().with_extension("gfa")).map(BufWriter::new)?;
    writeln!(asm_file, "{gfa}")
}

fn parse_json(filename: &Path) -> std::io::Result<DataSet> {
    debug!("RESUME\t{filename:?}");
    std::fs::File::open(filename)
        .map(std::io::BufReader::new)
        .map(serde_json::de::from_reader)
        .map(|x| x.unwrap())
}

fn log(ds: &DataSet, path: &Path) -> std::io::Result<()> {
    let mut wtr = std::fs::File::create(path).map(BufWriter::new)?;
    serde_json::ser::to_writer(&mut wtr, ds).unwrap();
    Ok(())
}

fn parse_input(input_file: &PathBuf, read_type: ReadType) -> std::io::Result<DataSet> {
    debug!("Opening {:?}", input_file);
    let reader = std::fs::File::open(input_file).map(BufReader::new)?;
    let extension = input_file.extension().unwrap().to_str().unwrap();
    let is_fasta = extension.ends_with('a');
    let is_fastq = extension.ends_with('q');
    let seqs: Vec<(String, Vec<u8>)> = if is_fasta {
        bio_utils::fasta::parse_into_vec_from(reader)?
            .into_iter()
            .map(|records| {
                let (id, _, seq) = records.into();
                let seq = seq.into_bytes();
                (id, seq)
            })
            .collect()
    } else if is_fastq {
        bio_utils::fastq::parse_into_vec_from(reader)?
            .into_iter()
            .map(|record| {
                if record.seq().iter().any(|x| !b"ACGT".contains(x)) {
                    debug!("{}", record);
                }
                let (id, seq, _) = record.into();
                (id, seq)
            })
            .collect()
    } else {
        panic!("file type:{:?} not supported", input_file)
    };
    Ok(DataSet::entry(input_file, seqs, read_type))
}

fn parse_si(input: &str) -> Result<usize, std::num::ParseFloatError> {
    assert!(!input.is_empty(), "Please specify genome size.");
    let mut input = input.to_string();
    let last = input.chars().last().unwrap();
    let mult = match last {
        'k' | 'K' => 1_000,
        'm' | 'M' => 1_000_000,
        'g' | 'G' => 1_000_000_000,
        '0'..='9' => 1,
        _ => panic!("si prefix {} is not supported yet.", last),
    };
    if last.is_ascii_alphabetic() {
        input.pop();
    }
    let number = input.parse::<f64>();
    number.map(|x| (x * mult as f64).round() as usize)
}
