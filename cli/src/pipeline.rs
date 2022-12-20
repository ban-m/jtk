use definitions::DataSet;
use serde::{Deserialize, Serialize};
extern crate log;
use log::*;
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct PipelineConfig {
    input_file: String,
    read_type: String,
    out_dir: String,
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
    let file_stem = format!("{out_dir}/{prefix}");
    let entry = format!("{file_stem}.encoded.json");
    let encoded = format!("{file_stem}.encoded.json");
    let clustered = format!("{file_stem}.clustered.json");
    let dense_encoded = format!("{file_stem}.de.json");
    let corrected = format!("{file_stem}.json");
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
    let take_num = genome_size / chunk_len;
    let repeat_mask_config = RepeatMaskConfig::new(kmersize, top_freq, min_count);
    let select_unit_config = DetermineUnitConfig::new(
        chunk_len,
        take_num,
        margin,
        threads,
        exclude,
        purge_copy_num,
        seed,
    );
    let pick_component_config = ComponentPickingConfig::new(component_num);
    let draft = format!("{file_stem}.draft.gfa");
    let multp_config = MultiplicityEstimationConfig::new(seed, Some(&draft));
    let purge_config = PurgeDivConfig::new();
    let de = format!("{file_stem}.draft2.gfa");
    let dense_encode_config = DenseEncodingConfig::new(compress_contig, Some(&de));
    let correction_config = CorrectionConfig::default();
    use haplotyper::determine_units::STDDEV_OR_ERROR;
    let dump = Some(file_stem.as_str());
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
    let mut ds = match resume && matches!(std::fs::try_exists(&entry), Ok(true)) {
        false => {
            let mut ds = parse_input(&input_file, &read_type)?;
            if let Some(hap) = haploid_coverage {
                ds.coverage = definitions::Coverage::Protected(hap);
            }
            log(&ds, &entry)?;
            ds
        }
        true => parse_input(&input_file, &read_type)?,
    };
    if resume && matches!(std::fs::try_exists(&encoded), Ok(true)) {
        ds = parse_json(&encoded)?
    } else {
        ds.mask_repeat(&repeat_mask_config);
        ds.select_chunks(&select_unit_config);
        ds.pick_top_n_component(&pick_component_config);
        ds.correct_deletion(&correct_deletion_config);
        ds.remove_erroneous_nodes();
        ds.estimate_multiplicity(&multp_config);
        ds.purge_multiplicity(purge_copy_num);
        log(&ds, &encoded)?;
    }
    if resume && matches!(std::fs::try_exists(&clustered), Ok(true)) {
        ds = parse_json(&clustered)?
    } else {
        ds.local_clustering();
        log(&ds, &clustered)?;
    }
    if resume && matches!(std::fs::try_exists(&dense_encoded), Ok(true)) {
        ds = parse_json(&dense_encoded)?;
    } else {
        ds.purge(&purge_config);
        ds.purge(&purge_config);
        ds.correct_deletion(&correct_deletion_config_recluster);
        ds.dense_encoding(&dense_encode_config);
        ds.correct_deletion(&correct_deletion_config_recluster);
        log(&ds, &dense_encoded)?;
    }
    if resume && matches!(std::fs::try_exists(&corrected), Ok(true)) {
        ds = parse_json(&corrected)?;
    } else {
        ds.squish_erroneous_clusters(&squish_config);
        ds.correct_clustering(&correction_config);
        log(&ds, &corrected)?;
    }
    // Flush the result.
    let gfa = ds.assemble(&assemble_config);
    let mut asm_file = std::fs::File::create(format!("{file_stem}.gfa")).map(BufWriter::new)?;
    writeln!(asm_file, "{gfa}")
}

fn parse_json(filename: &str) -> std::io::Result<DataSet> {
    debug!("RESUME\t{filename}");
    std::fs::File::open(filename)
        .map(std::io::BufReader::new)
        .map(serde_json::de::from_reader)
        .map(|x| x.unwrap())
}

fn log(ds: &DataSet, path: &str) -> std::io::Result<()> {
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
