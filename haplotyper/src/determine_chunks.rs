use crate::ALN_PARAMETER;
const MIN_REQ_NEW_UNIT: usize = 10;
use super::encode::Encode;
use super::polish_chunks::PolishChunk;
use super::polish_chunks::PolishChunkConfig;
use definitions::*;
use log::*;
use rand::prelude::*;
use rand_xoshiro::Xoroshiro128Plus;
use rayon::prelude::*;
use std::collections::HashMap;

use std::hash::Hash;
#[derive(Debug, Clone)]
pub struct DetermineUnitConfig {
    pub chunk_len: usize,
    pub chunk_num: usize,
    pub margin: usize,
    pub threads: usize,
    pub min_cluster: usize,
    pub exclude_repeats: f64,
    pub purge_copy_num: usize,
    pub seed: u64,
}

pub const STDDEV_OR_ERROR: f64 = 0.01;
const FIRST_RELAX: f64 = 2f64;

pub const TAKE_THR: f64 = 0.999;
impl DetermineUnitConfig {
    pub fn new(
        chunk_len: usize,
        chunk_num: usize,
        margin: usize,
        threads: usize,
        exclude_repeats: f64,
        purge_copy_num: usize,
        seed: u64,
    ) -> Self {
        Self {
            chunk_len,
            margin,
            chunk_num,
            threads,
            min_cluster: 2,
            exclude_repeats,
            purge_copy_num,
            seed,
        }
    }
}

pub trait DetermineUnit {
    fn select_chunks(&mut self, config: &DetermineUnitConfig);
}

fn show_minimap2_version() {
    let version = std::process::Command::new("minimap2")
        .args(["--version"])
        .output()
        .expect("Minimap2 is unavailable.");
    if !version.status.success() {
        error!("Minimap2 is not available. Please install minimap2 (https://github.com/lh3/minimap2) first.");
        panic!("Minimap2,{:?}", String::from_utf8_lossy(&version.stderr));
    } else {
        let version = std::str::from_utf8(&version.stdout).unwrap();
        debug!("MINIMAP2\tVERSION\t{}", version);
    }
}
fn update_get_coverage(ds: &mut DataSet) -> f64 {
    crate::misc::update_coverage(ds);
    ds.coverage.unwrap()
}

const FIRST_CONS_COV: usize = 20;
const CONS_COV: usize = 30;
const COPY_NUM_OFFSET: usize = 3;
const LOWER_FRAC: f64 = 0.1;
impl DetermineUnit for definitions::DataSet {
    fn select_chunks(&mut self, config: &DetermineUnitConfig) {
        show_minimap2_version();
        self.selected_chunks.clear();
        self.encoded_reads.clear();
        debug!("Select Unit: Configuration:{:?}", config);
        let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(config.seed);
        self.selected_chunks = pick_random(self, config, &mut rng);
        debug!("UNITNUM\t{}\tPICKED", self.selected_chunks.len());
        let overlap_identity_thr = self.read_type.overlap_identity_thr();
        remove_overlapping_chunks(self, overlap_identity_thr, config).unwrap();
        compaction_chunks(self);
        // 1st polishing.
        use crate::repeat_masking::RepeatMask;
        use crate::stats::Stats;
        let repetitive_kmer = self.get_repetitive_kmer();
        let mut sim_thr = self.read_type.sim_thr();
        {
            debug!("UNITNUM\t{}\tREMOVED", self.selected_chunks.len());
            self.encode(config.threads, FIRST_RELAX * sim_thr, STDDEV_OR_ERROR);
            debug!("ERRORRATE\t{}", self.error_rate());
            let haploid_coverage = update_get_coverage(self);
            let upper_count =
                haploid_coverage.ceil() as usize * (config.purge_copy_num + COPY_NUM_OFFSET);
            let filter_size = (haploid_coverage * LOWER_FRAC).ceil() as usize;
            remove_frequent_chunks(self, upper_count);
            dump_histogram(self);
            let polish_config = PolishChunkConfig::new(self.read_type, filter_size, FIRST_CONS_COV);
            self.consensus_chunk(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED\t1", self.selected_chunks.len());
        }
        // 2nd polishing.
        {
            self.encode(config.threads, sim_thr, STDDEV_OR_ERROR);
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            let fill_config = crate::encode::deletion_fill::CorrectDeletionConfig::new(
                false,
                Some(sim_thr),
                Some(STDDEV_OR_ERROR),
            );
            for _ in 0..10 {
                let new_chunk = fill_sparse_region(self, &repetitive_kmer, config)
                    + fill_tips(self, &repetitive_kmer, config);
                crate::encode::deletion_fill::correct_chunk_deletion(self, &fill_config);
                if new_chunk < MIN_REQ_NEW_UNIT {
                    break;
                }
            }
            compaction_chunks(self);
            // Here we usually see almost no errors in the reference chunk,
            // thus using very conservative and threshold for removing reference chunks.
            let ovlp_thr = 0.95;
            let haploid_coverage = update_get_coverage(self);
            let upper_count =
                haploid_coverage.ceil() as usize * (config.purge_copy_num + COPY_NUM_OFFSET);
            let filter_size = (haploid_coverage * LOWER_FRAC).ceil() as usize;
            remove_overlapping_chunks(self, ovlp_thr, config).unwrap();
            remove_frequent_chunks(self, upper_count);
            filter_chunk_by_ovlp(self, config);
            debug!("UNITNUM\t{}\tFILTERED\t1", self.selected_chunks.len());
            let polish_config = PolishChunkConfig::new(self.read_type, filter_size, CONS_COV);
            dump_histogram(self);
            self.polish_chunk(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED\t2", self.selected_chunks.len());
        };
        // Final polish
        {
            debug!("UNITNUM\t{}\tRAWUNIT", self.selected_chunks.len());
            use crate::encode::encode_by_mm2;
            encode_by_mm2(self, config.threads, sim_thr).unwrap();
            let haploid_coverage = update_get_coverage(self);
            let upper_count =
                haploid_coverage.ceil() as usize * (config.purge_copy_num + COPY_NUM_OFFSET);
            let filter_size = (haploid_coverage * LOWER_FRAC).ceil() as usize;
            remove_frequent_chunks(self, upper_count);
            filter_chunk_by_ovlp(self, config);
            self.encode(config.threads, sim_thr, self.read_type.sd_of_error());
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            remove_frequent_chunks(self, upper_count);
            filter_chunk_by_ovlp(self, config);
            compaction_chunks(self);
            debug!("UNITNUM\t{}\tFILTERED\t2", self.selected_chunks.len());
            remove_frequent_chunks(self, upper_count);
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            let polish_config = PolishChunkConfig::new(self.read_type, 2 * filter_size, CONS_COV);
            dump_histogram(self);
            self.polish_chunk(&polish_config);
            self.selected_chunks
                .retain(|c| repetitive_kmer.repetitiveness(c.seq()) < config.exclude_repeats);
            debug!("UNITNUM\t{}\tPOLISHED\t3", self.selected_chunks.len());
        }
        {
            crate::misc::update_coverage(self);
            let upper_count = (self.coverage.unwrap().ceil() as usize)
                * (config.purge_copy_num + COPY_NUM_OFFSET);
            self.encode(config.threads, sim_thr, self.read_type.sd_of_error());
            debug!("ERRORRATE\t{}", self.error_rate());
            remove_frequent_chunks(self, upper_count);
            dump_histogram(self);
        }
        // If half of the coverage supports large deletion, remove them.
        const OCCUPY_FRACTION: f64 = 0.5;
        use crate::purge_diverged::*;
        let p_config = PurgeLargeDelConfig::new(crate::MAX_ALLOWED_GAP, OCCUPY_FRACTION, true);
        self.purge_largeindel(&p_config);
        compaction_chunks(self);
    }
}

fn remove_frequent_chunks(ds: &mut DataSet, upper_count: usize) {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.chunk).or_default() += 1;
    }
    counts.retain(|_, occ| *occ > upper_count);
    ds.selected_chunks.retain(|u| !counts.contains_key(&u.id));
    for read in ds.encoded_reads.iter_mut() {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if counts.contains_key(&node.chunk) => read.remove(idx),
                Some(_) => idx += 1,
                None => break,
            }
        }
    }
}

// Make the id consective
fn compaction_chunks(ds: &mut DataSet) {
    let mut mapping: HashMap<u64, u64> = HashMap::new();
    for (idx, chunk) in ds.selected_chunks.iter_mut().enumerate() {
        mapping.insert(chunk.id, idx as u64);
        chunk.id = idx as u64;
    }
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        read.nodes
            .iter_mut()
            .for_each(|n| n.chunk = mapping[&n.chunk]);
        read.edges.iter_mut().for_each(|e| {
            e.from = mapping[&e.from];
            e.to = mapping[&e.to];
        });
    });
}

use rand::Rng;
fn pick_random<R: Rng>(ds: &DataSet, config: &DetermineUnitConfig, rng: &mut R) -> Vec<Chunk> {
    use crate::RepeatMask;
    use rand::prelude::*;
    let mask = ds.get_repetitive_kmer();
    let subseqs: Vec<_> = ds
        .raw_reads
        .iter()
        .flat_map(|r| split_into(r, config))
        .map(|u| (u, mask.repetitiveness(u)))
        .filter(|&(_, repetitiveness)| repetitiveness < config.exclude_repeats)
        .collect();
    const SMALL: f64 = 0.0000000001;
    subseqs
        .choose_multiple_weighted(rng, config.chunk_num, |(_, repet)| {
            (1f64 - repet).max(SMALL)
        })
        .unwrap()
        .enumerate()
        .map(|(idx, (seq, _))| {
            let mut seq = seq.to_vec();
            seq.iter_mut().for_each(u8::make_ascii_uppercase);
            Chunk::new(idx as u64, seq, config.min_cluster)
        })
        .collect()
}

fn mm2_chunk_overlap(ds: &DataSet, config: &DetermineUnitConfig) -> std::io::Result<Vec<u8>> {
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut c_dir = std::env::current_dir()?;
    c_dir.push(format!("{}", id));
    debug!("DETERMINE\tUnitOvlp\tCreating\t{:?}", c_dir);
    std::fs::create_dir(&c_dir)?;
    use std::io::{BufWriter, Write};
    let chunks = {
        let mut reference = c_dir.clone();
        reference.push("chunks.fa");
        let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
        for chunk in ds.selected_chunks.iter() {
            writeln!(wtr, ">{}\n{}", chunk.id, &chunk.seq)?;
        }
        wtr.flush()?;
        reference.into_os_string().into_string().unwrap()
    };
    use crate::minimap2;
    let threads = format!("{}", config.threads);
    let mut args = vec!["-t", &threads, "-c", "--eqx", "-P"];
    match ds.read_type {
        ReadType::CCS => args.extend(vec!["-H", "-k", "18"]),
        ReadType::CLR => args.extend(vec!["-H", "-k", "15"]),
        ReadType::ONT => args.extend(vec!["-k", "17"]),
        _ => {}
    };
    args.push("-X");
    let mm2 = minimap2::minimap2_args(&chunks, &chunks, &args);
    debug!("DETERMINE\tUnitOvlp\tRemoving\t{:?}", c_dir);
    std::fs::remove_dir_all(c_dir)?;
    Ok(mm2)
}

pub fn is_proper_overlap(paf: &bio_utils::paf::PAF) -> bool {
    const ALLOWED_END_GAP: usize = 25;
    let (qstart, qend, qlen) = (paf.qstart, paf.qend, paf.qlen);
    let (tstart, tend, tlen) = (paf.tstart, paf.tend, paf.tlen);
    let is_forward = paf.relstrand;
    // Q:-------->
    // T:    --------->
    let q_to_t_forward = (qlen - qend < ALLOWED_END_GAP) && tstart < ALLOWED_END_GAP && is_forward;
    // Q:<--------
    // T:    --------->
    let q_to_t_rev = qstart < ALLOWED_END_GAP && tstart < ALLOWED_END_GAP && !is_forward;
    // Q:     ------->
    // T: ------->
    let t_to_q_forward = (tlen - tend < ALLOWED_END_GAP) && qstart < ALLOWED_END_GAP && is_forward;
    // Q:    <--------
    // T: ------->
    let t_to_q_rev =
        (qlen - qend < ALLOWED_END_GAP) && (tlen - tend < ALLOWED_END_GAP) && !is_forward;
    q_to_t_forward || q_to_t_rev || t_to_q_forward || t_to_q_rev
}

fn remove_overlapping_chunks(
    ds: &mut DataSet,
    overlap_thr: f64,
    config: &DetermineUnitConfig,
) -> std::io::Result<()> {
    let chunk_len = ds.selected_chunks.len();
    // How long one overlap should be at least.
    let overlap_len = config.chunk_len / 2;
    let mm2 = mm2_chunk_overlap(ds, config)?;
    let alignments = std::str::from_utf8(&mm2).unwrap();
    let alignments = alignments
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(is_proper_overlap)
        .filter(|paf| {
            let identity = paf.matchnum as f64 / paf.blocklen as f64;
            overlap_thr < identity && overlap_len < paf.blocklen
        });
    let mut graph = vec![vec![]; chunk_len];
    for aln in alignments {
        let chunk1: usize = aln.tname.parse().unwrap();
        let chunk2: usize = aln.qname.parse().unwrap();
        if chunk1 != chunk2 {
            graph[chunk1].push(chunk2);
            graph[chunk2].push(chunk1);
        }
    }
    for edges in graph.iter_mut() {
        edges.sort_unstable();
        edges.dedup();
    }
    let to_be_removed = approx_vertex_cover(graph, ds.selected_chunks.len());
    ds.selected_chunks
        .retain(|chunk| !to_be_removed[chunk.id as usize]);
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if to_be_removed[node.chunk as usize] => read.remove(idx),
                Some(_) => idx += 1,
                None => return,
            }
        }
    });
    Ok(())
}

fn dump_histogram(ds: &DataSet) {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.chunk).or_default() += 1;
    }
    let mut counts: Vec<usize> = counts.iter().map(|x| *(x.1)).collect();
    counts.sort_unstable();
    let end_pos = counts.len().max(10) - 10;
    debug!("Hi-Freqs:{:?}", &counts[end_pos..]);
    let hist = histgram_viz::Histgram::new(&counts[..end_pos]);
    debug!("Histgrapm\n{}", hist.format(40, 20));
}

type NormedEdge = ((u64, bool), (u64, bool));
fn normalize_edge(w: &[Node]) -> (NormedEdge, bool) {
    let forward = ((w[0].chunk, w[0].is_forward), (w[1].chunk, w[1].is_forward));
    let reverse = (
        (w[1].chunk, !w[1].is_forward),
        (w[0].chunk, !w[0].is_forward),
    );
    if forward <= reverse {
        (forward, true)
    } else {
        (reverse, false)
    }
}

type FilledEdge = ((u64, bool), (u64, bool));
// Offset from the `from` position.
const SKIP_OFFSET: usize = 5;
type FilledEdges = HashMap<FilledEdge, Chunk>;
fn enumerate_filled_edges(
    ds: &DataSet,
    repetitive_kmers: &crate::repeat_masking::RepeatAnnot,
    config: &DetermineUnitConfig,
) -> HashMap<FilledEdge, Vec<u8>> {
    let mut edge_count: HashMap<_, Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter().filter(|r| !r.nodes.is_empty()) {
        assert_eq!(read.nodes.len(), read.edges.len() + 1);
        for (edge, w) in read.edges.iter().zip(read.nodes.windows(2)) {
            let label = edge.label();
            if config.chunk_len + SKIP_OFFSET < label.len() {
                let (edge, is_forward) = normalize_edge(w);
                let mut label: Vec<_> = if is_forward {
                    label
                        .iter()
                        .skip(SKIP_OFFSET)
                        .take(config.chunk_len)
                        .copied()
                        .collect()
                } else {
                    let start = label.len() - SKIP_OFFSET - config.chunk_len;
                    let end = label.len() - SKIP_OFFSET;
                    bio_utils::revcmp(&label[start..end])
                };
                label.iter_mut().for_each(u8::make_ascii_uppercase);
                edge_count.entry(edge).or_default().push(label);
            }
        }
    }
    // We fill chunks for each sparsed region.
    let count_thr = get_count_thr(ds, config);
    debug!("FillSparse\tEdge\tThreshold\t{}", count_thr);
    edge_count.retain(|_, seqs| {
        count_thr < seqs.len()
            && seqs
                .iter()
                .all(|s| repetitive_kmers.repetitiveness(s) < config.exclude_repeats)
    });
    take_consensus(&edge_count, &ds.read_type, config)
}

fn take_consensus<K: Hash + Clone + Eq + Sync + Send>(
    chunks: &HashMap<K, Vec<Vec<u8>>>,
    read_type: &ReadType,
    config: &DetermineUnitConfig,
) -> HashMap<K, Vec<u8>> {
    chunks
        .par_iter()
        .map(|(key, seqs)| {
            let radius = read_type.band_width(config.chunk_len);
            let draft = pick_median_length(seqs.as_slice());
            let consensus = kiley::bialignment::guided::polish_until_converge(&draft, seqs, radius);
            (key.clone(), consensus)
        })
        .collect()
}

fn pick_median_length(xs: &[Vec<u8>]) -> Vec<u8> {
    let mut lens: Vec<_> = xs.iter().map(|x| x.len()).collect();
    let len = xs.len() / 2;
    let (_, &mut median, _) = lens.select_nth_unstable(len);
    xs.iter().find(|x| x.len() == median).unwrap().clone()
}

fn get_count_thr(ds: &DataSet, _config: &DetermineUnitConfig) -> usize {
    let mut count: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *count.entry(node.chunk).or_default() += 1;
    }
    let mut count: Vec<_> = count.into_values().collect();
    let median = count.len() / 2;
    *count.select_nth_unstable(median).1 / 4
}
fn fill_edge(
    read: &mut EncodedRead,
    seq: &[u8],
    edge_chunks: &FilledEdges,
    readtype: ReadType,
    config: &DetermineUnitConfig,
) {
    let inserts = fill_sparse_edges_in_read(read, seq, edge_chunks, readtype, config);
    for (accum_inserts, (idx, node)) in inserts.into_iter().enumerate() {
        match idx + accum_inserts {
            pos if pos < read.nodes.len() => read.nodes.insert(idx + accum_inserts, node),
            _ => read.nodes.push(node),
        }
    }
    re_encode_read(read, seq);
}

fn fill_sparse_edges_in_read(
    read: &EncodedRead,
    seq: &[u8],
    edge_chunks: &FilledEdges,
    readtype: ReadType,
    _config: &DetermineUnitConfig,
) -> Vec<(usize, Node)> {
    let mut inserts = vec![];
    for (idx, w) in read.nodes.windows(2).enumerate() {
        let (edge, is_forward) = normalize_edge(w);
        let start = w[0].position_from_start + w[0].seq().len();
        let end = w[1].position_from_start;
        if end <= start {
            continue;
        }
        if let Some(chunk) = edge_chunks.get(&edge) {
            if let Some(node) = fill_gap(seq, start, end, is_forward, chunk, readtype) {
                inserts.push((idx + 1, node))
            }
        }
    }
    inserts
}

fn fill_gap(
    seq: &[u8],
    start: usize,
    end: usize,
    direction: bool,
    chunk: &Chunk,
    readtype: ReadType,
) -> Option<Node> {
    let chunk_len = chunk.seq().len();
    let (start, end) = match direction {
        true => (start, (start + chunk_len + SKIP_OFFSET).min(end)),
        false => (end.saturating_sub(chunk_len + SKIP_OFFSET), end),
    };
    let mut seq = match direction {
        true => seq[start..end].to_vec(),
        false => bio_utils::revcmp(&seq[start..end]),
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    let alignment = edlib_sys::align(chunk.seq(), &seq, mode, task);
    let (seq_start, seq_end) = alignment.location().unwrap();
    let seq_end = seq_end + 1;
    let position_from_start = match direction {
        true => start + seq_start,
        false => start + seq.len() - seq_end,
    };
    let seq = &seq[seq_start..seq_end];
    let band = readtype.band_width(chunk.seq().len());
    let edlib_to_op = {
        use kiley::Op::*;
        [Match, Del, Ins, Mismatch]
    };
    let alignment = alignment.operations().unwrap();
    let ops: Vec<_> = alignment.iter().map(|&x| edlib_to_op[x as usize]).collect();
    let (_, ops) =
        kiley::bialignment::guided::global_guided(chunk.seq(), seq, &ops, band, ALN_PARAMETER);
    let mat_num = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
    let identity = mat_num as f64 / ops.len() as f64;
    (1f64 - identity < readtype.sim_thr()).then(|| {
        let cigar = crate::misc::kiley_op_to_ops(&ops).0;
        let seq = seq.to_vec();
        Node::new(chunk.id, direction, seq, cigar, position_from_start, 2)
    })
}

fn re_encode_read(read: &mut EncodedRead, seq: &[u8]) {
    if !read.nodes.is_empty() {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use crate::encode::{
            nodes_to_encoded_read, remove_overlapping_encoding, remove_slippy_alignment,
        };
        nodes.sort_by_key(|n| n.chunk);
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        nodes = remove_overlapping_encoding(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
}

fn fill_sparse_region(
    ds: &mut DataSet,
    repetitive_kmers: &crate::repeat_masking::RepeatAnnot,
    config: &DetermineUnitConfig,
) -> usize {
    let max_idx: u64 = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    let edge_chunks: HashMap<_, _> = enumerate_filled_edges(ds, repetitive_kmers, config)
        .into_iter()
        .filter(|(_, seq)| repetitive_kmers.repetitiveness(seq) < config.exclude_repeats)
        .enumerate()
        .map(|(idx, (key, seq))| {
            let chunk = Chunk::new(max_idx + 1 + idx as u64, seq, 2);
            (key, chunk)
        })
        .collect();
    let rawseq: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let readtype = ds.read_type;
    let len = edge_chunks.len();
    debug!("FillSparse\tEdge\tCosed\t{len}");
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let rawseq = &rawseq[&read.id];
        fill_edge(read, rawseq, &edge_chunks, readtype, config);
    });
    debug!("FillSparse\tEdge\t{len}");
    ds.selected_chunks.extend(edge_chunks.into_values());
    len
}

type FilledTips = HashMap<(u64, bool), Chunk>;
fn fill_tips(
    ds: &mut DataSet,
    repetitive_kmers: &crate::repeat_masking::RepeatAnnot,
    config: &DetermineUnitConfig,
) -> usize {
    let max_idx = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    let tip_chunks: HashMap<_, _> = enumerate_filled_tips(ds, config)
        .into_iter()
        .filter(|(_, seq)| repetitive_kmers.repetitiveness(seq) < config.exclude_repeats)
        .enumerate()
        .map(|(idx, (key, seq))| (key, Chunk::new(max_idx + 1 + idx as u64, seq, 2)))
        .collect();
    let rawseq: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let readtype = ds.read_type;
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let rawseq = &rawseq[&read.id];
        fill_tip(read, rawseq, &tip_chunks, readtype, config);
    });
    let len = tip_chunks.len();
    debug!("FillSparse\tTip\t{len}");
    ds.selected_chunks.extend(tip_chunks.into_values());
    len
}

fn enumerate_filled_tips(
    ds: &DataSet,
    config: &DetermineUnitConfig,
) -> HashMap<(u64, bool), Vec<u8>> {
    let mut tail_counts: HashMap<_, Vec<_>> = HashMap::new();
    let take_len = config.chunk_len + SKIP_OFFSET;
    for read in ds.encoded_reads.iter() {
        if let Some(head) = read.nodes.first() {
            let tip = &read.leading_gap;
            if take_len < tip.len() {
                tail_counts
                    .entry((head.chunk, !head.is_forward))
                    .or_default()
                    .push((tip, false));
            }
        }
        if let Some(tail) = read.nodes.last() {
            let tip = &read.trailing_gap;
            if take_len < tip.len() {
                tail_counts
                    .entry((tail.chunk, tail.is_forward))
                    .or_default()
                    .push((tip, true));
            }
        }
    }
    // We fill chunks for each sparsed region.
    let count_thr = get_count_thr(ds, config);
    tail_counts.retain(|_, labels| count_thr.max(1) < labels.len());
    debug!("FillSparse\tTip\tThreshold\t{}", count_thr);
    let tail_counts: HashMap<_, _> = tail_counts
        .into_par_iter()
        .map(|(key, labels)| {
            let mut seqs: Vec<_> = labels
                .iter()
                .map(|(tip, is_forward)| match is_forward {
                    true => tip.as_slice()[SKIP_OFFSET..take_len].to_vec(),
                    false => {
                        let start = tip.len() - take_len;
                        let end = tip.len() - SKIP_OFFSET;
                        bio_utils::revcmp(&tip.as_slice()[start..end])
                    }
                })
                .collect();
            seqs.iter_mut()
                .for_each(|r| r.iter_mut().for_each(u8::make_ascii_uppercase));
            (key, seqs)
        })
        .collect();
    take_consensus(&tail_counts, &ds.read_type, config)
}

fn fill_tip(
    read: &mut EncodedRead,
    seq: &[u8],
    tip_chunks: &FilledTips,
    readtype: ReadType,
    _config: &DetermineUnitConfig,
) {
    let head_tip = read
        .nodes
        .first()
        .filter(|node| SKIP_OFFSET < node.position_from_start)
        .and_then(|node| {
            tip_chunks
                .get(&(node.chunk, !node.is_forward))
                .and_then(|chunk| {
                    let start = 0;
                    let end = node.position_from_start;
                    fill_gap(seq, start, end, false, chunk, readtype)
                })
        });
    let tail_tip = read
        .nodes
        .last()
        .filter(|n| n.position_from_start + n.seq().len() < seq.len().saturating_sub(SKIP_OFFSET))
        .and_then(|node| {
            tip_chunks
                .get(&(node.chunk, node.is_forward))
                .and_then(|chunk| {
                    let start = node.position_from_start + node.seq().len();
                    let end = seq.len();
                    fill_gap(seq, start, end, true, chunk, readtype)
                })
        });
    let is_updated = head_tip.is_some() | tail_tip.is_some();
    if let Some(head) = head_tip {
        read.nodes.reverse();
        read.nodes.push(head);
        read.nodes.reverse();
    }
    if let Some(tail) = tail_tip {
        read.nodes.push(tail)
    }
    if is_updated {
        re_encode_read(read, seq);
    }
}

// Or, something went wrong? Please check it out.
fn split_into<'a>(r: &'a RawRead, c: &DetermineUnitConfig) -> Vec<&'a [u8]> {
    let seq = r.seq();
    if seq.len() < c.margin * 2 {
        vec![]
    } else {
        let bound = seq.len() - c.margin;
        (0..)
            .map(|i| (c.margin + i * c.chunk_len, c.margin + (i + 1) * c.chunk_len))
            .take_while(|&(_, end)| end <= bound)
            .map(|(s, e)| &seq[s..e])
            .collect()
    }
}

fn filter_chunk_by_ovlp(ds: &mut DataSet, config: &DetermineUnitConfig) {
    let overlap_thr: usize = match ds.read_type {
        ReadType::CCS => config.chunk_len / 2,
        ReadType::CLR => config.chunk_len / 3,
        ReadType::ONT => config.chunk_len / 3,
        ReadType::None => config.chunk_len / 3,
    };
    let chunk_len = ds.selected_chunks.iter().map(|u| u.id).max().unwrap();
    let chunk_len = chunk_len as usize + 1;
    assert!(ds.selected_chunks.len() <= chunk_len);
    let mut graph = vec![vec![]; chunk_len];
    for read in ds.encoded_reads.iter() {
        for (i, node1) in read.nodes.iter().enumerate() {
            for node2 in read.nodes.iter().skip(i + 1) {
                let node_end = node1.position_from_start + node1.seq.as_slice().len();
                let mode_start = node2.position_from_start;
                let ovlp_len = node_end.max(mode_start) - mode_start;
                if overlap_thr < ovlp_len {
                    let (node1, node2) = (node1.chunk as usize, node2.chunk as usize);
                    if !graph[node1].contains(&node2) {
                        assert!(!graph[node2].contains(&node1));
                        graph[node1].push(node2);
                        graph[node2].push(node1);
                    }
                }
            }
        }
    }

    let to_be_removed = approx_vertex_cover(graph, chunk_len);
    ds.selected_chunks
        .retain(|chunk| !to_be_removed[chunk.id as usize]);
    debug!("UNITNUM\t{}\tVertexCovering", ds.selected_chunks.len());
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if to_be_removed[node.chunk as usize] => read.remove(idx),
                Some(_) => idx += 1,
                None => return,
            }
        }
    });
}

fn approx_vertex_cover(mut edges: Vec<Vec<usize>>, nodes: usize) -> Vec<bool> {
    let mut degrees: Vec<_> = edges.iter().map(|n| n.len()).collect();
    let mut to_be_removed = vec![false; nodes];
    loop {
        let (argmax, &max) = degrees.iter().enumerate().max_by_key(|x| x.1).unwrap();
        if max == 0 {
            break;
        }
        to_be_removed[argmax] = true;
        let removed: Vec<_> = edges[argmax].clone();
        edges[argmax].clear();
        for &to in removed.iter() {
            degrees[to] = degrees[to].saturating_sub(1);
            edges[to].retain(|&n| n != argmax);
        }
        degrees[argmax] = 0;
    }
    to_be_removed
}

fn error(node: &definitions::Node, ref_chunk: &Chunk) -> f64 {
    let (query, aln, refr) = node.recover(ref_chunk);
    let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
    let del = query.iter().filter(|&&x| x == b' ').count() as f64;
    let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
    let aln_len = aln.len() as f64;
    (mismat + del + ins) / aln_len
}

// Return the error rate of `quantile`-quantile.
pub fn calc_sim_thr(ds: &DataSet, quantile: f64) -> f64 {
    use rayon::prelude::*;
    let ref_chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let mut error_rates: Vec<f64> = ds
        .encoded_reads
        .par_iter()
        .flat_map(|r| {
            r.nodes
                .par_iter()
                .map(|n| error(n, ref_chunks[&n.chunk]))
                .collect::<Vec<_>>()
        })
        .collect();
    error_rates.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!(quantile <= 1f64);
    let idx = ((error_rates.len() as f64 * quantile).floor() as usize).min(error_rates.len() - 1);
    error_rates[idx]
}
