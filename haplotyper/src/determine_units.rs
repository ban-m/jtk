use crate::ALN_PARAMETER;

use super::encode::Encode;
use super::polish_units::PolishUnit;
use super::polish_units::PolishUnitConfig;
use definitions::*;
use rand::prelude::*;
use rand_xoshiro::Xoroshiro128Plus;
use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct UnitConfig {
    pub chunk_len: usize,
    pub skip_len: usize,
    pub unit_num: usize,
    pub margin: usize,
    pub threads: usize,
    pub min_cluster: usize,
    pub exclude_repeats: f64,
    pub upper_count: usize,
    pub lower_count: usize,
}

const FIRST_RELAX: f64 = 2f64;

const TAKE_THR: f64 = 0.999;
impl UnitConfig {
    #[allow(clippy::too_many_arguments)]
    pub fn new_ccs(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
        exclude_repeats: f64,
        upper_count: usize,
        lower_count: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            unit_num,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
    #[allow(clippy::too_many_arguments)]
    pub fn new_clr(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
        exclude_repeats: f64,
        upper_count: usize,
        lower_count: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            unit_num,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
    #[allow(clippy::too_many_arguments)]
    pub fn new_ont(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
        exclude_repeats: f64,
        upper_count: usize,
        lower_count: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            unit_num,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
}

pub trait DetermineUnit {
    fn select_chunks(&mut self, config: &UnitConfig);
}

impl DetermineUnit for definitions::DataSet {
    // TODO: We can make this process much faster, by just skipping the needless re-encoding.
    // TODO: Parametrize the number of reads used in consensus generation.
    // TOOD: Maybe we can remove some low-quality reads by
    fn select_chunks(&mut self, config: &UnitConfig) {
        let filter_size = match self.read_type {
            ReadType::CCS => 2,
            ReadType::None | ReadType::CLR => 5,
            ReadType::ONT => 3,
        };
        debug!("Select Unit: Configuration:{:?}", config);
        let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(342903);
        self.selected_chunks = {
            let mut config = config.clone();
            config.chunk_len = 12 * config.chunk_len / 10;
            pick_random(&self.raw_reads, &config, &mut rng)
        };
        debug!("UNITNUM\t{}\tPICKED", self.selected_chunks.len());
        remove_overlapping_units_dev(self, config).unwrap();
        compaction_units(self);
        // 1st polishing.
        use crate::stats::Stats;
        let mut sim_thr = self.read_type.sim_thr();
        {
            debug!("UNITNUM\t{}\tREMOVED", self.selected_chunks.len());
            self.encode(config.threads, FIRST_RELAX * sim_thr);
            debug!("ERRORRATE\t{}", self.error_rate());
            remove_frequent_units(self, config.upper_count);
            dump_histogram(self);
            let polish_config = PolishUnitConfig::new(self.read_type, filter_size, 30);
            self.consensus_unit(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED\t1", self.selected_chunks.len());
        }
        // 2nd polishing.
        {
            self.encode(config.threads, sim_thr);
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            for _ in 0..4 {
                fill_sparse_region_dev(self, config);
                fill_tips_dev(self, config);
                crate::encode::deletion_fill::correct_unit_deletion(self, sim_thr);
            }
            compaction_units(self);
            remove_overlapping_units_dev(self, config).unwrap();
            remove_frequent_units(self, config.upper_count);
            let polish_config = PolishUnitConfig::new(self.read_type, filter_size, 30);
            dump_histogram(self);
            self.polish_unit(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED\t2", self.selected_chunks.len());
        };
        // Final polish
        {
            debug!("UNITNUM\t{}\tRAWUNIT", self.selected_chunks.len());
            self.encode(config.threads, sim_thr);
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            remove_frequent_units(self, config.upper_count);
            filter_unit_by_ovlp(self, config);
            compaction_units(self);
            debug!("UNITNUM\t{}\tFILTERED", self.selected_chunks.len());
            remove_frequent_units(self, config.upper_count);
            sim_thr = calc_sim_thr(self, TAKE_THR).max(self.read_type.sim_thr());
            debug!("ERRORRATE\t{}\t{}", self.error_rate(), sim_thr);
            let polish_config = PolishUnitConfig::new(self.read_type, 2 * filter_size, 100);
            dump_histogram(self);
            self.polish_unit(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED\t3", self.selected_chunks.len());
        }
        // Tune the length, removing high-frequent units...
        self.selected_chunks
            .retain(|unit| config.chunk_len < unit.seq.len());
        self.selected_chunks
            .iter_mut()
            .for_each(|unit| unit.seq.truncate(config.chunk_len));
        {
            self.encode(config.threads, sim_thr);
            debug!("ERRORRATE\t{}", self.error_rate());
            remove_frequent_units(self, config.upper_count);
            dump_histogram(self);
        }
        compaction_units(self);
    }
}

// fn re_encode(read: &mut EncodedRead, convert_table: &HashMap<u64, u64>, seq: &[u8]) {
//     let mut nodes = Vec::with_capacity(read.nodes.len());
//     nodes.append(&mut read.nodes);
//     nodes
//         .iter_mut()
//         .for_each(|n| n.unit = convert_table[&n.unit]);
//     *read = crate::encode::nodes_to_encoded_read(read.id, nodes, seq).unwrap();
// }

fn remove_frequent_units(ds: &mut DataSet, upper_count: usize) {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    counts.retain(|_, occ| *occ > upper_count);
    ds.selected_chunks.retain(|u| !counts.contains_key(&u.id));
    for read in ds.encoded_reads.iter_mut() {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if counts.contains_key(&node.unit) => read.remove(idx),
                Some(_) => idx += 1,
                None => break,
            }
        }
    }
}

// Make the id consective
fn compaction_units(ds: &mut DataSet) {
    let mut mapping: HashMap<u64, u64> = HashMap::new();
    for (idx, unit) in ds.selected_chunks.iter_mut().enumerate() {
        mapping.insert(unit.id, idx as u64);
        unit.id = idx as u64;
    }
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        read.nodes
            .iter_mut()
            .for_each(|n| n.unit = mapping[&n.unit]);
        read.edges.iter_mut().for_each(|e| {
            e.from = mapping[&e.from];
            e.to = mapping[&e.to];
        });
    });
}

use rand::Rng;
fn pick_random<R: Rng>(reads: &[RawRead], config: &UnitConfig, rng: &mut R) -> Vec<Unit> {
    use rand::prelude::*;
    let subseqs: Vec<_> = reads
        .iter()
        .flat_map(|r| split_into(r, config))
        .filter(|u| !is_repetitive(u, config))
        .collect();
    subseqs
        .choose_multiple(rng, config.unit_num)
        .enumerate()
        .map(|(idx, seq)| {
            let seq = String::from_utf8_lossy(seq).to_string();
            Unit::new(idx as u64, seq, config.min_cluster)
        })
        .collect()
}

fn mm2_unit_overlap(ds: &DataSet, config: &UnitConfig) -> std::io::Result<Vec<u8>> {
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
        for unit in ds.selected_chunks.iter() {
            writeln!(&mut wtr, ">{}\n{}", unit.id, &unit.seq)?;
        }
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

fn remove_overlapping_units_dev(ds: &mut DataSet, config: &UnitConfig) -> std::io::Result<()> {
    let unit_len = ds.selected_chunks.len();
    // How long one overlap should be at least.
    let overlap_len = config.chunk_len / 2;
    // This is the percent identy.
    let overlap_thr: f64 = match ds.read_type {
        ReadType::CCS => 0.95,
        ReadType::CLR => 0.75,
        ReadType::ONT => 0.85,
        ReadType::None => 0.85,
    };
    let mm2 = mm2_unit_overlap(ds, config)?;
    let alignments = String::from_utf8_lossy(&mm2);
    let alignments = alignments
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(is_proper_overlap)
        .filter(|paf| {
            let identity = paf.matchnum as f64 / paf.blocklen as f64;
            overlap_thr < identity && overlap_len < paf.blocklen
        });
    // TODO: FIXME.
    let mut edges: Vec<_> = (0..unit_len).map(|i| vec![0; i]).collect();
    for aln in alignments {
        let node_unit: usize = aln.tname.parse().unwrap();
        let mode_unit: usize = aln.qname.parse().unwrap();
        let (i, j) = (node_unit.max(mode_unit), node_unit.min(mode_unit));
        if i != j {
            edges[i as usize][j as usize] += 1;
        }
    }
    let edges: Vec<Vec<_>> = edges
        .iter()
        .map(|xs| xs.iter().map(|&x| MIN_OCC < x).collect())
        .collect();
    let to_be_removed = approx_vertex_cover(edges, ds.selected_chunks.len());
    ds.selected_chunks
        .retain(|unit| !to_be_removed[unit.id as usize]);
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if to_be_removed[node.unit as usize] => read.remove(idx),
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
        *counts.entry(node.unit).or_default() += 1;
    }
    let mut counts: Vec<usize> = counts.iter().map(|x| *(x.1)).collect();
    counts.sort_unstable();
    let end_pos = counts.len().max(10) - 10;
    debug!("Hi-Freqs:{:?}", &counts[end_pos..]);
    let hist = histgram_viz::Histgram::new(&counts[..end_pos]);
    debug!("DUMP\t{:?}", &counts[..end_pos]);
    debug!("Histgrapm\n{}", hist.format(40, 20));
}

const MIN_OCC: usize = 5;
#[allow(dead_code)]
fn remove_overlapping_units(ds: &DataSet, thr: usize) -> std::io::Result<Vec<Unit>> {
    const ALLOWED_END_GAP: usize = 50;
    let mm2 = crate::encode::mm2_alignment(ds, thr)?;
    let alignments: Vec<_> = String::from_utf8_lossy(&mm2)
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(|aln| aln.tstart < ALLOWED_END_GAP && aln.tlen - aln.tend < ALLOWED_END_GAP)
        .collect();
    let mut buckets: HashMap<_, Vec<_>> = HashMap::new();
    for aln in alignments.iter() {
        buckets.entry(aln.qname.clone()).or_default().push(aln);
    }
    buckets
        .values_mut()
        .for_each(|xs| xs.sort_by_key(|aln| aln.qstart));
    let unit_len = ds.selected_chunks.len();
    let mut edges: Vec<_> = (0..unit_len).map(|i| vec![0; i]).collect();
    let overlap_thr: usize = match ds.read_type {
        ReadType::CCS => 150,
        ReadType::CLR => 500,
        ReadType::ONT => 400,
        ReadType::None => 500,
    };

    for alns in buckets.values() {
        for (i, node) in alns.iter().enumerate() {
            for mode in alns.iter().skip(i + 1) {
                let node_unit: usize = node.tname.parse().unwrap();
                let mode_unit: usize = mode.tname.parse().unwrap();
                let ovlp_len = node.qend.saturating_sub(mode.qstart);
                if 2 * overlap_thr < ovlp_len {
                    let (i, j) = (node_unit.max(mode_unit), node_unit.min(mode_unit));
                    if i != j {
                        edges[i as usize][j as usize] += 1;
                    }
                }
            }
        }
    }
    let edges: Vec<Vec<_>> = edges
        .iter()
        .map(|xs| xs.iter().map(|&x| MIN_OCC < x).collect())
        .collect();
    let to_be_removed = approx_vertex_cover(edges, ds.selected_chunks.len());
    let mut chunks = ds.selected_chunks.clone();
    chunks.retain(|unit| !to_be_removed[unit.id as usize]);
    for (idx, unit) in chunks.iter_mut().enumerate() {
        unit.id = idx as u64;
    }
    Ok(chunks)
}

type NormedEdge = ((u64, bool), (u64, bool));
fn normalize_edge(w: &[Node]) -> (NormedEdge, bool) {
    let forward = ((w[0].unit, w[0].is_forward), (w[1].unit, w[1].is_forward));
    let reverse = ((w[1].unit, !w[1].is_forward), (w[0].unit, !w[0].is_forward));
    if forward <= reverse {
        (forward, true)
    } else {
        (reverse, false)
    }
}

type FilledEdge = ((u64, bool), (u64, bool));
// Offset from the `from` position.
const SKIP_OFFSET: usize = 5;
type FilledEdges = HashMap<FilledEdge, Unit>;
fn enumerate_filled_edges(ds: &DataSet, config: &UnitConfig) -> FilledEdges {
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
    // We fill units for each sparsed region.
    let count_thr = get_count_thr(ds, config);
    debug!("FillSparse\tThreshold\t{}", count_thr);
    let mut edge_units: HashMap<_, _> = edge_count
        .into_par_iter()
        .filter(|x| count_thr < x.1.len())
        .map(|(key, ref seqs)| {
            let radius = ds.read_type.band_width(config.chunk_len);
            let cons = kiley::ternary_consensus_by_chunk(seqs, radius);
            let consensus = kiley::bialignment::guided::polish_until_converge(&cons, &seqs, radius);
            let seq = String::from_utf8(consensus).unwrap();
            let unit = Unit::new(0, seq, 2);
            (key, unit)
        })
        .collect();
    let max_idx: u64 = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    edge_units.values_mut().enumerate().for_each(|(idx, unit)| {
        unit.id = max_idx + 1 + idx as u64;
    });
    edge_units
}

fn get_count_thr(ds: &DataSet, _config: &UnitConfig) -> usize {
    let mut count: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *count.entry(node.unit).or_default() += 1;
    }
    let mut count: Vec<_> = count.into_values().collect();
    let median = count.len() / 2;
    *count.select_nth_unstable(median).1 / 4
}
fn fill_edge(
    read: &mut EncodedRead,
    seq: &[u8],
    edge_units: &FilledEdges,
    readtype: ReadType,
    config: &UnitConfig,
) {
    let inserts = fill_sparse_edges_in_read(read, seq, edge_units, readtype, config);
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
    edge_units: &FilledEdges,
    readtype: ReadType,
    _config: &UnitConfig,
) -> Vec<(usize, Node)> {
    let mut inserts = vec![];
    for (idx, w) in read.nodes.windows(2).enumerate() {
        let (edge, is_forward) = normalize_edge(w);
        let start = w[0].position_from_start + w[0].seq().len();
        let end = w[1].position_from_start;
        if let Some(unit) = edge_units.get(&edge) {
            if let Some(node) = fill_gap(seq, start, end, is_forward, unit, readtype) {
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
    unit: &Unit,
    readtype: ReadType,
) -> Option<Node> {
    // let (old_start, old_end) = (start, end);
    let unit_len = unit.seq().len();
    let (start, end) = match direction {
        true => (start, (start + unit_len + SKIP_OFFSET).min(end)),
        false => (end.saturating_sub(unit_len + SKIP_OFFSET), end),
    };
    let mut seq = match direction {
        true => seq[start..end].to_vec(),
        false => bio_utils::revcmp(&seq[start..end]),
    };
    seq.iter_mut().for_each(u8::make_ascii_uppercase);
    let mode = edlib_sys::AlignMode::Infix;
    let task = edlib_sys::AlignTask::Alignment;
    let alignment = edlib_sys::edlib_align(unit.seq(), &seq, mode, task);
    let (seq_start, seq_end) = alignment.locations.unwrap()[0];
    let seq_end = seq_end + 1;
    let position_from_start = match direction {
        true => start + seq_start,
        false => start + seq.len() - seq_end,
    };
    let seq = &seq[seq_start..seq_end];
    let band = readtype.band_width(unit.seq().len());
    let edlib_to_op = {
        use kiley::Op::*;
        [Match, Del, Ins, Mismatch]
    };
    let alignment = alignment.operations.unwrap();
    let ops: Vec<_> = alignment.iter().map(|&x| edlib_to_op[x as usize]).collect();
    let (_, ops) =
        kiley::bialignment::guided::global_guided(unit.seq(), seq, &ops, band, ALN_PARAMETER);
    let mat_num = ops.iter().filter(|&&op| op == kiley::Op::Match).count();
    let identity = mat_num as f64 / ops.len() as f64;
    (1f64 - identity < readtype.sim_thr()).then(|| {
        let cigar = crate::encode::compress_kiley_ops(&ops);
        Node::new(unit.id, direction, seq, cigar, position_from_start, 2)
    })
}

fn re_encode_read(read: &mut EncodedRead, seq: &[u8]) {
    if !read.nodes.is_empty() {
        let mut nodes = vec![];
        nodes.append(&mut read.nodes);
        use crate::encode::{nodes_to_encoded_read, remove_slippy_alignment};
        nodes.sort_by_key(|n| n.unit);
        nodes = remove_slippy_alignment(nodes);
        nodes.sort_by_key(|n| n.position_from_start);
        nodes = remove_slippy_alignment(nodes);
        *read = nodes_to_encoded_read(read.id, nodes, seq).unwrap();
    }
}

fn fill_sparse_region_dev(ds: &mut DataSet, config: &UnitConfig) {
    let edge_units = enumerate_filled_edges(ds, config);
    let rawseq: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let readtype = ds.read_type;
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let rawseq = &rawseq[&read.id];
        fill_edge(read, rawseq, &edge_units, readtype, config);
    });
    debug!("FillSparse\t{}", edge_units.len());
    ds.selected_chunks.extend(edge_units.into_values());
}

#[allow(dead_code)]
fn fill_sparse_region(ds: &mut DataSet, config: &UnitConfig) {
    let mut edge_count: HashMap<_, Vec<_>> = HashMap::new();
    for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
        let key = (edge.from.min(edge.to), edge.from.max(edge.to));
        edge_count.entry(key).or_default().push(edge);
    }
    // We fill units for each sparsed region.
    let count_thr = {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *count.entry(node.unit).or_default() += 1;
        }
        let mut count: Vec<_> = count.into_values().collect();
        let median = count.len() / 2;
        *count.select_nth_unstable(median).1 / 4
    };
    let mut picked_units = vec![];
    for (_, edges) in edge_count {
        let total_len: i64 = edges
            .iter()
            .map(|edge| match edge.label.is_empty() {
                true => edge.offset,
                false => edge.label().len() as i64,
            })
            .sum();
        let mean_len = total_len / edges.len() as i64;
        if (config.chunk_len as i64) < mean_len && count_thr < edges.len() {
            let edge = edges.iter().max_by_key(|e| e.label().len()).unwrap();
            let seq = edge.label();
            let new_units = (0..)
                .map(|i| (config.chunk_len * i, (i + 1) * config.chunk_len))
                .take_while(|&(_, y)| y < seq.len())
                .map(|(s, t)| &seq[s..t])
                .filter(|u| !is_repetitive(u, config));
            picked_units.extend(new_units);
        }
    }
    debug!("FillSparse\t{}", picked_units.len());
    let last_unit = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    ds.selected_chunks
        .extend(picked_units.iter().enumerate().map(|(i, seq)| {
            let id = i as u64 + last_unit + 1;
            let seq = String::from_utf8_lossy(seq).to_string();
            Unit::new(id, seq, config.min_cluster)
        }));
}

type FilledTips = HashMap<(u64, bool), Unit>;
fn fill_tips_dev(ds: &mut DataSet, config: &UnitConfig) {
    let tip_units = enumerate_filled_tips(ds, config);
    let rawseq: HashMap<u64, _> = ds.raw_reads.iter().map(|r| (r.id, r.seq())).collect();
    let readtype = ds.read_type;
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let rawseq = &rawseq[&read.id];
        fill_tip(read, rawseq, &tip_units, readtype, config);
    });
    debug!("FillTip\t{}", tip_units.len());
    ds.selected_chunks.extend(tip_units.into_values());
}

fn enumerate_filled_tips(ds: &DataSet, config: &UnitConfig) -> FilledTips {
    let mut tail_counts: HashMap<_, Vec<_>> = HashMap::new();
    let take_len = config.chunk_len + SKIP_OFFSET;
    for read in ds.encoded_reads.iter() {
        if let Some(head) = read.nodes.first() {
            let tip = &read.leading_gap;
            if take_len < tip.len() {
                tail_counts
                    .entry((head.unit, !head.is_forward))
                    .or_default()
                    .push((tip, false));
            }
        }
        if let Some(tail) = read.nodes.last() {
            let tip = &read.trailing_gap;
            if take_len < tip.len() {
                tail_counts
                    .entry((tail.unit, tail.is_forward))
                    .or_default()
                    .push((tip, true));
            }
        }
    }
    // We fill units for each sparsed region.
    let count_thr = get_count_thr(&ds, config);
    debug!("FillTip\tThreshold\t{}", count_thr);
    let mut tip_units: HashMap<_, _> = tail_counts
        .par_iter()
        .filter(|(_, labels)| count_thr.max(1) < labels.len())
        .map(|(key, labels)| {
            let mut seqs: Vec<_> = labels
                .iter()
                .map(|(tip, is_forward)| match is_forward {
                    true => tip[SKIP_OFFSET..take_len].to_vec(),
                    false => {
                        let start = tip.len() - take_len;
                        let end = tip.len() - SKIP_OFFSET;
                        bio_utils::revcmp(&tip[start..end])
                    }
                })
                .collect();
            seqs.iter_mut()
                .for_each(|r| r.iter_mut().for_each(u8::make_ascii_uppercase));
            let radius = ds.read_type.band_width(config.chunk_len);
            let cons = kiley::ternary_consensus_by_chunk(&seqs, radius);
            let consensus = kiley::bialignment::guided::polish_until_converge(&cons, &seqs, radius);
            let seq = String::from_utf8(consensus).unwrap();
            (*key, Unit::new(0, seq, 2))
        })
        .collect();
    let max_idx = ds.selected_chunks.iter().map(|c| c.id).max().unwrap();
    tip_units.values_mut().enumerate().for_each(|(idx, unit)| {
        unit.id = max_idx + 1 + idx as u64;
    });
    tip_units
}

fn fill_tip(
    read: &mut EncodedRead,
    seq: &[u8],
    tip_units: &FilledTips,
    readtype: ReadType,
    _config: &UnitConfig,
) {
    let head_tip = read
        .nodes
        .first()
        .filter(|node| SKIP_OFFSET < node.position_from_start)
        .and_then(|node| {
            tip_units
                .get(&(node.unit, !node.is_forward))
                .and_then(|unit| {
                    let start = 0;
                    let end = node.position_from_start;
                    fill_gap(seq, start, end, false, unit, readtype)
                })
        });
    let tail_tip = read
        .nodes
        .last()
        .filter(|n| n.position_from_start + n.seq().len() < seq.len().saturating_sub(SKIP_OFFSET))
        .and_then(|node| {
            tip_units
                .get(&(node.unit, node.is_forward))
                .and_then(|unit| {
                    let start = node.position_from_start + node.seq().len();
                    let end = seq.len();
                    fill_gap(seq, start, end, true, unit, readtype)
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

#[allow(dead_code)]
fn fill_tail_end(ds: &mut DataSet, config: &UnitConfig) {
    let mut tail_counts: HashMap<_, Vec<_>> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        if let Some(head) = read.nodes.first() {
            tail_counts
                .entry((head.unit, head.is_forward))
                .or_default()
                .push(&read.leading_gap);
        }
        if let Some(tail) = read.nodes.last() {
            tail_counts
                .entry((tail.unit, !tail.is_forward))
                .or_default()
                .push(&read.trailing_gap);
        }
    }
    tail_counts.values_mut().for_each(|tips| {
        tips.retain(|label| config.chunk_len < label.len());
    });
    // We fill units for each sparsed region.
    let count_thr = {
        let mut count: HashMap<_, usize> = HashMap::new();
        for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
            *count.entry(node.unit).or_default() += 1;
        }
        let mut count: Vec<_> = count.into_values().collect();
        let median = count.len() / 2;
        *count.select_nth_unstable(median).1 / 8
    };
    let mut picked_units = vec![];
    for ((unit, direction), labels) in tail_counts
        .into_iter()
        .filter(|(_, labels)| labels.len() > count_thr.max(1))
    {
        let seq = labels.iter().max_by_key(|xs| xs.len()).unwrap();
        trace!("FillSparse\t{}\t{}\t{}", unit, direction, labels.len());
        let new_units = seq
            .chunks_exact(config.chunk_len)
            .filter(|u| !is_repetitive(u, config));
        picked_units.extend(new_units);
    }
    debug!("FillSparse\t{}\t{}", picked_units.len(), count_thr);
    let last_unit = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    ds.selected_chunks
        .extend(picked_units.iter().enumerate().map(|(i, seq)| {
            let id = i as u64 + last_unit + 1;
            let seq = String::from_utf8_lossy(seq).to_string();
            Unit::new(id, seq, config.min_cluster)
        }));
}

fn is_repetitive(unit: &[u8], config: &UnitConfig) -> bool {
    let tot = unit.len();
    let lowercase = unit.iter().filter(|c| c.is_ascii_lowercase()).count();
    lowercase as f64 / tot as f64 > config.exclude_repeats
}

fn split_into<'a>(r: &'a RawRead, c: &UnitConfig) -> Vec<&'a [u8]> {
    let seq = r.seq();
    if seq.len() < c.margin * 2 {
        vec![]
    } else {
        let end = seq.len() - c.margin;
        let stride = c.chunk_len + c.skip_len;
        (0..)
            .map(|i| (stride * i, stride * i + c.chunk_len))
            .map(|(x, y)| (x + c.margin, y + c.margin))
            .take_while(|&(_, y)| y < end)
            .map(|(s, t)| &seq[s..t])
            .collect()
    }
}

fn filter_unit_by_ovlp(ds: &mut DataSet, config: &UnitConfig) {
    let overlap_thr: usize = match ds.read_type {
        ReadType::CCS => config.chunk_len / 2,
        ReadType::CLR => config.chunk_len / 4,
        ReadType::ONT => config.chunk_len / 3,
        ReadType::None => config.chunk_len / 4,
    };
    let unit_len = ds.selected_chunks.iter().map(|u| u.id).max().unwrap();
    let unit_len = unit_len as usize + 1;
    assert!(ds.selected_chunks.len() <= unit_len);
    // TODO: FIXME
    // Maybe it became infeasible due to O(N^2) allcation would be occured.
    let mut edges: Vec<_> = (0..unit_len).map(|i| vec![false; i]).collect();
    for read in ds.encoded_reads.iter() {
        for (i, node) in read.nodes.iter().enumerate() {
            for mode in read.nodes.iter().skip(i + 1) {
                let node_end = node.position_from_start + node.seq.as_bytes().len();
                let mode_start = mode.position_from_start;
                let ovlp_len = node_end.max(mode_start) - mode_start;
                if overlap_thr < ovlp_len {
                    let (i, j) = (node.unit.max(mode.unit), node.unit.min(mode.unit));
                    if i != j {
                        edges[i as usize][j as usize] = true;
                    }
                }
            }
        }
    }
    let to_be_removed = approx_vertex_cover(edges, unit_len);
    let remaining = to_be_removed.iter().filter(|&&b| !b).count();
    debug!("UNITNUM\t{}\tVertexCovering", remaining);
    // let mut count: HashMap<_, i32> = HashMap::new();
    // for read in ds.encoded_reads.iter() {
    //     for node in read.nodes.iter() {
    //         *count.entry(node.unit).or_default() += 1;
    //     }
    // }
    // let median = {
    //     let mut count: Vec<_> = count.values().copied().collect();
    //     let median = count.len() / 2;
    //     *count.select_nth_unstable(median).1 as i32
    // };
    // let (lower, upper) = match ds.read_type {
    //     ReadType::CCS => ((median / 10).max(2), 10 * median),
    //     ReadType::CLR => ((median / 10).max(4), 5 * median),
    //     ReadType::ONT => ((median / 10).max(4), 5 * median),
    //     ReadType::None => (median / 10, 200),
    // };
    ds.selected_chunks.retain(|unit| {
        // let coverage = match count.get(&unit.id) {
        //     Some(res) => *res,
        //     None => return false,
        // };
        // !to_be_removed[unit.id as usize] && lower < coverage && coverage < upper
        !to_be_removed[unit.id as usize]
    });
    ds.encoded_reads.par_iter_mut().for_each(|read| {
        let mut idx = 0;
        loop {
            match read.nodes.get(idx) {
                Some(node) if to_be_removed[node.unit as usize] => read.remove(idx),
                Some(_) => idx += 1,
                None => return,
            }
        }
    });
    // ds.selected_chunks.iter_mut().for_each(|unit| {
    //     unit.id = idx;
    //     idx += 1;
    // });
    // ds.encoded_reads.clear();
}

fn approx_vertex_cover(mut edges: Vec<Vec<bool>>, nodes: usize) -> Vec<bool> {
    let mut degrees = vec![0; nodes];
    for (i, es) in edges.iter().enumerate() {
        for (j, _) in es.iter().enumerate().filter(|&(_, &b)| b) {
            degrees[i] += 1;
            degrees[j] += 1;
        }
    }
    let mut to_be_removed = vec![false; nodes];
    loop {
        let (argmax, &max) = degrees.iter().enumerate().max_by_key(|x| x.1).unwrap();
        if max == 0 {
            break;
        }
        to_be_removed[argmax] = true;
        edges[argmax].iter_mut().enumerate().for_each(|(i, b)| {
            degrees[i] -= *b as usize;
            *b = false;
        });
        degrees[argmax] = 0;
        if argmax < nodes {
            edges
                .iter_mut()
                .enumerate()
                .skip(argmax + 1)
                .for_each(|(i, es)| {
                    degrees[i] -= es[argmax] as usize;
                    es[argmax] = false;
                });
        }
    }
    to_be_removed
}

fn error(node: &definitions::Node, ref_unit: &Unit) -> f64 {
    let (query, aln, refr) = node.recover(ref_unit);
    let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
    let del = query.iter().filter(|&&x| x == b' ').count() as f64;
    let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
    let aln_len = aln.len() as f64;
    (mismat + del + ins) / aln_len
}

// Return the error rate of `quantile`-quantile.
pub fn calc_sim_thr(ds: &DataSet, quantile: f64) -> f64 {
    use rayon::prelude::*;
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let mut error_rates: Vec<f64> = ds
        .encoded_reads
        .par_iter()
        .flat_map(|r| {
            r.nodes
                .par_iter()
                .map(|n| error(n, ref_units[&n.unit]))
                .collect::<Vec<_>>()
        })
        .collect();
    error_rates.sort_by(|x, y| x.partial_cmp(y).unwrap());
    assert!(quantile <= 1f64);
    let idx = ((error_rates.len() as f64 * quantile).floor() as usize).max(error_rates.len() - 1);
    error_rates[idx]
}
