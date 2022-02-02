use super::encode::Encode;
use super::polish_units::PolishUnit;
use super::polish_units::PolishUnitConfig;
// use super::Encode;
use definitions::*;
use rand::prelude::*;
use rand_xoshiro::Xoroshiro128Plus;
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

// const DEFAULT_RNG: usize = 200;
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
        //self.selected_chunks = remove_overlapping_units(self, config.threads).unwrap();
        remove_overlapping_units_dev(self, config).unwrap();
        // 1st polishing.
        debug!("UNITNUM\t{}\tREMOVED", self.selected_chunks.len());
        self.encode(config.threads, self.read_type.sim_thr());
        remove_frequent_units(self, config.upper_count);
        dump_histogram(self);
        let polish_config = PolishUnitConfig::new(self.read_type, filter_size, 30);
        self.consensus_unit(&polish_config);
        debug!("UNITNUM\t{}\tPOLISHED\t1", self.selected_chunks.len());
        self.encode(config.threads, self.read_type.sim_thr());
        fill_sparse_region(self, config);
        fill_tail_end(self, config);
        // 2nd polishing.
        self.encode(config.threads, self.read_type.sim_thr());
        remove_frequent_units(self, config.upper_count);
        let polish_config = PolishUnitConfig::new(self.read_type, filter_size, 30);
        dump_histogram(self);
        self.polish_unit(&polish_config);
        debug!("UNITNUM\t{}\tPOLISHED\t2", self.selected_chunks.len());
        debug!("UNITNUM\t{}\tRAWUNIT", self.selected_chunks.len());
        self.encode(config.threads, self.read_type.sim_thr());
        remove_frequent_units(self, config.upper_count);
        filter_unit_by_ovlp(self, config);
        debug!("UNITNUM\t{}\tFILTERED", self.selected_chunks.len());
        self.encode(config.threads, self.read_type.sim_thr());
        // Final polish
        let polish_config = PolishUnitConfig::new(self.read_type, 2 * filter_size, 100);
        dump_histogram(self);
        self.polish_unit(&polish_config);
        self.selected_chunks
            .retain(|unit| config.chunk_len < unit.seq.len());
        for (idx, unit) in self.selected_chunks.iter_mut().enumerate() {
            unit.seq.truncate(config.chunk_len);
            unit.id = idx as u64;
        }
        self.encode(config.threads, self.read_type.sim_thr());
    }
}

fn remove_frequent_units(ds: &mut DataSet, upper_count: usize) {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    counts.retain(|_, occ| *occ > upper_count);
    ds.selected_chunks.retain(|u| !counts.contains_key(&u.id));
    for read in ds.encoded_reads.iter_mut() {
        'outer: loop {
            let len = read.nodes.len();
            for i in 0..len {
                if counts.contains_key(&read.nodes[i].unit) {
                    read.remove(i);
                    continue 'outer;
                }
            }
            break;
        }
    }
}

use rand::Rng;
fn pick_random<R: Rng>(reads: &[RawRead], config: &UnitConfig, rng: &mut R) -> Vec<Unit> {
    use rand::prelude::*;
    let subseqs: Vec<_> = reads
        .iter()
        .flat_map(|r| split_into(r, &config))
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
        let reference = reference.into_os_string().into_string().unwrap();
        reference
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
    let mut edges: Vec<_> = (0..unit_len).map(|i| vec![0; i]).collect();
    // How long one overlap should be at least.
    let overlap_len = config.chunk_len * 3 / 4;
    // This is the percent identy.
    let overlap_thr: f64 = match ds.read_type {
        ReadType::CCS => 0.05,
        ReadType::CLR => 0.35,
        ReadType::ONT => 0.15,
        ReadType::None => 0.15,
    };
    let mm2 = mm2_unit_overlap(ds, config)?;
    let alignments = String::from_utf8_lossy(&mm2);
    let alignments = alignments
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .filter(|paf| is_proper_overlap(paf))
        .filter(|paf| {
            let identity = paf.matchnum as f64 / paf.blocklen as f64;
            overlap_thr < identity && overlap_len < paf.blocklen
        });
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
    for (idx, unit) in ds.selected_chunks.iter_mut().enumerate() {
        unit.id = idx as u64;
    }
    Ok(())
}

fn dump_histogram(ds: &DataSet) {
    let mut counts: HashMap<_, usize> = HashMap::new();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        *counts.entry(node.unit).or_default() += 1;
    }
    let mut counts: Vec<usize> = counts.iter().map(|x| *(x.1)).collect();
    counts.sort_unstable();
    let hist = histgram_viz::Histgram::new(&counts[..counts.len() - 10]);
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

// TODO: This should be much more efficient!
fn fill_sparse_region(ds: &mut DataSet, config: &UnitConfig) {
    let mut edge_count: HashMap<_, Vec<_>> = HashMap::new();
    for edge in ds.encoded_reads.iter().flat_map(|r| r.edges.iter()) {
        // let len = if edge.label.is_empty() {
        //     edge.offset
        // } else {
        //     edge.label().len() as i64
        // };
        let key = (edge.from.min(edge.to), edge.from.max(edge.to));
        edge_count.entry(key).or_default().push(edge);
        // let entry = edge_count.entry(edge).or_default();
        // entry.0 += 1;
        // entry.1 += len;
    }
    // use std::collections::HashSet;
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
    let sparse_thr = config.chunk_len as i64;
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
        if sparse_thr < mean_len && count_thr < edges.len() {
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
    // let mut sparse_edge: HashSet<_> = edge_count
    //     .into_iter()
    //     .filter_map(|(edge, (count, len))| (len / count as i64 > sparse_thr).then(|| edge))
    //     .collect();
    // let stride = config.chunk_len + config.skip_len;
    // let picked_units: Vec<_> = {
    //     let mut new_units = vec![];
    //     for read in ds.encoded_reads.iter() {
    //         for edge in read.edges.iter() {
    //             let seq = edge.label();
    //             let edge = (edge.from.min(edge.to), edge.from.max(edge.to));
    //             if sparse_edge.contains(&edge) {
    //                 // Pick new units.
    //                 new_units.extend(
    //                     (0..)
    //                         .map(|i| (stride * i, stride * i + config.chunk_len))
    //                         .take_while(|&(_, y)| y + config.skip_len < seq.len())
    //                         .map(|(s, t)| &seq[s..t])
    //                         .filter(|u| !is_repetitive(u, config)),
    //                 );
    //                 sparse_edge.remove(&edge);
    //             }
    //         }
    //         if sparse_edge.is_empty() {
    //             break;
    //         }
    //     }
    //     new_units
    // };
    debug!("FillSparse\t{}", picked_units.len());
    let last_unit = ds.selected_chunks.iter().map(|x| x.id).max().unwrap();
    ds.selected_chunks
        .extend(picked_units.iter().enumerate().map(|(i, seq)| {
            let id = i as u64 + last_unit + 1;
            let seq = String::from_utf8_lossy(seq).to_string();
            Unit::new(id, seq, config.min_cluster)
        }));
}

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
    for ((unit, direction), labels) in tail_counts
        .into_iter()
        .filter(|(_, labels)| labels.len() > count_thr / 2 && !labels.is_empty())
    {
        let seq = labels.iter().max_by_key(|xs| xs.len()).unwrap();
        trace!("FillSparse\t{}\t{}\t{}", unit, direction, labels.len());
        let new_units = (0..)
            .map(|i| (config.chunk_len * i, config.chunk_len * (i + 1)))
            .take_while(|&(_, y)| y < seq.len())
            .map(|(s, t)| &seq[s..t])
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
    let mut count: HashMap<_, i32> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *count.entry(node.unit).or_default() += 1;
        }
    }
    let median = {
        let mut count: Vec<_> = count.values().copied().collect();
        let median = count.len() / 2;
        *count.select_nth_unstable(median).1 as i32
    };
    let (lower, upper) = match ds.read_type {
        ReadType::CCS => ((median / 10).max(2), 10 * median),
        ReadType::CLR => ((median / 10).max(4), 5 * median),
        ReadType::ONT => ((median / 10).max(4), 5 * median),
        ReadType::None => (median / 10, 200),
    };
    ds.selected_chunks.retain(|unit| {
        let coverage = match count.get(&unit.id) {
            Some(res) => *res,
            None => return false,
        };
        !to_be_removed[unit.id as usize] && lower < coverage && coverage < upper
    });
    let mut idx = 0;
    ds.selected_chunks.iter_mut().for_each(|unit| {
        unit.id = idx;
        idx += 1;
    });
    ds.encoded_reads.clear();
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

// #[allow(dead_code)]
// fn filter_unit(mut ds: DataSet, config: &UnitConfig) -> DataSet {
//     let mut counts: HashMap<_, usize> = HashMap::new();
//     let (lower, upper) = (config.lower_count, config.upper_count);
//     for read in ds.encoded_reads.iter() {
//         for node in read.nodes.iter() {
//             *counts.entry(node.unit).or_default() += 1;
//         }
//     }
//     let original_len = ds.selected_chunks.len();
//     let filtered_unit: HashMap<u64, u64> = counts
//         .iter()
//         .filter_map(|(&key, &val)| match val {
//             x if lower < x && x < upper => Some(key),
//             _ => None,
//         })
//         .enumerate()
//         .map(|(idx, key)| (key, idx as u64))
//         .collect();
//     debug!("# of units:{}=>{}", original_len, filtered_unit.len());
//     ds.selected_chunks
//         .retain(|chunk| filtered_unit.contains_key(&chunk.id));
//     // Never panic.
//     ds.selected_chunks
//         .iter_mut()
//         .for_each(|chunk| chunk.id = *filtered_unit.get(&chunk.id).unwrap());
//     let prev = ds.encoded_reads.len();
//     let node_prev = ds
//         .encoded_reads
//         .iter()
//         .map(|r| r.nodes.len())
//         .sum::<usize>();
//     ds.encoded_reads.iter_mut().for_each(|read| {
//         read.nodes
//             .retain(|node| filtered_unit.contains_key(&node.unit));
//         read.nodes
//             .iter_mut()
//             .for_each(|n| n.unit = *filtered_unit.get(&n.unit).unwrap());
//     });
//     ds.encoded_reads.retain(|read| !read.nodes.is_empty());
//     let now = ds.encoded_reads.len();
//     let node_now = ds
//         .encoded_reads
//         .iter()
//         .map(|r| r.nodes.len())
//         .sum::<usize>();
//     debug!(
//         "Encoded Reads{}->{}(Raw reads{})",
//         prev,
//         now,
//         ds.raw_reads.len()
//     );
//     debug!("Number of Unit {}->{}", node_prev, node_now);
//     ds
// }

// #[allow(dead_code)]
// fn to_u64(kmer: &[u8]) -> u64 {
//     let mut res = 0;
//     for x in kmer {
//         res = res << 2;
//         res += match x {
//             b'A' | b'a' => 0,
//             b'C' | b'c' => 1,
//             b'G' | b'g' => 2,
//             b'T' | b't' => 3,
//             _ => 0,
//         };
//     }
//     res
// }

// fn to_index(kmer: &[u8]) -> usize {
//     kmer.iter().fold(0, |x, y| match y {
//         b'A' => (x << 2),
//         b'C' => (x << 2) + 1,
//         b'G' => (x << 2) + 2,
//         b'T' => (x << 2) + 3,
//         _ => panic!(),
//     })
// }

// fn remove_overlapping_units_mm2(
//     mut units: Vec<Vec<u8>>,
//     config: &UnitConfig,
//     followup: bool,
// ) -> Option<Vec<Vec<u8>>> {
//     use rand::{thread_rng, Rng};
//     let mut rng = thread_rng();
//     let id: u64 = rng.gen::<u64>() % 100_000_000;
//     let mut c_dir = std::env::current_dir().ok()?;
//     c_dir.push(format!("{}", id));
//     debug!("Creating {:?}.", c_dir);
//     std::fs::create_dir(&c_dir).ok()?;
//     let unit = {
//         use bio_utils::fasta::Writer;
//         let mut path = c_dir.clone();
//         path.push("units.fa");
//         let mut wtr = std::fs::File::create(&path).map(Writer::new).ok()?;
//         for (i, unit) in units.iter().enumerate() {
//             let id = format!("{}", i);
//             let record = bio_utils::fasta::Record::with_data(&id, &None, unit);
//             wtr.write_record(&record).ok()?;
//         }
//         path.into_os_string().into_string().ok()?
//     };
//     let preset = if followup { "asm10" } else { "map-pb" };
//     let mm2 = crate::minimap2::minimap2(&unit, &unit, config.threads, preset, true, false);
//     let mm2: Vec<_> = String::from_utf8(mm2)
//         .unwrap()
//         .lines()
//         .filter_map(bio_utils::paf::PAF::new)
//         .collect();
//     debug!("{} Alignments.", mm2.len());
//     debug!("Removing {:?}", c_dir);
//     std::fs::remove_dir_all(c_dir).ok()?;
//     let mut edges: Vec<_> = (0..units.len()).map(|i| vec![false; i]).collect();
//     // TODO: Make as parameters.
//     for aln in mm2
//         .iter()
//         .filter(|a| a.matchnum > 500 && a.qname != a.tname)
//     {
//         let qname: usize = aln.qname.parse().ok()?;
//         let rname: usize = aln.tname.parse().ok()?;
//         let (i, j) = if qname < rname {
//             (rname, qname)
//         } else {
//             (qname, rname)
//         };
//         if i == j {
//             debug!("{:?}", aln);
//         }
//         edges[i][j] = true;
//     }
//     let num_edges = edges
//         .iter()
//         .map(|e| e.iter().filter(|&&b| b).count())
//         .sum::<usize>();
//     debug!("Graph constructed. {} edges.", num_edges);
//     let to_be_removed = approx_vertex_cover(edges, units.len());
//     let mut idx = 0;
//     units.retain(|_| {
//         idx += 1;
//         !to_be_removed[idx - 1]
//     });
//     debug!("Resulting {} units.", units.len());
//     Some(units)
// }

// #[allow(dead_code)]
// fn remove_overlapping_units(
//     mut units: Vec<Vec<u8>>,
//     config: &UnitConfig,
//     followup: bool,
// ) -> Vec<Vec<u8>> {
//     debug!("Collected {} units", units.len());
//     let k = config.k;
//     let kmer_vec: Vec<(Vec<_>, Vec<_>)> = units
//         .par_iter()
//         .map(|unit| {
//             let mut forward: Vec<_> = unit.windows(k).map(to_u64).collect();
//             forward.sort();
//             forward.dedup();
//             let mut reverse: Vec<_> = bio_utils::revcmp(unit).windows(k).map(to_u64).collect();
//             reverse.sort();
//             reverse.dedup();
//             (forward, reverse)
//         })
//         .collect();
//     let edges: Vec<Vec<bool>> = (0..units.len())
//         .into_par_iter()
//         .map(|i| {
//             let kmer = &kmer_vec[i].0;
//             (0..i)
//                 .into_par_iter()
//                 .map(|j| {
//                     let &(ref f, ref r) = &kmer_vec[j];
//                     if !followup {
//                         compute_jaccard(kmer, f) > config.jaccard_thr
//                             || compute_jaccard(kmer, r) > config.jaccard_thr
//                     } else {
//                         (compute_jaccard(kmer, f) > config.jaccard_thr
//                             && alignment(&units[i], &units[j]) > config.alignment_thr)
//                             || (compute_jaccard(kmer, r) > config.jaccard_thr
//                                 && alignment(&units[i], &bio_utils::revcmp(&units[j]))
//                                     > config.alignment_thr)
//                     }
//                 })
//                 .collect()
//         })
//         .collect();
//     let num_edges = edges
//         .iter()
//         .map(|e| e.iter().filter(|&&b| b).count())
//         .sum::<usize>();
//     debug!("Graph constructed. {} edges.", num_edges);
//     let to_be_removed = approx_vertex_cover(edges, units.len());
//     let mut idx = 0;
//     units.retain(|_| {
//         idx += 1;
//         !to_be_removed[idx - 1]
//     });
//     debug!("Resulting {} units.", units.len());
//     units
// }

// #[allow(dead_code)]
// fn compute_jaccard(kmer1: &[u64], kmer2: &[u64]) -> f64 {
//     let mut union = 0;
//     let mut intersection = 0;
//     let (mut p1, mut p2) = (0, 0);
//     while p1 < kmer1.len() && p2 < kmer2.len() {
//         union += 1;
//         use std::cmp::Ordering::*;
//         match kmer1[p1].cmp(&kmer2[p2]) {
//             Equal => {
//                 intersection += 1;
//                 p1 += 1;
//                 p2 += 1;
//             }
//             Less => {
//                 p1 += 1;
//             }
//             Greater => {
//                 p2 += 1;
//             }
//         }
//     }
//     union += (kmer1.len() - p1 - 1) + (kmer2.len() - p2 - 1);
//     intersection as f64 / union as f64
// }

// fn compute_jaccard(kmer_vec1: &[bool], kmer_vec2: &[bool], thr: f64) -> bool {
//     let (mut union, mut intersection) = (0, 0);
//     for (x, y) in kmer_vec1.iter().zip(kmer_vec2.iter()) {
//         union += (x | y) as u32;
//         intersection += (x & y) as u32;
//     }
//     intersection as f64 / union as f64 > thr
// }

// Overlapping alignment.
// #[allow(dead_code)]
// fn alignment(seq1: &[u8], seq2: &[u8]) -> f64 {
//     let mut dp = vec![vec![0; seq2.len() + 1]; seq1.len() + 1];
//     for i in 1..seq1.len() + 1 {
//         for j in 1..seq2.len() + 1 {
//             let is_match = seq1[i - 1].to_ascii_uppercase() == seq2[j - 1].to_ascii_uppercase();
//             let m = if is_match { 1 } else { -1 };
//             dp[i][j] = (dp[i - 1][j] - 1)
//                 .max(dp[i][j - 1] - 1)
//                 .max(dp[i - 1][j - 1] + m);
//         }
//     }
//     let column_max = dp.iter().map(|x| x[seq2.len()]).max().unwrap();
//     let row_max = *dp[seq1.len()].iter().max().unwrap();
//     column_max.max(row_max) as f64 / ((seq1.len() * seq2.len()) as f64).sqrt()
// }
