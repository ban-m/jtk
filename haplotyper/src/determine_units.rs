use super::polish_units::PolishUnit;
use super::polish_units::PolishUnitConfig;
// use super::Encode;
use definitions::*;
// use rayon::prelude::*;
use std::collections::HashMap;
#[derive(Debug, Clone)]
pub struct UnitConfig {
    pub chunk_len: usize,
    pub skip_len: usize,
    pub unit_num: usize,
    pub margin: usize,
    pub k: usize,
    pub jaccard_thr: f64,
    pub alignment_thr: f64,
    pub threads: usize,
    pub min_cluster: usize,
    pub exclude_repeats: f64,
    pub upper_count: usize,
    pub lower_count: usize,
}

// const DEFAULT_RNG: usize = 200;
impl UnitConfig {
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
            k: 9,
            unit_num,
            jaccard_thr: 0.2,
            alignment_thr: 0.3,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
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
            k: 6,
            unit_num,
            jaccard_thr: 0.15,
            alignment_thr: 0.4,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
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
            k: 6,
            unit_num,
            jaccard_thr: 0.15,
            alignment_thr: 0.4,
            threads,
            min_cluster: 2,
            exclude_repeats,
            upper_count,
            lower_count,
        }
    }
}

pub trait DetermineUnit {
    fn select_chunks(self, config: &UnitConfig) -> Self;
}

impl DetermineUnit for definitions::DataSet {
    fn select_chunks(mut self, config: &UnitConfig) -> Self {
        let mut reads: Vec<&RawRead> = self.raw_reads.iter().collect();
        reads.sort_by_key(|r| r.seq().len());
        reads.reverse();
        debug!("Select Unit: Configuration:{:?}", config);
        let selected_chunks = if self.read_type == ReadType::CCS {
            reads
                .iter()
                .flat_map(|r| split_into(r, config))
                .filter(|u| !is_repetitive(u, config))
                .take(config.unit_num)
                .enumerate()
                .map(|(idx, seq)| Unit {
                    id: idx as u64,
                    seq: String::from_utf8_lossy(seq).to_string(),
                    cluster_num: config.min_cluster,
                })
                .collect()
        } else {
            let clr_config = {
                let mut temp = config.clone();
                temp.chunk_len = 12 * temp.chunk_len / 10;
                debug!("Calib length {}->{}", config.chunk_len, temp.chunk_len);
                temp
            };
            self.selected_chunks = reads
                .iter()
                .flat_map(|r| split_into(r, &clr_config))
                .filter(|u| !is_repetitive(u, config))
                .take(clr_config.unit_num)
                .enumerate()
                .map(|(idx, seq)| Unit {
                    id: idx as u64,
                    seq: String::from_utf8_lossy(seq).to_string(),
                    cluster_num: config.min_cluster,
                })
                .collect();
            let seed = 20;
            debug!("UNITNUM\t{}\tPICKED", self.selected_chunks.len());
            self = crate::encode::encode_by_mm2(self, config.threads).unwrap();
            let polish_config = PolishUnitConfig::new(ReadType::CLR, 3, 10, seed);
            self = self.polish_unit(&polish_config);
            self = crate::encode::encode_by_mm2(self, config.threads).unwrap();
            let polish_config = PolishUnitConfig::new(ReadType::CLR, 5, 10, seed);
            self = self.polish_unit(&polish_config);
            self = crate::encode::encode_by_mm2(self, config.threads).unwrap();
            let polish_config = PolishUnitConfig::new(ReadType::CLR, 10, 10, seed);
            self = self.polish_unit(&polish_config);
            debug!("UNITNUM\t{}\tPOLISHED", self.selected_chunks.len());
            self.selected_chunks
                .into_iter()
                .filter(|unit| unit.seq.len() > config.chunk_len)
                .map(|mut unit| {
                    while unit.seq.len() > config.chunk_len {
                        unit.seq.pop();
                    }
                    unit
                })
                .enumerate()
                .map(|(idx, unit)| Unit {
                    id: idx as u64,
                    ..unit
                })
                .collect()
        };
        self.selected_chunks = selected_chunks;
        debug!("UNITNUM\t{}\tRAWUNIT", self.selected_chunks.len());
        self = crate::encode::encode_by_mm2(self, config.threads).unwrap();
        self = filter_unit_by_ovlp(self, config);
        debug!("UNITNUM\t{}\tFILTERED", self.selected_chunks.len());
        self = crate::encode::encode_by_mm2(self, config.threads).unwrap();
        let polish_config = PolishUnitConfig::new(ReadType::CLR, 10, 10, 24);
        self = self.polish_unit(&polish_config);
        debug!("UNITNUM\t{}\tPOLISHED", self.selected_chunks.len());
        let mut idx = 0;
        self.selected_chunks.iter_mut().for_each(|mut unit| {
            while unit.seq.len() > config.chunk_len {
                unit.seq.pop();
            }
            unit.id = idx;
            idx += 1;
        });
        self
    }
}
// fn raw_read_to_unit(ds: &DataSet, p: usize) -> std::io::Result<Vec<bio_utils::paf::PAF>> {
//     use rand::{thread_rng, Rng};
//     let mut rng = thread_rng();
//     let id: u64 = rng.gen::<u64>() % 100_000_000;
//     let mut c_dir = std::env::current_dir()?;
//     c_dir.push(format!("{}", id));
//     debug!("Creating {:?}.", c_dir);
//     std::fs::create_dir(&c_dir)?;
//     use std::io::{BufWriter, Write};
//     // Create reference and reads.
//     let (reference, reads) = {
//         let mut reference = c_dir.clone();
//         reference.push("units.fa");
//         let mut wtr = std::fs::File::create(&reference).map(BufWriter::new)?;
//         for unit in ds.selected_chunks.iter() {
//             writeln!(&mut wtr, ">{}\n{}", unit.id, &unit.seq)?;
//         }
//         let mut reads = c_dir.clone();
//         reads.push("reads.fa");
//         let mut wtr = std::fs::File::create(&reads).map(BufWriter::new)?;
//         for read in ds.raw_reads.iter() {
//             writeln!(&mut wtr, ">{}\n{}", read.name, &read.seq)?;
//         }
//         let reference = reference.into_os_string().into_string().unwrap();
//         let reads = reads.into_os_string().into_string().unwrap();
//         (reference, reads)
//     };
//     use crate::minimap2;
//     let threads = format!("{}", p);
//     let args = vec!["-H", "-k", "10", "-t", &threads, "-c"];
//     let mm2 = minimap2::minimap2_args(&reference, &reads, &args);
//     debug!("Removing {:?}", c_dir);
//     std::fs::remove_dir_all(c_dir)?;
//     let thr = 100;
//     let alignments: Vec<_> = String::from_utf8_lossy(&mm2)
//         .lines()
//         .filter_map(|l| bio_utils::paf::PAF::new(&l))
//         .filter(|a| a.tstart < thr && a.tlen - a.tend < thr)
//         .collect();
//     Ok(alignments)
// }

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

fn filter_unit_by_ovlp(mut ds: DataSet, config: &UnitConfig) -> DataSet {
    let unit_len = ds.selected_chunks.len();
    // Maybe it became infeasible due to O(N^2) allcation would be occured.
    debug!("UNITNUM\t{}", unit_len);
    let mut edges: Vec<_> = (0..unit_len).map(|i| vec![false; i]).collect();
    for read in ds.encoded_reads.iter() {
        for (i, node) in read.nodes.iter().enumerate() {
            for mode in read.nodes.iter().skip(i + 1) {
                let node_end = node.position_from_start + node.seq.as_bytes().len();
                let mode_start = mode.position_from_start;
                let ovlp_len = node_end.max(mode_start) - mode_start;
                if 500 < ovlp_len {
                    let (i, j) = (node.unit.max(mode.unit), node.unit.min(mode.unit));
                    if i != j {
                        edges[i as usize][j as usize] = true;
                    }
                }
            }
        }
    }
    let to_be_removed = approx_vertex_cover(edges, ds.selected_chunks.len());
    let survived_unit = to_be_removed.iter().filter(|&&b| !b).count();
    debug!("UNITNUM\t{}", survived_unit);
    let mut count: HashMap<_, _> = HashMap::new();
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *count.entry(node.unit).or_default() += 1;
        }
    }
    ds.selected_chunks.retain(|unit| {
        let coverage = match count.get(&unit.id) {
            Some(res) => *res,
            None => return false,
        };
        !to_be_removed[unit.id as usize]
            && config.lower_count < coverage
            && coverage < config.upper_count
    });
    debug!("UNITNUM\t{}", ds.selected_chunks.len());
    let mut idx = 0;
    ds.selected_chunks.iter_mut().for_each(|unit| {
        unit.id = idx;
        idx += 1;
    });
    ds.encoded_reads.clear();
    ds
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

#[allow(dead_code)]
fn filter_unit(mut ds: DataSet, config: &UnitConfig) -> DataSet {
    let mut counts: HashMap<_, usize> = HashMap::new();
    let (lower, upper) = (config.lower_count, config.upper_count);
    for read in ds.encoded_reads.iter() {
        for node in read.nodes.iter() {
            *counts.entry(node.unit).or_default() += 1;
        }
    }
    let original_len = ds.selected_chunks.len();
    let filtered_unit: HashMap<u64, u64> = counts
        .iter()
        .filter_map(|(&key, &val)| match val {
            x if lower < x && x < upper => Some(key),
            _ => None,
        })
        .enumerate()
        .map(|(idx, key)| (key, idx as u64))
        .collect();
    debug!("# of units:{}=>{}", original_len, filtered_unit.len());
    ds.selected_chunks
        .retain(|chunk| filtered_unit.contains_key(&chunk.id));
    // Never panic.
    ds.selected_chunks
        .iter_mut()
        .for_each(|chunk| chunk.id = *filtered_unit.get(&chunk.id).unwrap());
    let prev = ds.encoded_reads.len();
    let node_prev = ds
        .encoded_reads
        .iter()
        .map(|r| r.nodes.len())
        .sum::<usize>();
    ds.encoded_reads.iter_mut().for_each(|read| {
        read.nodes
            .retain(|node| filtered_unit.contains_key(&node.unit));
        read.nodes
            .iter_mut()
            .for_each(|n| n.unit = *filtered_unit.get(&n.unit).unwrap());
    });
    ds.encoded_reads.retain(|read| !read.nodes.is_empty());
    let now = ds.encoded_reads.len();
    let node_now = ds
        .encoded_reads
        .iter()
        .map(|r| r.nodes.len())
        .sum::<usize>();
    debug!(
        "Encoded Reads{}->{}(Raw reads{})",
        prev,
        now,
        ds.raw_reads.len()
    );
    debug!("Number of Unit {}->{}", node_prev, node_now);
    ds
}

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
