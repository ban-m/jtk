use super::polish_units::PolishUnit;
use super::polish_units::PolishUnitConfig;
use super::Encode;
use definitions::*;
use rayon::prelude::*;
#[derive(Debug, Clone)]
pub struct UnitConfig {
    pub chunk_len: usize,
    pub skip_len: usize,
    pub unit_num: usize,
    pub margin: usize,
    pub k: usize,
    pub jaccard_thr: f64,
    pub alignment_thr: f64,
    pub hits_thr: f64,
    pub hits_range: usize,
    pub read_type: ReadType,
    pub threads: usize,
    pub min_cluster: usize,
}
const DEFAULT_RNG: usize = 200;
impl UnitConfig {
    pub fn new_ccs(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 9,
            unit_num,
            jaccard_thr: 0.2,
            alignment_thr: 0.3,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::CCS,
            threads,
            min_cluster: 2,
        }
    }
    pub fn new_clr(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 9,
            unit_num,
            jaccard_thr: 0.1,
            alignment_thr: 0.3,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::CLR,
            threads,
            min_cluster: 2,
        }
    }
    pub fn new_ont(
        chunk_len: usize,
        unit_num: usize,
        skip_len: usize,
        margin: usize,
        threads: usize,
    ) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 6,
            unit_num,
            jaccard_thr: 0.15,
            alignment_thr: 0.4,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::ONT,
            threads,
            min_cluster: 2,
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
        let units: Vec<_> = if config.read_type == ReadType::CCS {
            let units: Vec<_> = reads
                .iter()
                .flat_map(|r| split_into(r, config))
                .take(config.unit_num)
                .map(|e| e.to_vec())
                .collect();
            let units = {
                let mut config = config.clone();
                config.jaccard_thr /= 4.;
                //remove_overlapping_units(units, &config, false)
                remove_overlapping_units_mm2(units, &config, false).unwrap()
            };
            //remove_overlapping_units(units, config, true)
            remove_overlapping_units_mm2(units, config, true).unwrap()
        } else {
            let clr_config = {
                let mut temp = config.clone();
                temp.chunk_len = temp.chunk_len * 110 / 100;
                debug!("Calib length {}->{}", config.chunk_len, temp.chunk_len);
                temp
            };
            let units: Vec<_> = reads
                .iter()
                .flat_map(|r| split_into(r, &clr_config))
                .take(clr_config.unit_num)
                .map(|e| e.to_vec())
                .collect();
            debug!("Take {} units.", units.len());
            // let units = {
            //     let mut config = config.clone();
            //     config.jaccard_thr /= 4.;
            //     //remove_overlapping_units(units, &config, false)
            //     remove_overlapping_units_mm2(units, &config, false).unwrap()
            // };
            self.selected_chunks = units
                .iter()
                .enumerate()
                .map(|(idx, seq)| Unit {
                    id: idx as u64,
                    seq: String::from_utf8_lossy(seq).to_string(),
                    cluster_num: config.min_cluster,
                })
                .collect();
            // Minimap2 polishing.
            self = self.encode(config.threads, false);
            let polish_config = PolishUnitConfig::new(ReadType::CLR, 10, 10);
            self = self.polish_unit(&polish_config);
            debug!("Polished.");
            let mut units = vec![];
            while let Some(unit) = self.selected_chunks.pop() {
                if unit.seq.len() > config.chunk_len {
                    let Unit { seq, .. } = unit;
                    units.push(seq.into_bytes());
                }
            }
            //remove_overlapping_units(units, &config, true)
            remove_overlapping_units_mm2(units, &config, true).unwrap()
        };
        // TODO: Make as parameters.
        self.selected_chunks = units
            .iter()
            .enumerate()
            .map(|(idx, seq)| Unit {
                id: idx as u64,
                seq: String::from_utf8_lossy(seq).to_string(),
                cluster_num: config.min_cluster,
            })
            .collect();
        self = self.encode(config.threads, false);
        use std::collections::HashMap;
        let mut count: HashMap<_, u32> = HashMap::new();
        for r in self.encoded_reads.iter() {
            for n in r.nodes.iter() {
                *count.entry(n.unit).or_default() += 1;
            }
        }
        // Cut upper/lower 1 percent units.
        let mut covs: Vec<_> = count.values().copied().collect();
        covs.sort();
        let (lower, upper) = (covs[covs.len() / 100], covs[covs.len() * 99 / 100]);
        let lower = lower.max(5);
        debug!("Taking{}-{}", lower, upper);
        self.selected_chunks = self
            .selected_chunks
            .into_iter()
            .filter(|unit| match count.get(&unit.id) {
                Some(&count) => lower < count && count < upper,
                None => false,
            })
            .enumerate()
            .map(|(idx, unit)| Unit {
                id: idx as u64,
                ..unit
            })
            .collect();
        debug!("Final:{} units", self.selected_chunks.len());
        self
    }
}

#[allow(dead_code)]
fn to_u64(kmer: &[u8]) -> u64 {
    let mut res = 0;
    for x in kmer {
        res = res << 2;
        res += match x {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0,
        };
    }
    res
}

// fn to_index(kmer: &[u8]) -> usize {
//     kmer.iter().fold(0, |x, y| match y {
//         b'A' => (x << 2),
//         b'C' => (x << 2) + 1,
//         b'G' => (x << 2) + 2,
//         b'T' => (x << 2) + 3,
//         _ => panic!(),
//     })
// }

fn remove_overlapping_units_mm2(
    mut units: Vec<Vec<u8>>,
    config: &UnitConfig,
    followup: bool,
) -> Option<Vec<Vec<u8>>> {
    use rand::{thread_rng, Rng};
    let mut rng = thread_rng();
    let id: u64 = rng.gen::<u64>() % 100_000_000;
    let mut c_dir = std::env::current_dir().ok()?;
    c_dir.push(format!("{}", id));
    debug!("Creating {:?}.", c_dir);
    std::fs::create_dir(&c_dir).ok()?;
    let unit = {
        use bio_utils::fasta::Writer;
        let mut path = c_dir.clone();
        path.push("units.fa");
        let mut wtr = std::fs::File::create(&path).map(Writer::new).ok()?;
        for (i, unit) in units.iter().enumerate() {
            let id = format!("{}", i);
            let record = bio_utils::fasta::Record::with_data(&id, &None, unit);
            wtr.write_record(&record).ok()?;
        }
        path.into_os_string().into_string().ok()?
    };
    let preset = if config.read_type == ReadType::CCS || followup {
        "asm10"
    } else {
        "map-pb"
    };
    let mm2 = crate::minimap2::minimap2(&unit, &unit, config.threads, preset, true, false);
    let mm2: Vec<_> = String::from_utf8(mm2)
        .unwrap()
        .lines()
        .filter_map(bio_utils::paf::PAF::new)
        .collect();
    debug!("{} Alignments.", mm2.len());
    debug!("Removing {:?}", c_dir);
    std::fs::remove_dir_all(c_dir).ok()?;
    let mut edges: Vec<_> = (0..units.len()).map(|i| vec![false; i]).collect();
    // TODO: Make as parameters.
    for aln in mm2
        .iter()
        .filter(|a| a.matchnum > 1000 && a.qname != a.tname)
    {
        let qname: usize = aln.qname.parse().ok()?;
        let rname: usize = aln.tname.parse().ok()?;
        let thr = 100;
        let is_overlap = if aln.relstrand {
            (aln.qstart < thr && aln.tlen - aln.tend < thr)
                || (aln.qlen - aln.qend < thr && aln.tstart < thr)
        } else {
            (aln.qlen - aln.qend < thr && aln.tlen - aln.tend < thr)
                || (aln.qstart < thr && aln.tstart < thr)
        };
        if is_overlap {
            let (i, j) = if qname < rname {
                (rname, qname)
            } else {
                (qname, rname)
            };
            edges[i][j] = true;
        }
    }
    let num_edges = edges
        .iter()
        .map(|e| e.iter().filter(|&&b| b).count())
        .sum::<usize>();
    debug!("Graph constructed. {} edges.", num_edges);
    let to_be_removed = approx_vertex_cover(edges, units.len());
    let mut idx = 0;
    units.retain(|_| {
        idx += 1;
        !to_be_removed[idx - 1]
    });
    debug!("Resulting {} units.", units.len());
    Some(units)
}

#[allow(dead_code)]
fn remove_overlapping_units(
    mut units: Vec<Vec<u8>>,
    config: &UnitConfig,
    followup: bool,
) -> Vec<Vec<u8>> {
    debug!("Collected {} units", units.len());
    let k = config.k;
    let kmer_vec: Vec<(Vec<_>, Vec<_>)> = units
        .par_iter()
        .map(|unit| {
            let mut forward: Vec<_> = unit.windows(k).map(to_u64).collect();
            forward.sort();
            forward.dedup();
            let mut reverse: Vec<_> = bio_utils::revcmp(unit).windows(k).map(to_u64).collect();
            reverse.sort();
            reverse.dedup();
            (forward, reverse)
        })
        .collect();
    let edges: Vec<Vec<bool>> = (0..units.len())
        .into_par_iter()
        .map(|i| {
            let kmer = &kmer_vec[i].0;
            (0..i)
                .into_par_iter()
                .map(|j| {
                    let &(ref f, ref r) = &kmer_vec[j];
                    if !followup {
                        compute_jaccard(kmer, f) > config.jaccard_thr
                            || compute_jaccard(kmer, r) > config.jaccard_thr
                    } else {
                        (compute_jaccard(kmer, f) > config.jaccard_thr
                            && alignment(&units[i], &units[j]) > config.alignment_thr)
                            || (compute_jaccard(kmer, r) > config.jaccard_thr
                                && alignment(&units[i], &bio_utils::revcmp(&units[j]))
                                    > config.alignment_thr)
                    }
                })
                .collect()
        })
        .collect();
    let num_edges = edges
        .iter()
        .map(|e| e.iter().filter(|&&b| b).count())
        .sum::<usize>();
    debug!("Graph constructed. {} edges.", num_edges);
    let to_be_removed = approx_vertex_cover(edges, units.len());
    let mut idx = 0;
    units.retain(|_| {
        idx += 1;
        !to_be_removed[idx - 1]
    });
    debug!("Resulting {} units.", units.len());
    units
}

#[allow(dead_code)]
fn compute_jaccard(kmer1: &[u64], kmer2: &[u64]) -> f64 {
    let mut union = 0;
    let mut intersection = 0;
    let (mut p1, mut p2) = (0, 0);
    while p1 < kmer1.len() && p2 < kmer2.len() {
        union += 1;
        use std::cmp::Ordering::*;
        match kmer1[p1].cmp(&kmer2[p2]) {
            Equal => {
                intersection += 1;
                p1 += 1;
                p2 += 1;
            }
            Less => {
                p1 += 1;
            }
            Greater => {
                p2 += 1;
            }
        }
    }
    union += (kmer1.len() - p1 - 1) + (kmer2.len() - p2 - 1);
    intersection as f64 / union as f64
}

// fn compute_jaccard(kmer_vec1: &[bool], kmer_vec2: &[bool], thr: f64) -> bool {
//     let (mut union, mut intersection) = (0, 0);
//     for (x, y) in kmer_vec1.iter().zip(kmer_vec2.iter()) {
//         union += (x | y) as u32;
//         intersection += (x & y) as u32;
//     }
//     intersection as f64 / union as f64 > thr
// }

// Overlapping alignment.
#[allow(dead_code)]
fn alignment(seq1: &[u8], seq2: &[u8]) -> f64 {
    let mut dp = vec![vec![0; seq2.len() + 1]; seq1.len() + 1];
    for i in 1..seq1.len() + 1 {
        for j in 1..seq2.len() + 1 {
            let m = if seq1[i - 1] == seq2[j - 1] { 1 } else { -1 };
            dp[i][j] = (dp[i - 1][j] - 1)
                .max(dp[i][j - 1] - 1)
                .max(dp[i - 1][j - 1] + m);
        }
    }
    let column_max = dp.iter().map(|x| x[seq2.len()]).max().unwrap();
    let row_max = *dp[seq1.len()].iter().max().unwrap();
    column_max.max(row_max) as f64 / ((seq1.len() * seq2.len()) as f64).sqrt()
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
        degrees[argmax] = 0;
        edges[argmax].iter_mut().enumerate().for_each(|(i, b)| {
            degrees[i] -= *b as usize;
            *b = false;
        });
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
