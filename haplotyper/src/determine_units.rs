use super::polish_units::PolishUnit;
use super::polish_units::PolishUnitConfig;
use super::Encode;
use super::ReadType;
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
            k: 7,
            unit_num,
            jaccard_thr: 0.2,
            alignment_thr: 0.4,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::CCS,
            threads,
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
            k: 7,
            unit_num,
            jaccard_thr: 0.2,
            alignment_thr: 0.35,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::CLR,
            threads,
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
            alignment_thr: 0.1,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.3,
            read_type: ReadType::ONT,
            threads,
        }
    }
}

pub trait DetermineUnit {
    fn select_chunks_freq(self, config: &UnitConfig) -> Self;
    fn select_chunks(self, config: &UnitConfig) -> Self;
}

use definitions::*;
impl DetermineUnit for definitions::DataSet {
    fn select_chunks(mut self, config: &UnitConfig) -> Self {
        let mut reads: Vec<&RawRead> = self.raw_reads.iter().collect();
        reads.sort_by_key(|r| r.seq().len());
        reads.reverse();
        let mut units: Vec<_> = if config.read_type == ReadType::CCS {
            reads
                .iter()
                .flat_map(|r| split_into(r, config))
                .take(config.unit_num)
                .map(|e| e.to_vec())
                .collect()
        } else {
            let clr_config = {
                let mut temp = config.clone();
                temp.chunk_len = temp.chunk_len * 110 / 100;
                debug!("Calib length {}->{}", config.chunk_len, temp.chunk_len);
                temp
            };
            let polish_config = PolishUnitConfig::new("CLR", 7, 10);
            let units: Vec<_> = reads
                .iter()
                .flat_map(|r| split_into(r, &clr_config))
                .take(clr_config.unit_num)
                .map(|e| e.to_vec())
                .collect();
            self.selected_chunks = units
                .iter()
                .enumerate()
                .map(|(idx, seq)| {
                    let id = idx as u64;
                    let seq = String::from_utf8_lossy(seq).to_string();
                    Unit { id, seq }
                })
                .collect();
            self = self.encode(config.threads);
            self = self.polish_unit(&polish_config);
            self = self.encode(config.threads);
            self = self.polish_unit(&polish_config);
            debug!("Polished.");
            self.selected_chunks
                .iter()
                .filter(|e| e.seq.len() > config.chunk_len)
                .map(|e| e.seq.as_bytes()[..config.chunk_len].to_vec())
                .collect()
        };
        debug!("Collected {} units", units.len());
        let k = config.k;
        let kmer_vec: Vec<(Vec<_>, Vec<_>)> = units
            .par_iter()
            .map(|unit| {
                let mut forward_finger = vec![false; 4usize.pow(k as u32)];
                for kmer in unit.windows(k) {
                    forward_finger[to_index(kmer)] = true;
                }
                let mut reverse_finger = vec![false; 4usize.pow(k as u32)];
                let unit = bio_utils::revcmp(unit);
                for kmer in unit.windows(k) {
                    reverse_finger[to_index(kmer)] = true;
                }
                (forward_finger, reverse_finger)
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
                        if compute_jaccard(kmer, f, config.jaccard_thr) {
                            let aln = alignment(&units[i], &units[j]);
                            aln > config.alignment_thr
                        } else if compute_jaccard(kmer, r, config.jaccard_thr) {
                            let rev = bio_utils::revcmp(&units[j]);
                            let aln = alignment(&units[i], &rev);
                            aln > config.alignment_thr
                        } else {
                            false
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
        self.selected_chunks = units
            .into_iter()
            .enumerate()
            .map(|(idx, seq)| {
                let id = idx as u64;
                let seq = String::from_utf8(seq).unwrap();
                Unit { id, seq }
            })
            .collect();
        self
    }
    fn select_chunks_freq(mut self, config: &UnitConfig) -> Self {
        let mut reads: Vec<&RawRead> = self.raw_reads.iter().collect();
        reads.sort_by_key(|r| r.seq().len());
        reads.reverse();
        let mut units: Vec<Vec<u8>> = vec![];
        let mut forward_kmer_vec: Vec<Vec<bool>> = vec![];
        let k = config.k;
        let mut hits = std::collections::VecDeque::new();
        for unit in reads.iter().flat_map(|r| split_into(r, config)) {
            let mut forward_finger = vec![false; 4usize.pow(k as u32)];
            for kmer in unit.windows(k) {
                forward_finger[to_index(kmer)] = true;
            }
            let mut reverse_finger = vec![false; 4usize.pow(k as u32)];
            let unit = bio_utils::revcmp(unit);
            for kmer in unit.windows(k) {
                reverse_finger[to_index(kmer)] = true;
            }
            let is_duplicate =
                units
                    .par_iter()
                    .zip(forward_kmer_vec.par_iter())
                    .any(|(_unit2, finger)| {
                        let thr = config.jaccard_thr;
                        compute_jaccard(finger, &forward_finger, thr)
                            || compute_jaccard(finger, &reverse_finger, thr)
                    });
            hits.push_back(!is_duplicate);
            if hits.len() > config.hits_range {
                hits.pop_front();
            }
            if !is_duplicate {
                units.push(unit.to_vec());
                forward_kmer_vec.push(forward_finger);
            }
            let hitting_rate = hits.iter().filter(|&&b| b).count() as f64 / hits.len() as f64;
            if hitting_rate < config.hits_thr {
                break;
            }
        }
        debug!("Resulting {} units.", units.len());
        self.selected_chunks = units
            .into_iter()
            .enumerate()
            .map(|(idx, seq)| {
                let id = idx as u64;
                let seq = String::from_utf8(seq).unwrap();
                Unit { id, seq }
            })
            .collect();
        self
    }
}

fn compute_jaccard(kmer_vec1: &[bool], kmer_vec2: &[bool], thr: f64) -> bool {
    let (mut union, mut intersection) = (0, 0);
    for (x, y) in kmer_vec1.iter().zip(kmer_vec2.iter()) {
        union += (x | y) as u32;
        intersection += (x & y) as u32;
    }
    intersection as f64 / union as f64 > thr
}

fn to_index(kmer: &[u8]) -> usize {
    kmer.iter().fold(0, |x, y| match y {
        b'A' => (x << 2),
        b'C' => (x << 2) + 1,
        b'G' => (x << 2) + 2,
        b'T' => (x << 2) + 3,
        _ => panic!(),
    })
}

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
