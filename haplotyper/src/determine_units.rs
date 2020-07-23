use rayon::prelude::*;
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitConfig {
    pub chunk_len: usize,
    pub skip_len: usize,
    pub margin: usize,
    pub k: usize,
    pub jaccard_thr: f64,
    pub alignment_thr: f64,
    pub hits_thr: f64,
    pub hits_range: usize,
}

const DEFAULT_JCD: f64 = 0.7;
const DEFAULT_ALN: f64 = 0.95;
const DEFAULT_RNG: usize = 200;
const DEFAULT_HITS: f64 = 0.8;
impl UnitConfig {
    pub fn new_ccs(chunk_len: usize, skip_len: usize, margin: usize) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 7,
            jaccard_thr: DEFAULT_JCD,
            alignment_thr: DEFAULT_ALN,
            hits_range: DEFAULT_RNG,
            hits_thr: DEFAULT_HITS,
        }
    }
    pub fn new_clr(chunk_len: usize, skip_len: usize, margin: usize) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 4,
            jaccard_thr: 0.5,
            alignment_thr: 0.8,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.6,
        }
    }
    pub fn new_ont(chunk_len: usize, skip_len: usize, margin: usize) -> Self {
        Self {
            chunk_len,
            skip_len,
            margin,
            k: 4,
            jaccard_thr: 0.5,
            alignment_thr: 0.8,
            hits_range: DEFAULT_RNG,
            hits_thr: 0.6,
        }
    }
}

pub trait DetermineUnit {
    // fn select_chunks_num(self, config: &UnitConfig) -> Self;
    fn select_chunks_freq(self, config: &UnitConfig) -> Self;
}

use definitions::*;
impl DetermineUnit for definitions::DataSet {
    // fn select_chunks_num(mut self, config: &UnitConfig) -> Self {
    //     rayon::ThreadPoolBuilder::new()
    //         .num_threads(config.threads)
    //         .build_global()
    //         .unwrap();
    //     let mut reads: Vec<&RawRead> = self.raw_reads.iter().collect();
    //     reads.sort_by_key(|r| r.seq().len());
    //     reads.reverse();
    //     debug!("Reads are sorted. {} reads.", reads.len());
    //     let selected_chunks: Vec<Vec<u8>> = reads
    //         .iter()
    //         .flat_map(|r| split_into(r, config))
    //         .take(config.chunk_num)
    //         .map(|e| e.to_vec())
    //         .collect();
    //     debug!(
    //         "Units collected. {}/{} units.",
    //         selected_chunks.len(),
    //         config.chunk_num
    //     );
    //     debug!("Removing overlapping units");
    //     let selected_chunks = remove_overlapped_units(selected_chunks, config.k);
    //     debug!("Resulting {} units.", selected_chunks.len());
    //     self.selected_chunks = selected_chunks
    //         .into_iter()
    //         .enumerate()
    //         .map(|(idx, seq)| {
    //             let id = idx as u64;
    //             let seq = String::from_utf8(seq).unwrap();
    //             Unit { id, seq }
    //         })
    //         .collect();
    //     self
    // }
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
            let is_unique =
                units
                    .par_iter()
                    .zip(forward_kmer_vec.par_iter())
                    .all(|(unit2, finger)| {
                        let forward_jaccard = compute_jaccard(finger, &forward_finger);
                        let reverse_jaccard = compute_jaccard(finger, &reverse_finger);
                        if forward_jaccard > config.jaccard_thr {
                            alignment(&unit, unit2) < config.alignment_thr
                        } else if reverse_jaccard > config.jaccard_thr {
                            let unit2 = bio_utils::revcmp(unit2);
                            alignment(&unit, &unit2) < config.alignment_thr
                        } else {
                            true
                        }
                    });
            hits.push_back(is_unique);
            if hits.len() > config.hits_range {
                hits.pop_front();
            }
            if is_unique {
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

// fn remove_overlapped_units(mut units: Vec<Vec<u8>>, k: usize) -> Vec<Vec<u8>> {
//     let forward_kmer_vec: Vec<_> = units
//         .iter()
//         .map(|unit| {
//             let mut fingerprint = vec![false; 4usize.pow(k as u32)];
//             for kmer in unit.windows(k) {
//                 fingerprint[to_index(kmer)] = true;
//             }
//             fingerprint
//         })
//         .collect();
//     let reverse_kmer_vec: Vec<_> = units
//         .iter()
//         .map(|unit| {
//             let mut fingerprint = vec![false; 4usize.pow(k as u32)];
//             let unit = bio_utils::revcmp(unit);
//             for kmer in unit.windows(k) {
//                 fingerprint[to_index(kmer)] = true;
//             }
//             fingerprint
//         })
//         .collect();
//     let result: Vec<_> = units
//         .par_iter()
//         .enumerate()
//         .flat_map(|(idx, unit1)| {
//             units
//                 .iter()
//                 .enumerate()
//                 .skip(idx + 1)
//                 .filter_map(|(jdx, unit2)| {
//                     let forward_jaccard =
//                         compute_jaccard(&forward_kmer_vec[idx], &forward_kmer_vec[jdx]);
//                     let reverse_jaccard =
//                         compute_jaccard(&forward_kmer_vec[idx], &reverse_kmer_vec[jdx]);
//                     if forward_jaccard > 0.7 {
//                         let aln = alignment(unit1, unit2);
//                         Some((idx, jdx, aln))
//                     } else if reverse_jaccard > 0.7 {
//                         let unit2 = bio_utils::revcmp(unit2);
//                         let aln = alignment(unit1, &unit2);
//                         Some((idx, jdx, aln))
//                     } else {
//                         None
//                     }
//                 })
//                 .collect::<Vec<_>>()
//         })
//         .collect();

//     let retain_nodes = {
//         let mut temp = vec![true; units.len()];
//         for (idx, jdx, aln) in result.into_iter() {
//             debug!("ALN\t{}\t{}\t{}", idx, jdx, aln);
//             if aln > 0.95 {
//                 temp[idx.max(jdx)] = false;
//             }
//         }
//         temp
//     };
//     let mut idx = 0;
//     units.retain(|_| {
//         idx += 1;
//         retain_nodes[idx - 1]
//     });
//     units
// }

fn compute_jaccard(kmer_vec1: &[bool], kmer_vec2: &[bool]) -> f64 {
    let (union, intersection) =
        kmer_vec1
            .iter()
            .zip(kmer_vec2.iter())
            .fold((0, 0), |(union, inter), (&x, &y)| match (x, y) {
                (true, true) => (union + 1, inter + 1),
                (false, true) | (true, false) => (union + 1, inter),
                (false, false) => (union, inter),
            });
    intersection as f64 / union as f64
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
