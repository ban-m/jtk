use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
#[derive(Debug, Clone)]
pub struct RepeatMaskConfig {
    k: usize,
    freq: f64,
    min: u32,
}

impl RepeatMaskConfig {
    pub fn new(k: usize, freq: f64, min: u32) -> Self {
        Self { k, freq, min }
    }
}

pub trait RepeatMask {
    fn mask_repeat(&mut self, config: &RepeatMaskConfig);
    fn get_repetitive_kmer(&self) -> RepeatAnnot;
}

#[derive(Debug, Clone)]
pub struct RepeatAnnot {
    k: usize,
    kmers: HashSet<u64>,
}

impl RepeatAnnot {
    pub fn repetitiveness(&self, seq: &[u8]) -> f64 {
        let kmers = seq.len().saturating_sub(self.k) + 1;
        let rep_kmers = seq
            .windows(self.k)
            .map(to_idx)
            .filter(|x| self.kmers.contains(x))
            .count();
        rep_kmers as f64 / kmers as f64
    }
}

impl RepeatMask for definitions::DataSet {
    fn get_repetitive_kmer(&self) -> RepeatAnnot {
        let k = self.masked_kmers.k;
        let kmers: HashSet<_> = self
            .raw_reads
            .par_iter()
            .flat_map(|r| {
                r.seq()
                    .windows(k)
                    .filter(|kmer| kmer.iter().all(u8::is_ascii_lowercase))
                    .map(to_idx)
                    .collect::<Vec<_>>()
            })
            .collect();
        RepeatAnnot { k, kmers }
    }
    fn mask_repeat(&mut self, config: &RepeatMaskConfig) {
        let (mask, thr) = {
            debug!("Counting {}-mers", config.k);
            let kmer_count = kmer_counting(&self.raw_reads, config);
            debug!("Counted");
            create_mask(&kmer_count, config)
        };
        self.masked_kmers.k = config.k;
        self.masked_kmers.thr = thr;
        debug!("Constructed {}-mer filter. Size:{}", config.k, mask.len());
        self.raw_reads
            .par_iter_mut()
            .for_each(|read| mask_repeats(read.seq.seq_mut(), &mask, config.k));
        let num_bases = self.raw_reads.iter().map(|r| r.seq.len()).sum::<usize>();
        let num_lower_base = self
            .raw_reads
            .iter()
            .map(|r| r.seq.iter().filter(|x| x.is_ascii_lowercase()).count())
            .sum::<usize>();
        debug!(
            "Masked {} bases out of {} bases.",
            num_lower_base, num_bases,
        );
        debug!("{} Masked {}-mers", mask.len(), config.k);
    }
}

// pub fn mask_repeat_in_seq(seq: &mut [u8], config: &RepeatMaskConfig) {
//     let mut count: HashMap<_, _> = HashMap::with_capacity(1_000_000_000);
//     assert!(config.k <= 32);
//     for w in seq.windows(config.k) {
//         let idx = to_idx(w);
//         *count.entry(idx).or_default() += 1;
//     }
//     let (mask, _) = create_mask(&count, config);
//     mask_repeats(seq, &mask, config.k);
// }

// pub fn mask_repeats_in_reads(seqs: &mut [Vec<u8>], config: &RepeatMaskConfig) {
//     let mut count: HashMap<_, _> = HashMap::with_capacity(1_000_000_000);
//     assert!(config.k <= 32);
//     let kmers = seqs
//         .into_par_iter()
//         .map(|seq| seq.windows(config.k).map(to_idx))
//         .fold(Vec::new, |mut x, y| {
//             x.extend(y);
//             x
//         })
//         .reduce(Vec::new, |mut x, mut y| {
//             x.append(&mut y);
//             x
//         });
//     for idx in kmers {
//         *count.entry(idx).or_default() += 1;
//     }
//     let mask = create_mask(&count, config);
//     seqs.par_iter_mut()
//         .for_each(|read| mask_repeats(read, &mask, config.k));
// }

use definitions::*;
fn kmer_counting(reads: &[RawRead], config: &RepeatMaskConfig) -> HashMap<u64, u32> {
    let k = config.k;
    assert!(k <= 32);
    reads
        .into_par_iter()
        .map(|read| read.seq().windows(k).map(to_idx))
        .fold(HashMap::new, |mut x, kmers| {
            for kmer in kmers {
                *x.entry(kmer).or_default() += 1;
            }
            x
        })
        .reduce(HashMap::new, |mut x, kmercounts| {
            for (kmer, count) in kmercounts {
                *x.entry(kmer).or_default() += count;
            }
            x
        })
}

fn to_idx(w: &[u8]) -> u64 {
    // Determine if this k-mer is canonical.
    let is_canonical = {
        let mut idx = 0;
        while idx < w.len() / 2
            && w[idx].to_ascii_uppercase() == w[w.len() - idx - 1].to_ascii_uppercase()
        {
            idx += 1;
        }
        w[idx] <= w[w.len() - idx - 1]
    };
    if is_canonical {
        w.iter()
            .fold(0, |cum, &x| (cum << 2) | BASE2BIT[x as usize])
    } else {
        w.iter()
            .rev()
            .fold(0, |cum, &x| (cum << 2) | BASE2BITCMP[x as usize])
    }
}

const BASE2BITCMP: [u64; 256] = base2bitcmp();
const BASE2BIT: [u64; 256] = base2bit();

const fn base2bitcmp() -> [u64; 256] {
    let mut slots = [0; 256];
    slots[b'A' as usize] = 3;
    slots[b'a' as usize] = 3;
    slots[b'C' as usize] = 2;
    slots[b'c' as usize] = 2;
    slots[b'G' as usize] = 1;
    slots[b'g' as usize] = 1;
    slots
}

const fn base2bit() -> [u64; 256] {
    let mut slots = [0; 256];
    slots[b'C' as usize] = 1;
    slots[b'c' as usize] = 1;
    slots[b'G' as usize] = 2;
    slots[b'g' as usize] = 2;
    slots[b'T' as usize] = 3;
    slots[b't' as usize] = 3;
    slots
}

fn create_mask(kmercount: &HashMap<u64, u32>, config: &RepeatMaskConfig) -> (HashSet<u64>, u32) {
    // if the size of the hashmap is too large, scale /100.
    let mut counts: Vec<_> = if kmercount.len() > 1_000_000 {
        kmercount
            .values()
            .enumerate()
            .filter(|(i, _)| i % 100 == 0)
            .map(|(_, &c)| c)
            .collect()
    } else {
        kmercount.values().copied().collect()
    };
    counts.sort_unstable();
    counts.reverse();
    let thr = counts[((counts.len() as f64) * config.freq).floor() as usize];
    let thr = thr.max(config.min);
    debug!("Masking {}-mer occuring more than {} times", config.k, thr);
    let kmers: HashSet<_> = kmercount
        .iter()
        .filter_map(|(&key, &val)| if val > thr { Some(key) } else { None })
        .collect();
    (kmers, thr)
}

fn mask_repeats(seq: &mut [u8], mask: &HashSet<u64>, k: usize) {
    if seq.len() <= k {
        return;
    }
    let mut farthest = 0;
    for idx in 0..seq.len() - k + 1 {
        let kmer_position = to_idx(&seq[idx..idx + k]);
        if mask.contains(&kmer_position) {
            for seq in seq.iter_mut().take(idx + k).skip(farthest.max(idx)) {
                seq.make_ascii_lowercase();
            }
            farthest = idx + k - 1;
        }
    }
}
