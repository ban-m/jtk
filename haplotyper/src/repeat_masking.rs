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
}

impl RepeatMask for definitions::DataSet {
    fn mask_repeat(&mut self, config: &RepeatMaskConfig) {
        let mask = {
            debug!("Counting {}-mers", config.k);
            let kmer_count = kmer_counting(&self.raw_reads, config);
            debug!("Counted");
            create_mask(&kmer_count, config)
        };
        debug!("Constructed {}-mer filter. Size:{}", config.k, mask.len());
        self.raw_reads
            .par_iter_mut()
            .for_each(|read| mask_repeats(read.seq.as_mut(), &mask, config.k));
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

pub fn mask_repeat_in_seq(seq: &mut [u8], config: &RepeatMaskConfig) {
    let mut count: HashMap<_, _> = HashMap::with_capacity(1_000_000_000);
    assert!(config.k <= 32);
    for w in seq.windows(config.k) {
        let idx = to_idx(w);
        *count.entry(idx).or_default() += 1;
    }
    let mask = create_mask(&count, config);
    mask_repeats(seq, &mask, config.k);
}

pub fn mask_repeats_in_reads(seqs: &mut [Vec<u8>], config: &RepeatMaskConfig) {
    let mut count: HashMap<_, _> = HashMap::with_capacity(1_000_000_000);
    assert!(config.k <= 32);
    let kmers = seqs
        .into_par_iter()
        .map(|seq| seq.windows(config.k).map(to_idx))
        .fold(Vec::new, |mut x, y| {
            x.extend(y);
            x
        })
        .reduce(Vec::new, |mut x, mut y| {
            x.append(&mut y);
            x
        });
    for idx in kmers {
        *count.entry(idx).or_default() += 1;
    }
    let mask = create_mask(&count, config);
    seqs.par_iter_mut()
        .for_each(|read| mask_repeats(read, &mask, config.k));
}

use definitions::*;
fn kmer_counting(reads: &[RawRead], config: &RepeatMaskConfig) -> HashMap<u64, u32> {
    let k = config.k;
    assert!(k <= 32);
    let mut result = HashMap::with_capacity(1_000_000_000);
    let kmers = reads
        .into_par_iter()
        .map(|read| read.seq().windows(k).map(to_idx))
        .fold(Vec::new, |mut x, y| {
            x.extend(y);
            x
        })
        .reduce(Vec::new, |mut x, mut y| {
            x.append(&mut y);
            x
        });
    for idx in kmers {
        *result.entry(idx).or_default() += 1;
    }
    result
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
    let adder = |sum, &c| match c {
        b'A' | b'a' => sum << 2,
        b'C' | b'c' => (sum << 2) | 1u64,
        b'G' | b'g' => (sum << 2) | 2u64,
        b'T' | b't' => (sum << 2) | 3u64,
        _ => (sum << 2),
    };
    let adder_rev = |sum, &c| match c {
        b'A' | b'a' => (sum << 2) | 3u64,
        b'C' | b'c' => (sum << 2) | 2u64,
        b'G' | b'g' => (sum << 2) | 1u64,
        b'T' | b't' => sum << 2,
        _ => (sum << 2),
    };
    if is_canonical {
        w.iter().fold(0, adder)
    } else {
        w.iter().rev().fold(0, adder_rev)
    }
}

fn create_mask(kmercount: &HashMap<u64, u32>, config: &RepeatMaskConfig) -> HashSet<u64> {
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
    kmercount
        .iter()
        .filter_map(|(&key, &val)| if val > thr { Some(key) } else { None })
        .collect()
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
