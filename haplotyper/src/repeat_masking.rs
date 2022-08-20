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
    /// Return the fraction of the repetitive kmer.
    /// Here, the definition of the repetitiveness is both global (is it lowercase?) and local (is it occurred more than twice in this reads?)
    pub fn repetitiveness(&self, seq: &[u8]) -> f64 {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for kmer in seq
            .windows(self.k)
            .map(to_idx)
            .filter(|x| self.kmers.contains(x))
        {
            *counts.entry(kmer).or_default() += 1;
        }
        let rep_kmers: u32 = counts.values().filter(|&&count| 1 < count).sum();
        let kmers = seq.len().saturating_sub(self.k) + 1;
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
        self.masked_kmers.k = config.k;
        self.raw_reads.par_iter_mut().for_each(|r| {
            r.seq
                .seq_mut()
                .iter_mut()
                .for_each(u8::make_ascii_uppercase)
        });
        let (mask, _thr) = create_mask(&self.raw_reads, config);
        debug!("MASKREPEAT\tMaskLen\t{}\t{}", config.k, mask.len());
        self.raw_reads
            .par_iter_mut()
            .for_each(|read| mask_repeats(read.seq.seq_mut(), &mask, config.k));
        let num_bases = self.raw_reads.iter().map(|r| r.seq.len()).sum::<usize>();
        let num_lower_base = self
            .raw_reads
            .iter()
            .map(|r| r.seq.iter().filter(|x| x.is_ascii_lowercase()).count())
            .sum::<usize>();
        debug!("MASKREPEAT\tTotalMask\t{num_lower_base}\t{num_bases}");
        debug!("MASKREPEAT\tMaskedKmers\t{}\t{}", mask.len(), config.k);
    }
}

const BUCKET_SIZE: usize = 1000;
use definitions::*;
pub fn kmer_counting(reads: &[RawRead], k: usize) -> HashMap<u64, u32> {
    assert!(k <= 32);
    let mut counts: HashMap<u64, u32> = HashMap::new();
    for bucket in reads.chunks(BUCKET_SIZE) {
        let kmers: Vec<_> = bucket
            .into_par_iter()
            .fold(Vec::new, |mut x, read| {
                x.extend(read.seq().windows(k).map(to_idx));
                x
            })
            .flatten()
            .collect();
        for kmer in kmers {
            *counts.entry(kmer).or_default() += 1;
        }
    }
    counts
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

fn create_mask(reads: &[RawRead], config: &RepeatMaskConfig) -> (HashSet<u64>, u32) {
    let mut kmercount = kmer_counting(reads, config.k);
    let (mut num_singleton, total) = (0, kmercount.len());
    kmercount.retain(|_, &mut count| {
        num_singleton += (count == 1) as usize;
        count != 1
    });
    debug!("MASKREPEAT\tUNIQUE\t{num_singleton}\t{total}");
    let percentile = (kmercount.len() as f64 * (1f64 - config.freq)).floor() as usize;
    let below_min = kmercount.values().filter(|&&x| x <= config.min).count();
    if percentile < below_min {
        warn!("MASKREPEAT\tTakesAll\t{}\t{}", config.min, config.freq);
        let kmers: HashSet<_> = kmercount.keys().copied().collect();
        (kmers, 0)
    } else {
        let percentile = percentile - below_min;
        let mut counts: Vec<_> = kmercount
            .values()
            .filter(|&&count| config.min < count)
            .copied()
            .collect();
        counts.sort_unstable();
        let thr = counts[percentile];
        debug!("MASKREPEAT\tMaskMoreThan\t{thr}\t{}", config.k);
        let kmers: HashSet<_> = kmercount
            .iter()
            .filter_map(|(&key, &val)| if val > thr { Some(key) } else { None })
            .collect();
        (kmers, thr)
    }
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
