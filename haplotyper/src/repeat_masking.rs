use log::*;
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

#[derive(Debug, Clone)]
pub struct KMers<'a> {
    input: &'a [u8],
    idx: usize,
    forward: u64,
    reverse: u64,
    next: Option<u64>,
    k: usize,
}

impl<'a> KMers<'a> {
    fn new(input: &'a [u8], k: usize) -> Self {
        assert!(k <= 32);
        let forward: u64 = input
            .iter()
            .take(k)
            .enumerate()
            .map(|(i, &b)| BASE2BIT[b as usize] << (2 * i))
            .sum();
        let reverse = input
            .iter()
            .take(k)
            .rev()
            .enumerate()
            .map(|(i, &b)| BASE2BITCMP[b as usize] << (2 * i))
            .sum();
        let next = (k <= input.len()).then_some(forward.min(reverse));
        Self {
            input,
            idx: k,
            forward,
            reverse,
            next,
            k,
        }
    }
}

impl<'a> std::iter::Iterator for KMers<'a> {
    type Item = u64;
    fn next(&mut self) -> Option<Self::Item> {
        let ret_value = self.next;
        // Set up next value, if possible.
        if let Some(&base) = self.input.get(self.idx) {
            self.idx += 1;
            self.forward = (self.forward >> 2) | (BASE2BIT[base as usize] << (2 * (self.k - 1)));
            self.reverse = {
                let mask = (1 << (2 * self.k)) - 1;
                (mask & (self.reverse << 2)) | BASE2BITCMP[base as usize]
            };
            self.next = Some(self.forward.min(self.reverse));
        } else {
            self.next = None;
        }
        ret_value
    }
}

impl RepeatAnnot {
    /// Return the fraction of the repetitive kmer.
    /// Here, the definition of the repetitiveness is both global (is it lowercase?) and local (is it occurred more than twice in this reads?)
    pub fn repetitiveness(&self, seq: &[u8]) -> f64 {
        let mut counts: HashMap<_, u32> = HashMap::new();
        for kmer in KMers::new(seq, self.k).filter(|x| self.kmers.contains(x)) {
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
        let num_lower_base: usize = self
            .raw_reads
            .par_iter_mut()
            .map(|read| mask_repeats(read.seq.seq_mut(), &mask, config.k))
            .sum();
        let num_bases = self.raw_reads.iter().map(|r| r.seq.len()).sum::<usize>();
        debug!("MASKREPEAT\tTotalMask\t{num_lower_base}\t{num_bases}");
        debug!("MASKREPEAT\tMaskedKmers\t{}\t{}", mask.len(), config.k);
    }
}

use definitions::*;
pub fn kmer_counting(reads: &[RawRead], k: usize) -> HashMap<u64, u32> {
    assert!(k <= 32);
    reads
        .par_iter()
        .fold(HashMap::new, |mut counts, read| {
            for kmer in KMers::new(read.seq(), k) {
                *counts.entry(kmer).or_default() += 1;
            }
            counts
        })
        .reduce(HashMap::new, |mut x, y| {
            for (key, val) in y {
                *x.entry(key).or_default() += val;
            }
            x
        })
}

pub fn to_idx(w: &[u8]) -> u64 {
    let forward: u64 = w
        .iter()
        .enumerate()
        .map(|(i, &b)| BASE2BIT[b as usize] << (2 * i))
        .sum();
    let reverse = w
        .iter()
        .rev()
        .enumerate()
        .map(|(i, &b)| BASE2BITCMP[b as usize] << (2 * i))
        .sum();
    forward.min(reverse)
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

fn mask_repeats(seq: &mut [u8], mask: &HashSet<u64>, k: usize) -> usize {
    let (mut start, mut end) = (0, 0);
    let mut mask_ranges = vec![];
    for (idx, kmer) in KMers::new(seq, k).enumerate() {
        if mask.contains(&kmer) {
            if end < idx {
                mask_ranges.push((start, end));
                start = idx;
                end = idx + k;
            } else {
                end = idx + k;
            }
        }
    }
    mask_ranges.push((start, end));
    let num_lower = mask_ranges.iter().map(|(s, e)| e - s).sum();
    for (start, end) in mask_ranges {
        seq.iter_mut()
            .take(end)
            .skip(start)
            .for_each(|s| s.make_ascii_lowercase());
    }
    num_lower
}

#[cfg(test)]
pub mod test {
    use super::*;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    use rand_xoshiro::Xoroshiro128Plus;
    #[test]
    fn kmers() {
        let kmers = b"AAAA";
        let mut kmers_iter = KMers::new(kmers, 1);
        assert_eq!(kmers_iter.clone().count(), 4);
        assert!(kmers_iter.all(|x| x == 0));
        assert!(KMers::new(kmers, 2).all(|x| x == 0));
        let kmers = b"TTTT";
        let mut kmers_iter = KMers::new(kmers, 1);
        assert_eq!(kmers_iter.clone().count(), 4);
        assert!(kmers_iter.all(|x| x == 0));
        assert!(KMers::new(kmers, 2).all(|x| x == 0));

        let kmers = b"CAGTGCAT";
        let k = 2;
        let kmers_iter = KMers::new(kmers, k);
        println!("INIT\t{}\t{}", kmers_iter.forward, kmers_iter.reverse);
        for (i, kmer) in kmers_iter.enumerate() {
            assert!(kmer < (1 << (2 * k)));
            let input = &kmers[i..i + k];
            let seq = String::from_utf8_lossy(input);
            assert_eq!(kmer, to_idx(input), "{},{:0b},{}", i, kmer, seq);
        }

        let mut rng: Xoroshiro128Plus = SeedableRng::seed_from_u64(832904);
        let len = 1000;
        let kmers: Vec<u8> = (0..len)
            .map(|_| *b"ACGT".choose(&mut rng).unwrap())
            .collect();
        let k = 10;
        let kmers_iter = KMers::new(&kmers, k);
        assert_eq!(kmers_iter.clone().count(), len - k + 1);
        for (i, kmer) in kmers_iter.enumerate() {
            assert!(kmer < (1 << (2 * k)));
            let input = &kmers[i..i + k];
            let seq = String::from_utf8_lossy(input);
            assert_eq!(kmer, to_idx(input), "{},{:0b},{}", i, kmer, seq);
        }
    }
}
