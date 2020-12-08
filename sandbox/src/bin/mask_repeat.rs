const K: usize = 20;
const THR: u32 = 2000;
fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let reads = bio_utils::fasta::parse_into_vec(&args[1])?;
    let kmer_count = kmer_counting(&reads, K);
    let mask = create_mask(&kmer_count, THR);
    let reads: Vec<_> = reads
        .into_iter()
        .map(|read| mask_repeats(read, &mask, K))
        .collect();
    let out = std::io::stdout();
    let mut wtr = bio_utils::fasta::Writer::new(out.lock());
    for read in reads {
        wtr.write_record(&read)?;
    }
    Ok(())
}
use bio_utils::fasta::Record;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
fn kmer_counting(reads: &[Record], k: usize) -> HashMap<u64, u32> {
    assert!(k <= 32);
    let mut result = HashMap::with_capacity(1_000_000_000);
    let kmers = reads
        .into_par_iter()
        .map(|read| read.seq().windows(k).map(|w| to_idx(w)))
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
        b'A' | b'a' => (sum << 2) | 0u64,
        b'C' | b'c' => (sum << 2) | 1u64,
        b'G' | b'g' => (sum << 2) | 2u64,
        b'T' | b't' => (sum << 2) | 3u64,
        _ => (sum << 2),
    };
    let adder_rev = |sum, &c| match c {
        b'A' | b'a' => (sum << 2) | 3u64,
        b'C' | b'c' => (sum << 2) | 2u64,
        b'G' | b'g' => (sum << 2) | 1u64,
        b'T' | b't' => (sum << 2) | 0u64,
        _ => (sum << 2),
    };
    if is_canonical {
        w.iter().fold(0, adder)
    } else {
        w.iter().rev().fold(0, adder_rev)
    }
}

fn create_mask(kmercount: &HashMap<u64, u32>, thr: u32) -> HashSet<u64> {
    kmercount
        .iter()
        .inspect(|(&key, &val)| {
            if val > thr {
                eprintln!("{:b}\t{}", key, val);
            }
        })
        .filter_map(|(&key, &val)| if val > thr { Some(key) } else { None })
        .collect()
}

fn mask_repeats(read: Record, mask: &HashSet<u64>, k: usize) -> Record {
    let mut seq = vec![];
    let mut idx = 0;
    while idx <= read.seq().len() - k {
        let kmer = &read.seq()[idx..idx + k];
        let kmer_position = to_idx(kmer);
        if mask.contains(&kmer_position) {
            seq.extend(kmer.iter().map(|&x| x.to_ascii_lowercase()));
            idx += k;
        } else {
            seq.push(read.seq()[idx]);
            idx += 1;
        }
    }
    while idx < read.seq().len() {
        seq.push(read.seq()[idx]);
        idx += 1;
    }
    let desc = read.desc().map(|x| x.clone());
    Record::with_data(read.id(), &desc, &seq)
}
