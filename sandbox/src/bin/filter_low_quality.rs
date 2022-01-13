// calc k-mer entropy on each read.
extern crate bio_utils;
extern crate rayon;
use rayon::prelude::*;
use std::ops::Shl;
const K: usize = 6;
const THR: f64 = 7.0;
const LEN_THR: usize = 4_000;
fn main() -> std::io::Result<()> {
    let stdin = std::io::stdin();
    let reads = bio_utils::fasta::parse_into_vec_from(stdin.lock())?;
    let stdout = std::io::stdout();
    let mut stdout = bio_utils::fasta::Writer::new(stdout.lock());
    let reads: Vec<_> = reads
        .into_par_iter()
        .filter(|read| calc_entropy(read.seq(), K) > THR)
        .filter(|read| read.seq().len() > LEN_THR)
        .collect();
    eprintln!("Records:{}", reads.len());
    for read in reads {
        stdout.write_record(&read)?;
    }
    Ok(())
}

fn calc_entropy(read: &[u8], k: usize) -> f64 {
    if read.len() < k {
        0.
    } else {
        let mut slots: Vec<u32> = vec![0; 4usize.pow(k as u32)];
        let mask = 4usize.pow(k as u32) - 1;
        let mut current = calc_index(&read[..k - 1]);
        let total = (read.len() - k + 1) as f64;
        for base in &read[k..] {
            current = current.shl(2) & mask;
            current += match base {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => unreachable!(),
            };
            slots[current] += 1;
        }
        slots
            .into_iter()
            .filter(|&count| count != 0)
            .map(|e| e as f64 / total)
            .map(|e| -e * e.log2())
            .sum::<f64>()
    }
}

#[inline]
fn calc_index(seq: &[u8]) -> usize {
    seq.iter()
        .map(|base| match base {
            b'A' | b'a' => 0usize,
            b'C' | b'c' => 1usize,
            b'G' | b'g' => 2usize,
            b'T' | b't' => 3usize,
            _ => unreachable!(),
        })
        .fold(0, |sum, b| sum.shl(2) + b)
}
