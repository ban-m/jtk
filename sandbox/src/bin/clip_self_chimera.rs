extern crate bio;
extern crate rayon;
use bio::alignment::pairwise;
use bio::io::fasta;
use rayon::prelude::*;
const SCORE_THR: i32 = 80;
const MAP_DIFF_THR: isize = 100;
fn main() -> std::io::Result<()> {
    let start = std::time::Instant::now();
    let stdin = std::io::stdin();
    let reads: Vec<_> = fasta::Reader::new(stdin.lock())
        .records()
        .filter_map(|e| e.ok())
        .collect();
    let read_fin = std::time::Instant::now();
    let result: Vec<_> = reads.into_par_iter().map(trim_self_chimera).collect();
    let end = std::time::Instant::now();
    let stdout = std::io::stdout();
    let mut wtr = fasta::Writer::new(stdout.lock());
    result
        .into_iter()
        .for_each(|read| wtr.write_record(&read).unwrap());
    let write_end = std::time::Instant::now();
    eprintln!("Read Fin: {:?}", read_fin - start);
    eprintln!("Convert: {:?}", end - read_fin);
    eprintln!("Output: {:?}", write_end - end);
    eprintln!("Total: {:?}", write_end - start);
    Ok(())
}

#[inline]
fn trim_self_chimera(read: fasta::Record) -> fasta::Record {
    let seq = read.seq();
    let seq2 = bio::alphabets::dna::revcomp(seq);
    let w = seq.len() * seq.len() / 12_000_000;
    let mut aligner =
        pairwise::banded::Aligner::new(-4, -1, |a, b| if a == b { 1 } else { -1 }, 5, w);
    let result = aligner.local(seq, &seq2);
    if let Some((start, end)) = determine(&result) {
        let seq = &read.seq()[start..end];
        fasta::Record::with_attrs(read.id(), read.desc(), seq)
    } else {
        read
    }
}

#[inline]
fn determine(align: &bio::alignment::Alignment) -> Option<(usize, usize)> {
    if align.score < SCORE_THR {
        return None;
    }
    let len = align.xlen;
    let tstart = align.xstart;
    let tend = align.xend;
    // Y is revcomp of X.
    let rstart = align.ylen - align.yend;
    let rend = align.ylen - align.ystart;
    if (tstart as isize - rstart as isize).abs() < MAP_DIFF_THR
        && (tend as isize - rend as isize) < MAP_DIFF_THR
    {
        // Usual self chimeta
        let start = (tstart + rstart) / 2;
        let end = (tend + rend) / 2;
        let break_point = (start + end) / 2;
        if start < (len - end) {
            return Some((break_point, len));
        } else {
            return Some((0, break_point));
        }
    };
    // At this closure, the alignment seems to be splitted
    // into two part of the read.
    // First, we check wheter or not the two aligned region interconnected.
    // To this end, it is suffice to compare the bigger start position
    // and the smaller end position.
    let bigger_start = tstart.max(rstart);
    let smaller_end = tend.min(rend);
    if smaller_end < bigger_start {
        // The two section do not overlap. Thus, the read is considered as
        // split chimera.
        let start = tstart.min(rstart);
        let end = tend.max(rstart);
        if start < len - end {
            return Some((bigger_start, len));
        } else {
            return Some((0, smaller_end));
        }
    }
    None
}
