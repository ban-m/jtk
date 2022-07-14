pub trait Entry {
    fn entry(input_file: &str, raw_data: Vec<(String, Vec<u8>)>, rt: &str) -> Self;
}

impl Entry for definitions::DataSet {
    fn entry(input_file: &str, raw_data: Vec<(String, Vec<u8>)>, rt: &str) -> Self {
        use definitions::ReadType;
        let read_type = match rt {
            "CLR" => ReadType::CLR,
            "CCS" => ReadType::CCS,
            "ONT" => ReadType::ONT,
            _ => ReadType::None,
        };
        let compress_thr = match read_type {
            ReadType::CCS => 100,
            ReadType::CLR => 40,
            ReadType::ONT => 100,
            ReadType::None => 100,
        };
        let raw_reads: Vec<_> = raw_data
            .into_iter()
            .enumerate()
            //             .filter(|(_, (_, seq))| seq.len() > 4_000) // TODO: Why this filtering improve accuracy?
            .map(|(idx, (name, seq))| {
                let seq: definitions::DNASeq = compress_homopolymer(&seq, compress_thr).into();
                let id = idx as u64;
                definitions::RawRead {
                    name,
                    desc: String::new(),
                    seq,
                    id,
                }
            })
            .collect();
        debug!("Input\tReadNum\t{}", raw_reads.len());
        let sum: usize = raw_reads.iter().map(|r| r.seq().len()).sum();
        debug!("Input\tBasePair(Mbp)\t{}", sum / 1_000_000);
        use definitions::DataSet;
        for read in raw_reads.iter() {
            if read.seq().iter().any(|x| !b"ACGT".contains(x)) {
                let base = read.seq().iter().find(|x| !b"ACGT".contains(x)).unwrap();
                panic!("{}", base);
            }
        }
        DataSet::with_minimum_data(input_file, raw_reads, read_type)
    }
}

// TODO:Implement this.

// const SCORE_THR: i32 = 80;
// const MAP_DIFF_THR: isize = 100;

// #[inline]
// fn trim_self_chimera(read: &mut bio_utils::fasta::Record) -> usize {
//     let seq = read.seq();
//     let seq2 = bio_utils::revcmp(seq);
//     let w = seq.len() * seq.len() / 12_000_000;
//     let mut aligner =
//         pairwise::banded::Aligner::new(-4, -1, |a, b| if a == b { 1 } else { -1 }, 5, w);
//     let result = aligner.local(seq, &seq2);
//     if let Some((start, end)) = determine(&result) {
//         let seq = &read.seq()[start..end];
//         fasta::Record::with_attrs(read.id(), read.desc(), seq)
//     } else {
//         read
//     }
// }

// #[inline]
// fn determine(align: &bio::alignment::Alignment) -> Option<(usize, usize)> {
//     if align.score < SCORE_THR {
//         return None;
//     }
//     let len = align.xlen;
//     let tstart = align.xstart;
//     let tend = align.xend;
//     // Y is revcomp of X.
//     let rstart = align.ylen - align.yend;
//     let rend = align.ylen - align.ystart;
//     if (tstart as isize - rstart as isize).abs() < MAP_DIFF_THR
//         && (tend as isize - rend as isize) < MAP_DIFF_THR
//     {
//         // Usual self chimeta
//         let start = (tstart + rstart) / 2;
//         let end = (tend + rend) / 2;
//         let break_point = (start + end) / 2;
//         if start < (len - end) {
//             return Some((break_point, len));
//         } else {
//             return Some((0, break_point));
//         }
//     };
//     // At this closure, the alignment seems to be splitted
//     // into two part of the read.
//     // First, we check wheter or not the two aligned region interconnected.
//     // To this end, it is suffice to compare the bigger start position
//     // and the smaller end position.
//     let bigger_start = tstart.max(rstart);
//     let smaller_end = tend.min(rend);
//     if smaller_end < bigger_start {
//         // The two section do not overlap. Thus, the read is considered as
//         // split chimera.
//         let start = tstart.min(rstart);
//         let end = tend.max(rstart);
//         if start < len - end {
//             return Some((bigger_start, len));
//         } else {
//             return Some((0, smaller_end));
//         }
//     }
//     None
// }

// Compress homopolymer longer than x base into x base.
fn compress_homopolymer(seq: &[u8], len: usize) -> Vec<u8> {
    let mut compressed = Vec::with_capacity(seq.len());
    let mut seq = seq.iter().map(|x| x.to_ascii_uppercase()).peekable();
    while let Some(base) = seq.next() {
        let mut length = 1;
        while Some(base) == seq.peek().cloned() {
            seq.next().unwrap();
            length += 1;
        }
        compressed.extend(std::iter::repeat(base).take(length.min(len)));
    }
    compressed
}

// use std::ops::Shl;
// fn calc_entropy(read: &[u8], k: usize) -> f64 {
//     if read.len() < k {
//         0.
//     } else {
//         let mut slots: Vec<u32> = vec![0; 4usize.pow(k as u32)];
//         let mask = 4usize.pow(k as u32) - 1;
//         let mut current = calc_index(&read[..k - 1]);
//         let total = (read.len() - k + 1) as f64;
//         for base in &read[k..] {
//             current = current.shl(2) & mask;
//             current += match base {
//                 b'A' | b'a' => 0,
//                 b'C' | b'c' => 1,
//                 b'G' | b'g' => 2,
//                 b'T' | b't' => 3,
//                 _ => unreachable!(),
//             };
//             slots[current] += 1;
//         }
//         slots
//             .into_iter()
//             .filter(|&count| count != 0)
//             .map(|e| e as f64 / total)
//             .map(|e| -e * e.log2())
//             .sum::<f64>()
//     }
// }

// #[inline]
// fn calc_index(seq: &[u8]) -> usize {
//     seq.iter()
//         .map(|base| match base {
//             b'A' | b'a' => 0usize,
//             b'C' | b'c' => 1usize,
//             b'G' | b'g' => 2usize,
//             b'T' | b't' => 3usize,
//             _ => unreachable!(),
//         })
//         .fold(0, |sum, b| sum.shl(2) + b)
// }

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn compress_test() {
        let seq = b"AAACCAAAAAAAAATTGGGCTTT";
        let compressed = compress_homopolymer(seq, 2);
        assert_eq!(compressed, b"AACCAATTGGCTT");
        let compressed = compress_homopolymer(seq, 3);
        assert_eq!(compressed, b"AAACCAAATTGGGCTTT");
    }
}
