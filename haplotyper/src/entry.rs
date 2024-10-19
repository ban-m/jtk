use definitions::ReadType;
use std::path::Path;
pub trait Entry {
    fn entry(input_file: &Path, raw_data: Vec<(String, Vec<u8>)>, read_type: ReadType) -> Self;
}

impl Entry for definitions::DataSet {
    fn entry(input_file: &Path, raw_data: Vec<(String, Vec<u8>)>, read_type: ReadType) -> Self {
        let compress_thr = match read_type {
            ReadType::CCS => 100,
            ReadType::CLR => 40,
            ReadType::ONT => 100,
            ReadType::None => 100,
        };
        let raw_reads: Vec<_> = raw_data
            .into_iter()
            .enumerate()
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
        DataSet::new(input_file, raw_reads, read_type)
    }
}

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
