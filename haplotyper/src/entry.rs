pub trait Entry {
    fn entry(input_file: &str, raw_data: &[bio_utils::fasta::Record], rt: &str) -> Self;
}

impl Entry for definitions::DataSet {
    fn entry(input_file: &str, raw_data: &[bio_utils::fasta::Record], rt: &str) -> Self {
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
            .iter()
            .enumerate()
            .map(|(idx, read)| definitions::RawRead {
                name: read.id().to_string(),
                desc: read.desc().unwrap_or(&String::new()).clone(),
                id: idx as u64,
                seq: String::from_utf8(compress_homopolymer(read.seq(), compress_thr)).unwrap(),
            })
            .collect();
        use definitions::DataSet;
        DataSet::with_minimum_data(input_file, raw_reads, read_type)
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
