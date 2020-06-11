pub trait Entry {
    fn new(seqs: Vec<bio_utils::fasta::Record>) -> Self;
}

impl Entry for definitions::DataSet {
    fn new(seqs: Vec<bio_utils::fasta::Record>) -> Self {
        let raw_reads: Vec<_> = seqs
            .into_iter()
            .enumerate()
            .map(|(idx, read)| definitions::RawRead {
                name: read.id().to_string(),
                desc: read.desc().unwrap_or(&String::new()).clone(),
                id: idx as u64,
                seq: String::from_utf8_lossy(read.seq()).to_string(),
            })
            .collect();
        definitions::DataSet::with_param(raw_reads, vec![], vec![], vec![], vec![], vec![])
    }
}
