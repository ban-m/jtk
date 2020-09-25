pub trait Entry {
    fn new(input_file: &str, raw_data: &[bio_utils::fasta::Record]) -> Self;
}

impl Entry for definitions::DataSet {
    fn new(input_file: &str, raw_data: &[bio_utils::fasta::Record]) -> Self {
        let raw_reads: Vec<_> = raw_data
            .iter()
            .enumerate()
            .map(|(idx, read)| definitions::RawRead {
                name: read.id().to_string(),
                desc: read.desc().unwrap_or(&String::new()).clone(),
                id: idx as u64,
                seq: read
                    .seq()
                    .iter()
                    .map(|&x| (x as char).to_ascii_uppercase())
                    .collect::<String>(),
            })
            .collect();
        use definitions::DataSet;
        DataSet::with_minimum_data(input_file, raw_reads)
    }
}
