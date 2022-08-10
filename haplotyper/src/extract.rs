pub trait Extract {
    fn extract<W: std::io::Write>(&self, file: &mut W) -> std::io::Result<()>;
}

impl Extract for definitions::DataSet {
    fn extract<W: std::io::Write>(&self, file: &mut W) -> std::io::Result<()> {
        for read in self.raw_reads.iter() {
            let seq = std::str::from_utf8(read.seq()).unwrap();
            let (name, desc, id) = (&read.name, &read.desc, read.id);
            writeln!(file, "READ\t{name}\t{desc}\t{id}\t{seq}")?;
        }
        for unit in self.selected_chunks.iter() {
            let seq = std::str::from_utf8(unit.seq()).unwrap();
            writeln!(file, "UNIT\t{}\t{}\t{}", unit.id, unit.copy_num, seq)?;
        }
        for _read in self.encoded_reads.iter() {}
        Ok(())
    }
}
