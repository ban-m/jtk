pub trait Stats {
    fn stats<W: std::io::Write>(&self, wtr: W) -> std::io::Result<()>;
}

impl Stats for definitions::DataSet {
    fn stats<W: std::io::Write>(&self, mut wtr: W) -> std::io::Result<()> {
        // raw reads.
        if !self.raw_reads.is_empty() {
            let lens = self.raw_reads.iter().map(|r| r.seq().len());
            let sum = lens.clone().sum::<usize>();
            let min = lens.clone().min().unwrap_or(0);
            let max = lens.clone().max().unwrap_or(0);
            let len = self.raw_reads.len();
            let ave = sum / len;
            writeln!(&mut wtr, "Raw Reads")?;
            writeln!(
                &mut wtr,
                "Total Length:{}\n# of Read:{}\nMean Length:{}",
                sum, len, ave,
            )?;
            writeln!(&mut wtr, "Max Length:{}\nMin Length:{}", max, min)?;
        }
        // hic pairs
        if !self.hic_pairs.is_empty() {
            let lens = self
                .hic_pairs
                .iter()
                .map(|r| r.seq1().len() + r.seq2().len());
            let sum = lens.clone().sum::<usize>();
            let len = self.raw_reads.len() * 2;
            let ave = sum / len;
            writeln!(&mut wtr, "HiC Reads")?;
            writeln!(
                &mut wtr,
                "Total Length:{}\n# of Pairs:{}\nMean Length:{}",
                sum, len, ave
            )?;
        }
        // Selected chunks
        if !self.selected_chunks.is_empty() {
            let lens = self.selected_chunks.iter().map(|u| u.seq().len());
            let sum = lens.clone().sum::<usize>();
            let len = self.selected_chunks.len();
            let ave = sum / len;
            writeln!(&mut wtr, "Chunks")?;
            writeln!(
                &mut wtr,
                "Total Length:{}\n# of Units:{}\nMean Length:{}",
                sum, len, ave
            )?;
        }
        // Encoded Reads
        if !self.encoded_reads.is_empty() {
            let gap_read = self.encoded_reads.iter().filter(|e| e.is_gappy()).count();
            let covered_length = self
                .encoded_reads
                .iter()
                .map(|e| e.encoded_length())
                .sum::<usize>();
            let total_length = self
                .encoded_reads
                .iter()
                .map(|e| e.original_length)
                .sum::<usize>();
            let cover_rate = covered_length as f64 / total_length as f64;
            writeln!(&mut wtr, "EncodedRead")?;
            writeln!(
                &mut wtr,
                "Gappy read:{}\nEncodedRate:{:.4}(%)",
                gap_read, cover_rate
            )?;
        }
        Ok(())
    }
}
