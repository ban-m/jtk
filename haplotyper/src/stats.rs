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
            use std::collections::HashSet;
            let reads: HashSet<_> = self.encoded_reads.iter().map(|r| r.id).collect();
            let gap_read = self
                .raw_reads
                .iter()
                .filter(|e| !reads.contains(&e.id))
                .inspect(|e| writeln!(&mut wtr, "Gaped:{}", e.name).unwrap())
                .count();
            let gap_mean = self
                .raw_reads
                .iter()
                .filter(|e| !reads.contains(&e.id))
                .map(|e| e.seq().len())
                .sum::<usize>()
                / gap_read;
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
            writeln!(&mut wtr, "Gappy read:{}\nGapMean:{}", gap_read, gap_mean)?;
            writeln!(&mut wtr, "EncodedRate:{:.4}(%)", cover_rate)?;
        }
        // Unit statistics
        if !self.encoded_reads.is_empty() {
            let units: Vec<(u64, usize)> = {
                use std::collections::HashMap;
                let mut count: HashMap<u64, usize> = HashMap::new();
                for read in self.encoded_reads.iter() {
                    for node in read.nodes.iter() {
                        *count.entry(node.unit).or_default() += 1;
                    }
                }
                let mut count: Vec<(u64, usize)> = count.into_iter().collect();
                count.sort_by_key(|e| e.0);
                count
            };
            let max = units.iter().map(|e| e.1).max().unwrap();
            let min = units.iter().map(|e| e.1).max().unwrap();
            let sum = units.iter().map(|e| e.1).sum::<usize>();
            let ave = sum as f64 / units.len() as f64;
            writeln!(&mut wtr, "Encoding summary")?;
            writeln!(&mut wtr, "Min:{}\tMax:{}\tAve:{:.2}", max, min, ave)?;
            let units: Vec<usize> = self
                .encoded_reads
                .iter()
                .flat_map(|read| read.nodes.iter().map(|node| node.unit as usize))
                .collect();
            let hist = histgram_viz::Histgram::new(&units);
            writeln!(&mut wtr, "Unit Histgram\n{}", hist.format(20, 20))?;
        }
        Ok(())
    }
}
