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
            let lens: Vec<_> = lens.collect();
            let hist = histgram_viz::Histgram::new(&lens);
            writeln!(&mut wtr, "{}", hist.format(20, 40))?;
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
                //.inspect(|e| writeln!(&mut wtr, "Gaped:{}", e.name).unwrap())
                .count();
            let gap_sum = self
                .raw_reads
                .iter()
                .filter(|e| !reads.contains(&e.id))
                .map(|e| e.seq().len())
                .sum::<usize>();
            let gap_mean = match gap_read {
                0 => 0,
                gap_read => gap_sum / gap_read,
            };
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
            let lens: Vec<_> = self.encoded_reads.iter().map(|e| e.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&lens);
            writeln!(&mut wtr, "{}", hist.format(20, 40))?;
        }
        // Unit statistics
        if !self.encoded_reads.is_empty() {
            let mut units: Vec<(u64, usize)> = {
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
            units.sort_by_key(|e| e.1);
            let (argmax, max) = *units.last().unwrap_or(&(0, 0));
            let (argmin, min) = *units.first().unwrap_or(&(0, 0));
            let sum = units.iter().map(|e| e.1).sum::<usize>();
            let ave = sum as f64 / units.len() as f64;
            writeln!(&mut wtr, "Encoding summary")?;
            writeln!(
                &mut wtr,
                "Min:({},{})\tMax:({},{})\tAve:{:.2}",
                min, argmin, max, argmax, ave
            )?;
            let top_20: Vec<_> = units.iter().rev().take(20).copied().collect();
            let take_len = units.len() - 20.min(units.len());
            let units: Vec<_> = units.iter().take(take_len).map(|x| x.1).collect();
            let hist = histgram_viz::Histgram::new(&units);
            writeln!(&mut wtr, "Top 20 Occurences:{:?}", top_20)?;
            writeln!(&mut wtr, "The rest of the Units\n{}", hist.format(40, 20))?;
        }
        Ok(())
    }
}
