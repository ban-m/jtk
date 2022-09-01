//! This module includes light statistics metrics, such as error rates and coverages...
use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::collections::HashSet;

pub trait Stats {
    fn error_rate(&self) -> ErrorRate;
    fn stats<W: std::io::Write>(&self, wtr: W) -> std::io::Result<()>;
}

impl Stats for definitions::DataSet {
    fn stats<W: std::io::Write>(&self, mut wtr: W) -> std::io::Result<()> {
        // raw reads.
        if !self.raw_reads.is_empty() {
            let lens: Vec<_> = self.raw_reads.iter().map(|r| r.seq().len()).collect();
            let sum = lens.iter().sum::<usize>();
            let min = lens.iter().min().unwrap_or(&0);
            let max = lens.iter().max().unwrap_or(&0);
            let len = self.raw_reads.len();
            let ave = sum / len;
            let nfif = up_to(&lens, 0.5);
            let masked: usize = self
                .raw_reads
                .iter()
                .map(|r| r.seq().iter().filter(|&x| x.is_ascii_lowercase()).count())
                .sum();
            writeln!(wtr, "RAWREADS\tTotalLength\t{sum}")?;
            writeln!(wtr, "RAWREADS\tMaskedLength\t{masked}")?;
            writeln!(wtr, "RAWREADS\tNumOfRead\t{len}")?;
            writeln!(wtr, "RAWREADS\tMeanLength\t{ave}")?;
            writeln!(wtr, "RAWREADS\tMaxLength\t{max}")?;
            writeln!(wtr, "RAWREADS\tMinLength\t{min}")?;
            writeln!(wtr, "RAWREADS\tN50\t{nfif}")?;
            let mut lens: Vec<_> = self.raw_reads.iter().map(|r| r.seq().len()).collect();
            lens.sort_unstable();
            lens.reverse();
            let top_20: Vec<_> = lens
                .iter()
                .take(20)
                .map(|&x| format!("{:.1}K", x as f64 / 1_000f64))
                .collect();
            let lens = &lens[lens.len().min(20)..];
            let hist = histgram_viz::Histgram::new(lens);
            writeln!(wtr, "Top 20 Occurences:{}", top_20.join("\t"))?;
            writeln!(wtr, "{}", hist.format(20, 40))?;
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
            writeln!(wtr, "HICREADS\tNumOfReads\t{sum}")?;
            writeln!(wtr, "HICREADS\tTotalLength\t{len}")?;
            writeln!(wtr, "HICREADS\tMeanLength\t{ave}")?;
        }
        // Selected chunks
        if !self.selected_chunks.is_empty() {
            let lens = self.selected_chunks.iter().map(|u| u.seq().len());
            let sum = lens.clone().sum::<usize>();
            let len = self.selected_chunks.len();
            let ave = sum / len;
            writeln!(wtr, "CHUNKS\tTotalLength\t{sum}")?;
            writeln!(wtr, "CHUNKS\tNumOfUnits\t{len}")?;
            writeln!(wtr, "CHUNKS\tMeanLength\t{ave}")?;
        }
        // Encoded Reads
        if !self.encoded_reads.is_empty() {
            let reads: HashSet<_> = self
                .encoded_reads
                .iter()
                .filter(|r| !r.is_gappy())
                .map(|r| r.id)
                .collect();
            let gap_read = self
                .raw_reads
                .iter()
                .filter(|e| !reads.contains(&e.id))
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
            let num_nodes: usize = self.encoded_reads.iter().map(|r| r.nodes.len()).sum();
            let enc_num = self.encoded_reads.len();
            writeln!(wtr, "EncodedRead\tNumEncoded\t{enc_num}",)?;
            let enc_len = self
                .encoded_reads
                .iter()
                .flat_map(|r| r.nodes.iter())
                .map(|n| n.query_length())
                .sum::<usize>();
            writeln!(wtr, "ENCODEDREAD\tLenEncoded\t{enc_len}",)?;
            writeln!(wtr, "ENCODEDREAD\tNumGappy\t{gap_read}")?;
            writeln!(wtr, "ENCODEDREAD\tGapMean\t{gap_mean}")?;
            writeln!(wtr, "ENCODEDRATE\t{:.4}%", cover_rate)?;
            writeln!(wtr, "ENCODEDNODE\t{:.4}", num_nodes)?;
            let lens: Vec<_> = self.encoded_reads.iter().map(|e| e.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&lens);
            writeln!(wtr, "{}", hist.format(20, 40))?;
        }
        // Unit statistics
        if !self.encoded_reads.is_empty() {
            let mut units: Vec<(u64, usize)> = {
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
            writeln!(wtr, "ENCODING\tMin\t{min}\t{argmin}")?;
            writeln!(wtr, "ENCODING\tMax\t{max}\t{argmax}")?;
            writeln!(wtr, "ENCODING\tAve\t{ave:.2}")?;
            let top_20: Vec<_> = units.iter().rev().take(20).copied().collect();
            let take_len = units.len() - 20.min(units.len());
            let units: Vec<_> = units.iter().take(take_len).map(|x| x.1).collect();
            let hist = histgram_viz::Histgram::new(&units);
            writeln!(wtr, "Top 20 Occurences:{:?}", top_20)?;
            writeln!(wtr, "The rest of the Units\n{}", hist.format(40, 20))?;
        }
        // Encoding errors
        if !self.encoded_reads.is_empty() {
            let error = self.error_rate();
            writeln!(wtr, "ErrorRate\n{}", error)?;
        }
        Ok(())
    }
    fn error_rate(&self) -> ErrorRate {
        let ref_units: HashMap<_, _> = self.selected_chunks.iter().map(|c| (c.id, c)).collect();
        let summaries: Vec<_> = self
            .encoded_reads
            .par_iter()
            .flat_map(|r| r.nodes.par_iter())
            .map(|node| {
                let ref_unit = ref_units[&node.unit];
                let (query, aln, refr) = node.recover(ref_unit);
                let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
                let del = query.iter().filter(|&&x| x == b' ').count() as f64;
                let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
                let aln_len = aln.len() as f64;
                (del / aln_len, ins / aln_len, mismat / aln_len)
            })
            .collect();
        debug!("Summary:{}", summaries.len());
        let del_summary = summarize(summaries.iter().map(|x| x.0));
        let ins_summary = summarize(summaries.iter().map(|x| x.1));
        let mism_summary = summarize(summaries.iter().map(|x| x.2));
        let total_summary = summarize(summaries.iter().map(|(x, y, z)| x + y + z));
        ErrorRate::new(del_summary, ins_summary, mism_summary, total_summary)
    }
}

//Return median and median absolute deviation.
fn summarize<I: std::iter::Iterator<Item = f64>>(error_rates: I) -> (f64, f64) {
    let mut errors: Vec<_> = error_rates.collect();
    let idx = errors.len() / 2;
    let median = *errors
        .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap())
        .1;
    // Convert to abs deviation.
    errors.iter_mut().for_each(|x| *x = (*x - median).abs());
    let mad = *errors
        .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap())
        .1;
    (median, mad)
}

fn up_to(lens: &[usize], frac: f64) -> usize {
    let mut lens = lens.to_vec();
    lens.sort_unstable();
    lens.reverse();
    let sum: usize = lens.iter().sum();
    let target = (frac * sum as f64).floor() as usize;
    let mut acc = 0;
    for len in lens.iter() {
        acc += len;
        if target <= acc {
            return *len;
        }
    }
    panic!()
}
