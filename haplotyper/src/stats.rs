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
            let lens = self.raw_reads.iter().map(|r| r.seq().len());
            let sum = lens.clone().sum::<usize>();
            let min = lens.clone().min().unwrap_or(0);
            let max = lens.clone().max().unwrap_or(0);
            let len = self.raw_reads.len();
            let ave = sum / len;
            let masked: usize = self
                .raw_reads
                .iter()
                .map(|r| r.seq().iter().filter(|&x| x.is_ascii_lowercase()).count())
                .sum();
            writeln!(&mut wtr, "RAWREADS\tTotalLength\t{sum}")?;
            writeln!(&mut wtr, "RAWREADS\tMaskedLength\t{masked}")?;
            writeln!(&mut wtr, "RAWREADS\tNumOfRead\t{len}")?;
            writeln!(&mut wtr, "RAWREADS\tMeanLength\t{ave}")?;
            writeln!(&mut wtr, "RAWREADS\tMaxLength\t{max}")?;
            writeln!(&mut wtr, "RAWREADS\tMinLength\t{min}")?;
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
            writeln!(&mut wtr, "Top 20 Occurences:{}", top_20.join("\t"))?;
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
            writeln!(&mut wtr, "HICREADS\tNumOfReads\t{sum}")?;
            writeln!(&mut wtr, "HICREADS\tTotalLength\t{len}")?;
            writeln!(&mut wtr, "HICREADS\tMeanLength\t{ave}")?;
        }
        // Selected chunks
        if !self.selected_chunks.is_empty() {
            let lens = self.selected_chunks.iter().map(|u| u.seq().len());
            let sum = lens.clone().sum::<usize>();
            let len = self.selected_chunks.len();
            let ave = sum / len;
            writeln!(&mut wtr, "CHUNKS\tTotalLength\t{sum}")?;
            writeln!(&mut wtr, "CHUNKS\tNumOfUnits\t{len}")?;
            writeln!(&mut wtr, "CHUNKS\tMeanLength\t{ave}")?;
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
            writeln!(&mut wtr, "EncodedRead\tNumEncoded\t{enc_num}",)?;
            let enc_len = self
                .encoded_reads
                .iter()
                .flat_map(|r| r.nodes.iter())
                .map(|n| n.query_length())
                .sum::<usize>();
            writeln!(&mut wtr, "ENCODEDREAD\tLenEncoded\t{enc_len}",)?;
            writeln!(&mut wtr, "ENCODEDREAD\tNumGappy\t{gap_read}")?;
            writeln!(&mut wtr, "ENCODEDREAD\tGapMean\t{gap_mean}")?;
            writeln!(&mut wtr, "ENCODEDRATE\t{:.4}%", cover_rate)?;
            writeln!(&mut wtr, "ENCODEDNODE\t{:.4}", num_nodes)?;
            let lens: Vec<_> = self.encoded_reads.iter().map(|e| e.nodes.len()).collect();
            let hist = histgram_viz::Histgram::new(&lens);
            writeln!(&mut wtr, "{}", hist.format(20, 40))?;
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
            writeln!(&mut wtr, "ENCODING\tMin\t{min}\t{argmin}")?;
            writeln!(&mut wtr, "ENCODING\tMax\t{max}\t{argmax}")?;
            writeln!(&mut wtr, "ENCODING\tAve\t{ave:.2}")?;
            let top_20: Vec<_> = units.iter().rev().take(20).copied().collect();
            let take_len = units.len() - 20.min(units.len());
            let units: Vec<_> = units.iter().take(take_len).map(|x| x.1).collect();
            let hist = histgram_viz::Histgram::new(&units);
            writeln!(&mut wtr, "Top 20 Occurences:{:?}", top_20)?;
            writeln!(&mut wtr, "The rest of the Units\n{}", hist.format(40, 20))?;
        }
        // Encoding errors
        if !self.encoded_reads.is_empty() {
            let error = self.error_rate();
            writeln!(&mut wtr, "ErrorRate\n{}", error)?;
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
    // let (mut sum, mut sumsq, mut count) = (0f64, 0f64, 0);
    // for x in error_rates {
    //     sum += x;
    //     sumsq += x * x;
    //     count += 1;
    // }
    // let mean = sum / count as f64;
    // let variance = sumsq / count as f64 - mean * mean;
    // assert!(0f64 <= variance);
    // (mean, variance.sqrt())
}
