//! This module includes light statistics metrics, such as error rates and coverages...

use definitions::*;
use rayon::prelude::*;
use std::collections::HashMap;
pub fn error_rate(ds: &DataSet) -> ErrorRate {
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let summaries: Vec<_> = ds
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
    let del_summary = summarize(summaries.iter().map(|x| x.0));
    let ins_summary = summarize(summaries.iter().map(|x| x.1));
    let mism_summary = summarize(summaries.iter().map(|x| x.2));
    let total_summary = summarize(summaries.iter().map(|(x, y, z)| x + y + z));
    ErrorRate::new(del_summary, ins_summary, mism_summary, total_summary)
}

fn summarize<I: std::iter::Iterator<Item = f64>>(error_rates: I) -> (f64, f64) {
    let (mut sum, mut sumsq, mut count) = (0f64, 0f64, 0);
    for x in error_rates {
        sum += x;
        sumsq += x * x;
        count += 1;
    }
    let mean = sum / count as f64;
    let variance = sumsq / count as f64 - mean * mean;
    assert!(0f64 <= variance);
    (mean, variance.sqrt())
}
