//! This module includes light statistics metrics, such as error rates and coverages...

use std::collections::HashMap;

use definitions::*;
pub fn error_rate(ds: &DataSet) -> ErrorRate {
    let (mut dels, mut inss, mut mismatches) = (vec![], vec![], vec![]);
    let ref_units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    for node in ds.encoded_reads.iter().flat_map(|r| r.nodes.iter()) {
        let ref_unit = ref_units[&node.unit];
        let (query, aln, refr) = node.recover(ref_unit);
        let mismat = aln.iter().filter(|&&x| x == b'X').count() as f64;
        let del = query.iter().filter(|&&x| x == b' ').count() as f64;
        let ins = refr.iter().filter(|&&x| x == b' ').count() as f64;
        let aln_len = aln.len() as f64;
        dels.push(del / aln_len);
        inss.push(ins / aln_len);
        mismatches.push(mismat / aln_len);
    }
    let del_summary = summarize(dels.iter().copied());
    let ins_summary = summarize(inss.iter().copied());
    let mism_summary = summarize(mismatches.iter().copied());
    let totals = dels
        .iter()
        .zip(inss.iter())
        .zip(mismatches.iter())
        .map(|((x, y), z)| x + y + z);
    let total_summary = summarize(totals);
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
