use definitions::*;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct ErrorRate {
    pub read_error_rate: Vec<f64>,
    pub unit_error_rate: Vec<Vec<f64>>,
    pub median_of_sqrt_err: f64,
}
impl ErrorRate {
    pub fn read(&self, read: u64) -> f64 {
        self.read_error_rate[read as usize]
    }
    pub fn unit(&self, (unit, cluster): (u64, u64)) -> f64 {
        self.unit_error_rate[unit as usize][cluster as usize]
    }
}

type Read = (usize, Vec<(usize, usize, f64)>);
fn residual(errors: &[Read], reads: &[f64], units: &[Vec<f64>]) -> f64 {
    let residual: f64 = errors
        .iter()
        .map(|&(readid, ref errors)| -> f64 {
            let read = reads[readid];
            let data_error: f64 = errors
                .iter()
                .map(|&(unit, cluster, error)| (error - read - units[unit][cluster]).powi(2))
                .sum();
            data_error
        })
        .sum();
    let reg_term: f64 = units.iter().flatten().map(|x| x * x).sum();
    residual + reg_term
}

// Error Rate of the reads, error rate of the units, and the sqrt of the median of the squared error.
pub fn estimate_error_rate(ds: &DataSet, fallback: f64) -> ErrorRate {
    let max_read_id = ds.raw_reads.iter().map(|r| r.id).max().unwrap() as usize;
    let max_unit_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap() as usize;
    let (mut unit_error_rate, unit_counts) = {
        let mut cluster_num = vec![0; max_unit_id + 1];
        for u in ds.selected_chunks.iter() {
            cluster_num[u.id as usize] = u.cluster_num;
        }
        let unit_errors: Vec<_> = cluster_num.iter().map(|&x| vec![0f64; x]).collect();
        let mut counts: Vec<_> = cluster_num.iter().map(|&x| vec![0; x]).collect();
        for read in ds.encoded_reads.iter() {
            for node in read.nodes.iter() {
                let (unit, cluster) = (node.unit as usize, node.cluster as usize);
                counts[unit][cluster] += 1;
            }
        }
        (unit_errors, counts)
    };
    let mut read_error_rate = vec![0f64; max_read_id + 1];
    for read in ds.encoded_reads.iter() {
        read_error_rate[read.id as usize] = fallback;
    }
    let units: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let errors: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|read| {
            let errors: Vec<_> = read
                .nodes
                .iter()
                .map(|node| {
                    let aln = node.recover(units[&node.unit]).1;
                    let errors = aln.iter().filter(|&&x| x != b'|').count();
                    let error = errors as f64 / aln.len() as f64;
                    (node.unit as usize, node.cluster as usize, error)
                })
                .collect::<Vec<_>>();
            (read.id as usize, errors)
        })
        .collect();
    let mut current_resid = residual(&errors, &read_error_rate, &unit_error_rate);
    loop {
        // Re-estimation of unit error rate
        unit_error_rate.iter_mut().flatten().for_each(|x| *x = 0f64);
        for &(readid, ref errors) in errors.iter() {
            let read_error = read_error_rate[readid];
            for &(unit, cluster, error) in errors.iter() {
                unit_error_rate[unit][cluster] += error - read_error;
            }
        }
        for (resid, counts) in unit_error_rate.iter_mut().zip(unit_counts.iter()) {
            for (err, count) in resid.iter_mut().zip(counts) {
                *err = err.max(0f64) / (*count as f64 + 0.1f64);
            }
        }
        // Re-estimate read error rate
        for &(readid, ref errors) in errors.iter() {
            let residual: f64 = errors
                .iter()
                .map(|&(unit, cluster, error)| error - unit_error_rate[unit][cluster])
                .sum();
            read_error_rate[readid] = residual / errors.len() as f64;
        }
        let resid = residual(&errors, &read_error_rate, &unit_error_rate);
        if (current_resid - resid).abs() < 0.00001 {
            break;
        }
        current_resid = resid;
    }
    let mut residuals: Vec<f64> = errors
        .iter()
        .flat_map(|(readid, errors)| {
            let readerror = read_error_rate[*readid as usize];
            errors
                .iter()
                .map(|&(unit, cluster, error)| {
                    let expect = unit_error_rate[unit][cluster] + readerror;
                    error - expect
                })
                .map(|residual| residual.powi(2))
                .collect::<Vec<f64>>()
        })
        .collect();
    let idx = residuals.len() / 2;
    let median = residuals
        .select_nth_unstable_by(idx, |x, y| x.partial_cmp(y).unwrap())
        .1
        .sqrt();
    ErrorRate {
        read_error_rate,
        unit_error_rate,
        median_of_sqrt_err: median,
    }
}
