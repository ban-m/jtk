use definitions::*;
use std::collections::HashMap;
// ToDo What's this?

#[derive(Debug, Clone)]
pub struct ErrorRate {
    pub read_error_rate: Vec<f64>,
    pub chunk_error_rate: Vec<Vec<f64>>,
    pub median_of_sqrt_err: f64,
}
impl ErrorRate {
    pub fn read(&self, read: u64) -> f64 {
        self.read_error_rate[read as usize]
    }
    pub fn chunk(&self, (chunk, cluster): (u64, u64)) -> f64 {
        self.chunk_error_rate[chunk as usize][cluster as usize]
    }
}

type Read = (usize, Vec<(usize, usize, f64)>);
fn residual(errors: &[Read], reads: &[f64], chunks: &[Vec<f64>]) -> f64 {
    let residual: f64 = errors
        .iter()
        .map(|&(readid, ref errors)| -> f64 {
            let read = reads[readid];
            let data_error: f64 = errors
                .iter()
                .map(|&(chunk, cluster, error)| (error - read - chunks[chunk][cluster]).powi(2))
                .sum();
            data_error
        })
        .sum();
    let reg_term: f64 = chunks.iter().flatten().map(|x| x * x).sum();
    residual + reg_term
}

// Error Rate of the reads, error rate of the chunks, and the sqrt of the median of the squared error.
pub fn estimate_error_rate(ds: &DataSet, fallback: f64) -> ErrorRate {
    let max_read_id = ds.raw_reads.iter().map(|r| r.id).max().unwrap() as usize;
    let max_chunk_id = ds.selected_chunks.iter().map(|c| c.id).max().unwrap() as usize;
    let (mut chunk_error_rate, chunk_counts) = {
        let mut cluster_num = vec![0; max_chunk_id + 1];
        for u in ds.selected_chunks.iter() {
            cluster_num[u.id as usize] = u.cluster_num;
        }
        let chunk_errors: Vec<_> = cluster_num.iter().map(|&x| vec![0f64; x]).collect();
        let mut counts: Vec<_> = cluster_num.iter().map(|&x| vec![0; x]).collect();
        for read in ds.encoded_reads.iter() {
            for node in read.nodes.iter() {
                let (chunk, cluster) = (node.chunk as usize, node.cluster as usize);
                counts[chunk][cluster] += 1;
            }
        }
        (chunk_errors, counts)
    };
    let mut read_error_rate = vec![0f64; max_read_id + 1];
    for read in ds.encoded_reads.iter() {
        read_error_rate[read.id as usize] = fallback;
    }
    let chunks: HashMap<_, _> = ds.selected_chunks.iter().map(|c| (c.id, c)).collect();
    let errors: Vec<_> = ds
        .encoded_reads
        .iter()
        .map(|read| {
            let errors: Vec<_> = read
                .nodes
                .iter()
                .map(|node| {
                    let aln = node.recover(chunks[&node.chunk]).1;
                    let errors = aln.iter().filter(|&&x| x != b'|').count();
                    let error = errors as f64 / aln.len() as f64;
                    (node.chunk as usize, node.cluster as usize, error)
                })
                .collect::<Vec<_>>();
            (read.id as usize, errors)
        })
        .collect();
    let mut current_resid = residual(&errors, &read_error_rate, &chunk_error_rate);
    loop {
        // Re-estimation of a chunk's error rate
        chunk_error_rate
            .iter_mut()
            .flatten()
            .for_each(|x| *x = 0f64);
        for &(readid, ref errors) in errors.iter() {
            let read_error = read_error_rate[readid];
            for &(chunk, cluster, error) in errors.iter() {
                chunk_error_rate[chunk][cluster] += error - read_error;
            }
        }
        for (resid, counts) in chunk_error_rate.iter_mut().zip(chunk_counts.iter()) {
            for (err, count) in resid.iter_mut().zip(counts) {
                *err = err.max(0f64) / (*count as f64 + 0.1f64);
            }
        }
        // Re-estimate read error rate
        for &(readid, ref errors) in errors.iter() {
            let residual: f64 = errors
                .iter()
                .map(|&(chunk, cluster, error)| error - chunk_error_rate[chunk][cluster])
                .sum();
            read_error_rate[readid] = residual / errors.len() as f64;
        }
        let resid = residual(&errors, &read_error_rate, &chunk_error_rate);
        if (current_resid - resid).abs() < 0.00001 {
            break;
        }
        current_resid = resid;
    }
    let mut residuals: Vec<f64> = errors
        .iter()
        .flat_map(|(readid, errors)| {
            let readerror = read_error_rate[*readid];
            errors
                .iter()
                .map(|&(chunk, cluster, error)| {
                    let expect = chunk_error_rate[chunk][cluster] + readerror;
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
        chunk_error_rate,
        median_of_sqrt_err: median,
    }
}
