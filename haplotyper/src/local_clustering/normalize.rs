//! Normalize local clusteirng.
//! Specifically, by calling "normalize_local_clustering", on each chunk,
//! # of reads on clusters would be sorted in descending order with respect to the cluster index.
//! In other words, # of the 0-th cluster would be larger than that of 1-th, # of 1-st would be larger than 2-nd, and so on.

use definitions::*;
pub fn normalize_local_clustering(ds: &mut DataSet) {
    let max_chunk = ds.selected_chunks.iter().map(|x| x.id).max().unwrap() + 1;
    let cov = ds.coverage.unwrap_or(10f64).ceil() as usize;
    let mut pileups: Vec<Vec<&mut _>> = std::iter::repeat_with(|| Vec::with_capacity(cov))
        .take(max_chunk as usize)
        .collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        pileups[node.unit as usize].push(node);
    }
    for pileup in pileups {
        let max_cluster = pileup.iter().map(|n| n.cluster).max().unwrap();
        let mut counts: Vec<(u64, u32)> = (0..=max_cluster).map(|c| (c, 0)).collect();
        for node in pileup.iter() {
            counts[node.cluster as usize].1 += 1;
        }
        counts.sort_by_key(|x| x.1);
        counts.reverse();
        let mut mapsto = vec![0; max_cluster as usize + 1];
        for (to, &(from, _)) in counts.iter().enumerate() {
            mapsto[from as usize] = to as u64;
        }
        for node in pileup {
            node.cluster = mapsto[node.cluster as usize];
        }
    }
}
