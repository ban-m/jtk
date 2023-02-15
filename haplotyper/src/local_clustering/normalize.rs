//! Normalize local clusteirng.
//! Specifically, by calling "normalize_local_clustering", on each chunk,
//! # of reads on clusters would be sorted in descending order with respect to the cluster index.
//! In other words, # of the 0-th cluster would be larger than that of 1-th, # of 1-st would be larger than 2-nd, and so on.
use definitions::*;
pub fn normalize_local_clustering(ds: &mut DataSet) {
    let max_chunk = ds.selected_chunks.iter().map(|x| x.id).max().unwrap() + 1;
    let cov = ds.coverage.unwrap().ceil() as usize;
    let mut pileups: Vec<Vec<&mut _>> = std::iter::repeat_with(|| Vec::with_capacity(cov))
        .take(max_chunk as usize)
        .collect();
    for node in ds.encoded_reads.iter_mut().flat_map(|r| r.nodes.iter_mut()) {
        pileups[node.unit as usize].push(node);
    }
    use std::collections::HashMap;
    let cluster_num: HashMap<_, _> = ds
        .selected_chunks
        .iter()
        .map(|c| (c.id, c.cluster_num))
        .collect();
    for (chunkid, pileup) in pileups
        .into_iter()
        .enumerate()
        .filter(|(_, pu)| !pu.is_empty())
    {
        let max_cluster = cluster_num[&(chunkid as u64)] as u64;
        pileup
            .iter()
            .for_each(|n| assert_eq!(n.posterior.len(), max_cluster as usize));
        let mut counts: Vec<(u64, u32)> = (0..max_cluster).map(|c| (c, 0)).collect();
        for node in pileup.iter() {
            counts[node.cluster as usize].1 += 1;
        }
        counts.sort_by_key(|x| x.1);
        counts.reverse();
        let mut mapsto = vec![0; max_cluster as usize];
        for (to, &(from, _)) in counts.iter().enumerate() {
            mapsto[from as usize] = to as u64;
        }
        let mut indices = mapsto.clone();
        for node in pileup {
            indices
                .iter_mut()
                .zip(mapsto.iter())
                .for_each(|(x, &y)| *x = y);
            node.cluster = mapsto[node.cluster as usize];
            reorder(&mut node.posterior, &mut indices);
            assert!(indices.is_sorted());
        }
    }
}

// Reorder xs according to the indices.
fn reorder<T>(xs: &mut [T], indices: &mut [u64]) {
    let len = xs.len();
    for i in 0..len {
        while indices[i] as usize != i {
            let to = indices[i] as usize;
            xs.swap(i, to);
            indices.swap(i, to);
        }
    }
}

#[cfg(test)]
pub mod tests {
    #[test]
    fn reorder_test() {
        let mut arr = vec![50, 40, 70, 60, 90];
        let mut indices = vec![3, 0, 4, 1, 2];
        super::reorder(&mut arr, &mut indices);
        assert!(indices.is_sorted());
        assert_eq!(arr, vec![40, 60, 90, 50, 70]);
    }
}
